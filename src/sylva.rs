use crate::errors::{MatrixParseError, RuntimeError};
use crate::{errors, utils::*};
use anyhow::*;
use colored::Colorize;
use identity_hash::{IntMap, IntSet};
use itertools::Itertools;
use log::*;
use newick::{Newick, NewickTree};
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;
use syntesuite::genebook::GeneBook;

const RELAXED_SYNTENY_THRESHOLD: OrderedFloat<f32> = OrderedFloat(0.15);
const MIN_INFORMATIVE_SYNTENY: usize = 3;
const CORE_THRESHOLD: usize = 20;

type SpeciesID = usize;
type GeneID = usize;
type PolytomicGeneTree = PTree<usize, SpeciesID>;

use crate::dede::*;
use crate::polytomic_tree::*;

fn round(x: f32, d: i8) -> f32 {
    let ten = 10f32.powi(d as i32);
    (x * ten).round() / ten
}

pub struct Settings {
    pub logs: String,
    pub window: usize,
    pub merge_tandems: bool,
    pub dups_from_sequence: bool,
    pub synteny_threshold: Option<f32>,
}

struct Register<'a> {
    settings: Settings,
    /// A list of the core genes
    core: Vec<GeneID>,
    /// A set  "  "   "    "
    core_set: IntSet<GeneID>,
    /// An ID-to-name mapping of this family
    genes: Vec<String>,
    /// An ID-to-SpeciesID of this family
    species: Vec<SpeciesID>,
    /// An ID-to-lanscape size of this family
    landscape_size: Vec<usize>,
    /// A set of all the species ID found in this family
    all_species: HashSet<SpeciesID>,
    /// The synteny distance matrix of this family
    synteny: VecMatrix<f32>,
    /// The divergence matrix of this family
    divergence: VecMatrix<f32>,
    /// The global species tree
    species_tree: &'a NewickTree,

    extended: Vec<GeneID>,
    solos: Vec<GeneID>,
    fan_out: IntMap<GeneID, Vec<GeneID>>,
    tracked: Vec<String>,
}
impl Register<'_> {
    fn size(&self) -> usize {
        self.genes.len()
    }
}

#[derive(Debug)]
struct Duplication {
    root: SpeciesID,
    content: Vec<GeneID>,
    species: HashSet<SpeciesID>,
}
impl Duplication {
    pub fn pretty(&self, register: &Register) {
        println!(
            "  {:?}",
            view(&register.genes, &self.content).cloned().collect::<Vec<String>>().join(" ")
        );
        println!(
            "  {:?}",
            self.species
                .iter()
                .map(|x| register.species_name(*x))
                .collect::<Vec<String>>()
                .join(" ")
        );
    }
}

type Duplications = Vec<Duplication>;

impl<'st> Register<'st> {
    pub fn species_name(&self, x: SpeciesID) -> String {
        if x > 0 {
            self.species_tree.name(x).map(|x| x.to_string()).unwrap_or_else(|| "UNKNWN".to_string())
        } else {
            "UNKNWN".into()
        }
    }
    pub fn span(&self, mrca: SpeciesID) -> &IntSet<SpeciesID> {
        self.species_tree.cached_leaves_of(mrca)
    }

    pub fn mrca_span(&self, species: impl IntoIterator<Item = usize>) -> &IntSet<SpeciesID> {
        self.span(self.species_tree.mrca(species).unwrap())
    }

    fn rec_elc(&self, me: SpeciesID, missings: &mut IntSet<SpeciesID>, only_large: bool) -> i64 {
        let my_span = self.span(me);
        if self.species_tree[me].is_leaf() {
            0
        } else if my_span.is_subset(missings) {
            if !only_large || my_span.len() > 1 {
                1
            } else {
                0
            }
        } else {
            self.species_tree
                .children(me)
                .unwrap()
                .iter()
                .map(|n| self.rec_elc(*n, missings, only_large))
                .sum()
        }
    }

    pub fn elc<Iter>(&self, species: Iter) -> i64
    where
        Iter: IntoIterator<Item = usize>,
        Iter::IntoIter: Clone,
    {
        let species = species.into_iter();

        let mrca = self.species_tree.mrca(species.clone()).unwrap_or(0);
        self.elc_from(species.clone(), mrca)
    }

    pub fn elc_from(&self, species: impl IntoIterator<Item = usize>, mrca: SpeciesID) -> i64 {
        let species: IntSet<SpeciesID> = species.into_iter().collect();
        if species.is_empty() {
            return 0;
        }

        let mut missings = self.span(mrca).difference(&species).copied().collect::<IntSet<_>>();
        if missings.is_empty() {
            0
        } else {
            self.rec_elc(mrca, &mut missings, false)
        }
    }

    #[allow(dead_code)]
    pub fn ellc<'a>(&self, species: impl IntoIterator<Item = &'a usize>) -> i64 {
        let species: Vec<SpeciesID> = species.into_iter().copied().collect();
        let mrca = self.species_tree.mrca(species.iter().cloned()).unwrap(); // FIXME
        let mut missings = IntSet::<SpeciesID>::from_iter(self.span(mrca).iter().copied())
            .difference(&IntSet::from_iter(species.iter().copied()))
            .copied()
            .collect::<IntSet<_>>();
        if missings.is_empty() {
            0
        } else {
            self.rec_elc(mrca, &mut missings, true)
        }
    }

    pub fn make_label(&self, x: GeneID) -> String {
        format!("{}[&&NHX:S={}]", &self.genes[x], self.species_name(self.species[x]))
    }
}

fn jaccard<T: Eq + Hash>(a: &HashSet<T>, b: &HashSet<T>) -> f32 {
    a.intersection(b).count() as f32 / a.union(b).count() as f32
}

fn parse_dist_matrix<S: AsRef<str>>(filename: &str, ids: &[S]) -> Result<VecMatrix<f32>> {
    let n = ids.len();
    let mut r = VecMatrix::<f32>::with_elem(n, n, 0.);
    let mut tmp = VecMatrix::<f32>::with_elem(n, n, 0.);
    let mut r2g = HashMap::<String, usize>::new();

    let mut lines = BufReader::new(File::open(filename)?).lines();
    let nn = lines.next().ok_or(MatrixParseError::SizeMissing)??.parse::<usize>()?;
    ensure!(
        nn == n,
        format!("{} does not have the expected size (expected {}, found {})", filename, n, nn)
    );
    for (i, l) in lines.enumerate() {
        let l = l?;
        let mut s = l.split('\t');
        let gene = s.next().ok_or(MatrixParseError::ErroneousLine)?;
        r2g.insert(gene.into(), i);
        for (j, x) in s.enumerate() {
            tmp[(i, j)] = x.parse::<f32>()?;
        }
    }

    for i in 0..n {
        let gi = ids[i].as_ref();
        for j in 0..n {
            let gj = ids[j].as_ref();
            r[(i, j)] = tmp[(
                *r2g.get(gi).with_context(|| format!("`{}` not found in {}", gi, filename))?,
                *r2g.get(gj).with_context(|| format!("`{}` not found in {}", gj, filename))?,
            )];
        }
    }

    Ok(r)
}

fn make_register<'a>(
    id: &str,
    genes: &[String],
    book: &GeneBook,
    species_tree: &'a NewickTree,
    syntenies: &str,
    divergences: &str,
    tracked: &[String],
    settings: Settings,
) -> Result<Register<'a>> {
    info!("Building register");

    info!("Parsing synteny matrix");
    let synteny_matrix = &format!("{}/{}.dist", syntenies, id);
    let synteny = parse_dist_matrix(synteny_matrix, genes).map_err(|e| {
        RuntimeError::FailedToReadMatrix { source: e, filename: synteny_matrix.to_owned() }
    })?;

    info!("Parsing divergence matrix");
    let divergence_matrix = &format!("{}/{}.dist", divergences, id);
    let divergence = parse_dist_matrix(divergence_matrix, genes).map_err(|e| {
        RuntimeError::FailedToReadMatrix { source: e, filename: divergence_matrix.to_owned() }
    })?;

    let species2id = species_tree
        .leaves()
        .map(|l| (species_tree.name(l).expect("Found nameless leaf in species tree").to_owned(), l))
        .collect::<HashMap<_, _>>();

    let fan_out = genes
        .iter()
        .map(|g_id| book.get(g_id))
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .enumerate()
        .sorted_by_cached_key(|(_, g)| (g.species.clone(), g.chr.clone(), g.pos))
        .fold(vec![], |mut ax, (i, g)| {
            if settings.merge_tandems {
                if ax.is_empty() {
                    ax.push(vec![i]);
                } else {
                    let same_last =
                        g.left_landscape.last().map(|f| f.family == g.family).unwrap_or(false);
                    if same_last {
                        ax.last_mut().unwrap().push(i);
                    } else {
                        ax.push(vec![i])
                    }
                }
                ax
            } else {
                ax.push(vec![i]);
                ax
            }
        })
        .into_iter()
        .filter_map(|g| if g.len() > 1 { Some((g[0], g[1..].to_vec())) } else { None })
        .collect::<IntMap<_, _>>();
    let secondary_tandems = fan_out.values().flat_map(|g| g.iter()).cloned().collect::<IntSet<_>>();

    let landscapes = genes
        .iter()
        .map(|p| {
            book.get(p)
                .map(|x| x.landscape().collect::<Vec<_>>())
                .map_err(|_| RuntimeError::IdNotFound(p.to_string()).into())
        })
        .collect::<Result<Vec<_>>>()?;
    let landscape_sizes = landscapes.iter().map(Vec::len).collect::<Vec<_>>();

    let mut core = genes
        .iter()
        .enumerate()
        .filter_map(|(i, _)| if landscape_sizes[i] > CORE_THRESHOLD { Some(i) } else { None })
        .collect::<IntSet<_>>();
    if core.is_empty() {
        warn!("Core is empty; using all non-solo instead");
        core = (0..genes.len()).filter(|&p| landscape_sizes[p] > 1).collect();
    }
    let core = core.difference(&secondary_tandems).cloned().collect::<IntSet<_>>();

    let species = genes
        .iter()
        .map(|p| {
            book.get(p).map_err(|_| RuntimeError::IdNotFound(p.to_owned()).into()).and_then(
                |protein| {
                    species2id.get(&protein.species).cloned().ok_or_else(|| {
                        RuntimeError::SpeciesNotFound(protein.species.to_owned()).into()
                    })
                },
            )
        })
        .collect::<Result<Vec<SpeciesID>>>()?;
    let all_species = HashSet::from_iter(species.iter().cloned());

    let extended = genes
        .iter()
        .enumerate()
        .filter_map(|(i, _)| {
            if 1 < landscape_sizes[i] && landscape_sizes[i] <= CORE_THRESHOLD {
                Some(i)
            } else {
                None
            }
        })
        .collect::<IntSet<_>>()
        .difference(&(&core | &secondary_tandems))
        .cloned()
        .collect::<Vec<_>>();
    let solos = Vec::from_iter(
        (0..genes.len())
            .filter(|&p| landscape_sizes[p] == 1)
            .filter(|i| !secondary_tandems.contains(i)),
    );
    let register = Register {
        settings,
        landscape_size: landscape_sizes,
        core: core.iter().cloned().collect_vec(),
        core_set: core,
        genes: genes.to_owned(),
        species,
        all_species,
        synteny,
        divergence,
        species_tree,
        extended,
        solos,
        fan_out,
        tracked: tracked.to_owned(),
    };
    Ok(register)
}

fn mean(x: &[f32]) -> f32 {
    x.iter().sum::<f32>() / x.len() as f32
}

fn find_threshold(register: &Register) -> f32 {
    let all_ts = register
        .synteny
        .masked(&register.core, &register.core)
        .iter()
        .cloned()
        .collect::<Vec<f32>>();
    let mut tts = all_ts
        .iter()
        .cloned()
        .filter(|&x| RELAXED_SYNTENY_THRESHOLD < OrderedFloat(x) && x < 0.7)
        .collect::<Vec<f32>>();
    tts.sort_by(|a, b| a.partial_cmp(b).unwrap());
    tts.dedup();

    if tts.is_empty() {
        let t = mean(&all_ts);
        warn!("No suitable threshold found. Using {}", t);
        return t;
    }
    if tts.len() == 1 {
        warn!("Only one threshold found: using {}", tts[0]);
    }

    let mut avg_comp = Vec::new();
    let mut counts = Vec::<i64>::new();
    let mut reds = Vec::new();

    info!("Pre-computing clusters...");
    let clusterss = tts
        .par_iter()
        .map(|&t| {
            let mut clusters =
                clusterize(&register.synteny.masked(&register.core, &register.core), t)
                    .into_iter()
                    .filter(|c| !c.is_empty())
                    .collect::<Vec<_>>();
            clusters.iter_mut().for_each(|c| c.iter_mut().for_each(|x| *x = register.core[*x]));

            clusters
        })
        .collect::<Vec<_>>();

    info!("Computing optimization targets");
    for (_t, clusters) in tts.iter().zip(clusterss.into_iter()) {
        let redundancies = clusters.iter().map(|c| c.len()).sum::<usize>() as f32
            / clusters
                .iter()
                .map(|c| view(&register.species, c).copied().collect::<HashSet<SpeciesID>>().len())
                .sum::<usize>() as f32;

        let compacities = mean(
            &clusters
                .iter()
                .map(|c| {
                    let my_species =
                        IntSet::<SpeciesID>::from_iter(view(&register.species, c).copied());
                    my_species.len() as f32
                        / register.mrca_span(my_species.iter().cloned()).len() as f32
                })
                .collect::<Vec<f32>>(),
        );

        reds.push(round(redundancies, 2));
        counts.push(clusters.len() as i64);
        avg_comp.push(round(compacities, 2));
    }

    info!("Finding optimal threshold");
    let mut cmpct_maxima = (1..tts.len() - 1)
        .filter(|&i| avg_comp[i - 1] < avg_comp[i] && avg_comp[i] >= avg_comp[i + 1])
        .collect::<Vec<_>>();

    if cmpct_maxima.is_empty() {
        let t = mean(&tts);
        warn!("No compactness peak found in {:?}; using {}", tts, t);
        return t;
    }

    for &i in cmpct_maxima.iter() {
        debug!(
            "THR: {:.3} CNT: {:3} CMPCT: {:.2} RED: {:.2}",
            tts[i], counts[i], avg_comp[i], reds[i]
        );
    }
    debug!("\n");

    let best_red = view(&reds, &cmpct_maxima).fold(f32::INFINITY, |a, &b| a.min(b));
    cmpct_maxima.retain(|i| reds[*i] <= 1.1 * best_red);

    let best_comp = view(&avg_comp, &cmpct_maxima).fold(0f32, |a, &b| a.max(b));
    cmpct_maxima.retain(|i| avg_comp[*i] >= 0.95 * best_comp);

    let i = cmpct_maxima[0];
    debug!("SELECTED:");
    debug!("THR: {:.3} CNT: {:3} CMPCT: {:.2} RED: {:.2}", tts[i], counts[i], avg_comp[i], reds[i]);

    tts[i]
}

fn clusterize<T: PartialOrd>(m: &dyn Matrix<T>, threshold: T) -> Vec<Vec<usize>> {
    assert!(m.is_square());
    let n = m.nrows();
    let mut r = vec![0; n];
    let mut current_clus = 0;

    for i in 0..n {
        let mut my_clus = if r[i] == 0 {
            current_clus += 1;
            r[i] = current_clus;
            current_clus
        } else {
            r[i]
        };

        for j in 0..n {
            if m[(i, j)] >= threshold {
                if r[j] == 0 || r[j] == my_clus {
                    r[j] = my_clus;
                } else {
                    for k in 0..j {
                        if r[k] == my_clus {
                            r[k] = r[j];
                        }
                    }
                    my_clus = r[j];
                }
            }
        }
    }

    let mut clusters = HashMap::<usize, Vec<usize>>::new();
    for (i, r_i) in r.iter().enumerate() {
        clusters.entry(*r_i).or_default().push(i);
    }

    clusters.into_values().collect::<Vec<Vec<usize>>>()
}

fn remove_solos_clusters(t: &mut PolytomicGeneTree) -> Vec<NodeID> {
    let solo_clusters = t[1]
        .children
        .iter()
        .copied()
        .filter(|&n| t[n].children.is_empty() && t[n].content_len() == 1)
        .collect::<Vec<_>>();
    let r = solo_clusters
        .iter()
        .copied()
        .flat_map(|n| t[n].content_iter())
        .copied()
        .collect::<Vec<_>>();
    t.delete_nodes(&solo_clusters);
    r
}

fn inject_extended(t: &mut PolytomicGeneTree, register: &Register, root: usize) {
    let extended_clusters =
        register.extended.iter().fold(Vec::<Vec<GeneID>>::new(), |mut ax, &g| {
            // let target = ax
            //     .iter()
            //     .enumerate()
            //     .find(|(_, gs)| gs.iter().any(|&g0| register.are_satellites(g0, g)))
            //     .map(|x| x.0);

            // if let Some(t) = target {
            //     ax[t].push(g);
            // } else {
            ax.push(vec![g]);
            // };
            ax
        });

    let mut extended = Vec::new();
    for cluster in extended_clusters.into_iter() {
        if cluster.len() == 1 {
            extended.push(cluster[0]);
        } else {
            t.add_node(&cluster, 0, Some(root));
        }
    }

    let mut local_synteny = register.synteny.clone();
    for i in extended.iter() {
        for j in 0..register.size() {
            if j != *i {
                let score = local_synteny[(*i, j)]
                    * register.landscape_size[*i].max(register.landscape_size[j]) as f32;
                let new_score = score / register.landscape_size[*i] as f32;
                local_synteny[(*i, j)] = new_score;
            }
        }
    }
    extended.sort_by_cached_key(|&id| {
        (
            -OrderedFloat(register.synteny.masked(&[id], &register.core).max()),
            OrderedFloat(register.divergence.masked(&[id], &register.core).min()),
        )
    });

    let core_content = t
        .nodes()
        .copied()
        .map(|c| {
            (
                c,
                HashSet::<usize>::from_iter(t[c].content_iter().copied())
                    .intersection(&HashSet::from_iter(register.core.iter().copied()))
                    .copied()
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<HashMap<_, _>>();
    let mut cached_elcs = t
        .nodes()
        .copied()
        .par_bridge()
        .filter(|c| !t[*c].content_is_empty())
        .map(|c| (c, register.elc(clonable_view_cloned(&register.species, t[c].content_slice()))))
        .collect::<HashMap<_, _>>();

    for id in extended.iter() {
        let log = register.tracked.contains(&register.genes[*id]);

        let mut candidate_clusters = t
            .nodes()
            .copied()
            .par_bridge()
            .filter(|c| !t[*c].content_is_empty())
            .map(|c| {
                let c_content = t[c].content_slice();
                let delta_elc = register.elc(
                    clonable_view_cloned(&register.species, c_content)
                        .chain([register.species[*id]].iter().cloned()),
                ) - cached_elcs[&c];
                let divergence = round(register.divergence.masked(&[*id], c_content).min(), 2);
                let synteny = -round(local_synteny.masked(&[*id], &core_content[&c]).max(), 2);
                (c, (delta_elc, OrderedFloat(divergence), OrderedFloat(synteny)))
            })
            .collect::<Vec<_>>();

        let syntenies = candidate_clusters.iter().map(|c| c.1 .2).collect::<Vec<_>>();
        let elcs = candidate_clusters.iter().map(|c| c.1 .0).collect::<Vec<_>>();

        const ELC_THRESHOLD: i64 = 6;
        if syntenies.iter().any(|&s| s <= -RELAXED_SYNTENY_THRESHOLD) {
            // Synteny is negated
            if elcs.iter().any(|&e| e < ELC_THRESHOLD) {
                candidate_clusters.retain(|c| c.1 .0 < ELC_THRESHOLD);
            }
            if register.landscape_size[*id] > MIN_INFORMATIVE_SYNTENY {
                if log {
                    info!("Sorting by synteny/ELC/div");
                }
                candidate_clusters.sort_by_key(|c| (c.1 .2, c.1 .0, c.1 .1));
            } else {
                if log {
                    info!("Sorting by ELC/synteny/div");
                }
                candidate_clusters.sort_by_key(|c| (c.1 .0, c.1 .2, c.1 .1));
            }
        } else {
            if log {
                info!("Sorting by ELC/div/synteny");
            }
            candidate_clusters.sort_by_key(|c| c.1);
        }

        if log {
            info!("EXTENDED {}", register.genes[*id],);
            for cc in &candidate_clusters {
                let c = cc.1;
                info!(
                    "    {:20} --> SYN: {:.2} DV: {:2.2} ΔELC: {}",
                    register.genes[t[cc.0].content_slice()[0]],
                    c.2,
                    c.1,
                    c.0,
                );
            }
        }

        if let Some(cluster) = candidate_clusters.first() {
            let parent = cluster.0;
            t.add_leave(parent, *id);
            t[parent].tag = register
                .species_tree
                .mrca(clonable_view_cloned(&register.species, t[parent].content_slice()))
                .unwrap();
            cached_elcs.insert(
                parent,
                register.elc(clonable_view_cloned(&register.species, t[parent].content_slice())),
            );
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum Metric {
    Synteny(OrderedFloat<f32>),
    Sequence(OrderedFloat<f32>),
}
impl Metric {
    fn furthest(&self) -> OrderedFloat<f32> {
        match self {
            Metric::Synteny(synteny) => -synteny,
            Metric::Sequence(sequence) => *sequence,
        }
    }
}

fn grow_duplication(
    sources: &mut HashMap<SpeciesID, Vec<GeneID>>,
    seed_species: SpeciesID,
    register: &Register,
) -> Duplications {
    let mut ds = sources[&seed_species]
        .iter()
        .map(|&seed| Duplication {
            root: register.species[seed],
            content: vec![seed],
            species: HashSet::from_iter([seed_species]),
        })
        .collect::<Vec<_>>();
    sources.get_mut(&seed_species).unwrap().clear();
    let log = view(&register.genes, ds.iter().flat_map(|d| d.content.iter()))
        .any(|s| register.tracked.contains(s));
    if log {
        debug!("Seed: {}", register.species_name(seed_species));
        ds.iter().for_each(|d| d.pretty(register))
    }
    let mut n = ds[0].root;

    let mut history = Vec::new();
    while register.species_tree.parent(n).is_some() {
        n = register.species_tree.parent(n).unwrap();
        let mut history_entry = vec![];
        let total_span = register.species_tree.leaves_of(n);

        struct Link {
            arm: usize,
            gene: GeneID,
            metric: Metric,
            species: SpeciesID,
        }
        let mut filled = HashMap::<(usize, SpeciesID), bool>::new();
        let mut links = Vec::new();
        let mut putatives = ds.iter().map(|_| Vec::new()).collect::<Vec<_>>();

        for (i, d) in ds.iter().enumerate() {
            let new_span = total_span
                .iter()
                .copied()
                .filter(|x| !d.species.contains(x) && sources.contains_key(x))
                .collect::<Vec<_>>();

            for species in new_span {
                let candidates = sources.entry(species).or_default();
                candidates.sort_by_cached_key(|&c| {
                    OrderedFloat(-register.synteny.masked(&[c], &d.content).max())
                });

                for c in candidates {
                    let metric = if register.settings.dups_from_sequence {
                        Metric::Sequence(OrderedFloat(
                            register.divergence.masked(&[*c], &d.content).min(),
                        ))
                    } else {
                        Metric::Synteny(OrderedFloat(
                            register.synteny.masked(&[*c], &d.content).max(),
                        ))
                    };
                    links.push(Link { arm: i, gene: *c, metric, species });
                }
            }
        }
        links.sort_by_key(|l| {
            (-i64::try_from(ds[l.arm].content.len()).unwrap(), l.metric.furthest())
        });

        for l in links {
            if !filled.get(&(l.arm, l.species)).unwrap_or(&false)
                && sources[&l.species].contains(&l.gene)
            {
                putatives[l.arm].push((l.gene, l.species, l.metric));
                sources.get_mut(&l.species).unwrap().retain(|x| *x != l.gene);
                filled.insert((l.arm, l.species), true);
            }
        }

        for (i, d) in ds.iter_mut().enumerate() {
            for (new, species, _) in putatives[i].iter() {
                d.content.push(*new);
                d.species.insert(*species);
            }
        }

        if log {
            for (i, pp) in putatives.iter().enumerate() {
                for p in pp.iter() {
                    debug!(
                        "Arm #{} gets {}/{} SYN: {:?}",
                        i,
                        register.species_name(p.1),
                        register.genes[p.0],
                        p.2
                    );
                }
            }
        }

        let mut dcss = VecMatrix::<f32>::with_elem(ds.len(), ds.len(), 0.);
        for i in 0..ds.len() {
            for j in i..ds.len() {
                dcss[(i, j)] = ds[i].species.intersection(&ds[i].species).count() as f32
                    / ds[j].species.union(&ds[i].species).count() as f32;
                dcss[(j, i)] = ds[j].species.intersection(&ds[j].species).count() as f32
                    / ds[i].species.union(&ds[j].species).count() as f32;
            }
        }
        if log {
            debug!("DCSS:");
            debug!("{}", &dcss);
        }

        let lengths = ds.iter().map(|d| d.species.len()).collect::<Vec<_>>();
        for (i, d) in ds.iter_mut().enumerate() {
            let putative_dcss = (0..dcss.ncols())
                .filter(|j| *j != i)
                .filter(|j| lengths[*j] > 1)
                .map(|j| dcss[(i, j)])
                .filter(|&x| (0.6..=1.5).contains(&x))
                .map(OrderedFloat)
                .collect::<Vec<_>>();
            let best_dcs = if putative_dcss.is_empty() {
                OrderedFloat(-1.0)
            } else {
                *putative_dcss.iter().max().unwrap()
            };
            let filling = d.species.len() as f32
                / register
                    .mrca_span(d.species.iter().cloned())
                    .iter()
                    .copied()
                    .collect::<HashSet<_>>()
                    .intersection(&register.all_species)
                    .count() as f32;
            let filling_threshold = if total_span.len() <= 4 { 0.5 } else { 0.5 };

            if log {
                debug!(
                    "Arm #{} THR {:.3} FIL {:.3} PUT {:?} BST {:.3}",
                    i, filling_threshold, filling, putative_dcss, best_dcs
                );
            }
            if (filling > filling_threshold)
                && (!putatives[i].is_empty())
                && best_dcs > OrderedFloat(0.)
            {
                d.root = n;
                history_entry.push((i, putatives[i].clone()));
            } else {
                for (new, species, _synteny) in putatives[i].iter() {
                    d.content.retain(|c| c != new);
                    d.species.remove(species);
                    sources.get_mut(species).unwrap().push(*new);
                }
            }
        }
        history.push(history_entry);
    }

    'rewind: for i in (0..history.len()).rev() {
        let l = history[i].len();
        if l == 0 {
        } else if l == 1 {
            let arm = history[i][0].0;
            for (new, species, _synteny) in history[i][0].1.iter() {
                ds[arm].content.retain(|c| c != new);
                ds[arm].species.retain(|s| s != species);
                sources.get_mut(species).unwrap().push(*new);
            }
        } else {
            break 'rewind;
        }
    }

    for d in ds.iter_mut() {
        d.root = register.species_tree.mrca(d.species.iter().cloned()).unwrap();
    }

    ds
}

fn create_duplications(
    register: &Register,
    tree: &PolytomicGeneTree,
    reference: usize,
) -> Vec<Duplications> {
    let mut dups = Vec::new();
    let mut sources = HashMap::<SpeciesID, Vec<GeneID>>::new();
    for i in tree[reference].content_slice() {
        sources.entry(register.species[*i]).or_default().push(*i);
    }

    let mut seed_speciess = sources
        .keys()
        .copied()
        .filter(|s| sources[s].len() > 1)
        .sorted_by_cached_key(|s| {
            (-(sources[s].len() as i64), -register.species_tree.node_topological_depth(*s).unwrap())
        })
        .collect::<Vec<_>>();
    let duplicated_species =
        sources.keys().copied().filter(|s| sources[s].len() > 1).collect::<HashSet<_>>();

    while let Some(seed_species) = seed_speciess.first() {
        let new_family = grow_duplication(&mut sources, *seed_species, register);
        dups.push(new_family);
        sources.remove(seed_species);
        seed_speciess.retain(|s| sources.get(s).map(|v| v.len() > 1).unwrap_or(false));
    }

    assert!(sources.values().all(|l| l.len() <= 1));

    // Remaining singletons from previously duplicated species should be treated independently
    let remaining_species = HashSet::<SpeciesID>::from_iter(
        sources.iter().filter_map(|(k, v)| if !v.is_empty() { Some(k) } else { None }).copied(),
    );
    if !remaining_species.is_empty() {
        for s in remaining_species.intersection(&duplicated_species) {
            dups.push(vec![Duplication {
                root: register.species_tree.mrca([*s]).unwrap(),
                content: sources[s].clone(),
                species: HashSet::from_iter([*s].into_iter()),
            }]);
            sources.remove(s);
        }
    }

    // If the remaining species make for a subset of the existing duplications, they should go under
    let remaining_species = HashSet::from_iter(
        sources.iter().filter_map(|(k, v)| if !v.is_empty() { Some(k) } else { None }).copied(),
    );
    if !remaining_species.is_empty()
        && dups.iter().flat_map(|d| d.iter()).any(|d| {
            register
                .species_tree
                .descendants(d.root)
                .into_iter()
                .collect::<HashSet<_>>()
                .is_superset(&remaining_species)
        })
    {
        let remaining_genes =
            sources.values().filter(|v| !v.is_empty()).map(|v| v[0]).collect::<Vec<_>>();
        dups.push(vec![Duplication {
            root: register.species_tree.mrca(remaining_species.iter().cloned()).unwrap(),
            content: remaining_genes,
            species: remaining_species.into_iter().collect(),
        }]);
    }

    dups
}

fn resolve_duplications(t: &mut PolytomicGeneTree, register: &Register) {
    let roots = t[1].children.clone();

    for &i in &roots {
        let mut dups = create_duplications(register, t, i);
        dups.sort_by_cached_key(|f| {
            let all_species = f.iter().flat_map(|d| d.content.iter().map(|&x| register.species[x]));
            let unique_species = HashSet::<usize>::from_iter(all_species).len() as i64;
            -unique_species
        });

        let cluster_mrca = register
            .species_tree
            .mrca(clonable_view_cloned(&register.species, t[i].content_slice()))
            .unwrap();
        let my_owns = t[i]
            .content_slice()
            .iter()
            .copied()
            .filter(|g| dups.iter().flat_map(|ds| ds.iter()).all(|d| !d.content.contains(g)))
            .collect::<Vec<_>>();
        t.clear_leaves(i);
        reconcile(t, &my_owns, register, Some(i), Some(cluster_mrca), true);

        for f in dups.iter_mut() {
            // We use pop, i.e. implicit reverse
            f.sort_by_cached_key(|d| d.content.iter().unique_by(|g| register.species[**g]).count());

            while let Some(d) = f.pop() {
                let log = view(&register.genes, &d.content).any(|x| register.tracked.contains(x));
                if log {
                    info!(
                        "Current arm: {} {}",
                        register.genes[d.content[0]],
                        register.species_name(d.root)
                    );
                }

                let root_candidates = t
                    .descendants(i)
                    .iter()
                    .copied()
                    .filter(|n| t[*n].tag == d.root)
                    .sorted_by_cached_key(|n| {
                        let leaves = t.full_descendant_leaves(*n);
                        let synteny = if leaves.is_empty() {
                            OrderedFloat(0.0)
                        } else {
                            OrderedFloat(
                                register
                                    .synteny
                                    .masked_from_iter(d.content.iter().cloned(), leaves.into_iter())
                                    .max(),
                            )
                        };
                        let t_depth = t.topo_depth(*n) as i64;
                        (-synteny, -t_depth)
                    })
                    .collect::<Vec<_>>();

                if log {
                    for &c in root_candidates.iter() {
                        println!(
                            "TOPODPTH: {:3} SYN: {:2.2}",
                            t.topo_depth(c),
                            register
                                .synteny
                                .masked_from_iter(
                                    d.content.iter().cloned(),
                                    t.full_descendant_leaves(c).iter().cloned()
                                )
                                .max()
                        );
                    }
                }

                let current_root = if let Some(root) = root_candidates.first() { *root } else { i };

                // If the putative root is a leaf node, make it a non-leaf node
                if !t[current_root].content_is_empty() {
                    let leaves = t[current_root].content_slice().to_vec();
                    t.clear_leaves(current_root);
                    let _ = t.add_node(
                        &leaves,
                        register
                            .species_tree
                            .mrca(leaves.iter().map(|g| &register.species[*g]).cloned())
                            .unwrap(),
                        Some(current_root),
                    );
                }

                // If the putative root has a single child, we can reconcile directly in it
                // Otherwise, we need to create an intermediate node
                if t[current_root].children.len() > 1 {
                    let node_alpha = t.add_node(&[], t[current_root].tag, None);
                    for c in t[current_root].children.clone().into_iter() {
                        t.move_node(c, node_alpha);
                    }
                    t.move_node(node_alpha, current_root);
                    reconcile(
                        t,
                        &d.content,
                        register,
                        Some(current_root),
                        Some(t[current_root].tag),
                        true,
                    );
                } else {
                    reconcile(t, &d.content, register, Some(current_root), None, true);
                }
            }
        }
    }
}

fn expand_meta(t: &mut PolytomicGeneTree, r: &Register) -> Result<()> {
    const MIN_CS: f32 = 0.5;
    let mut todos = r
        .fan_out
        .iter()
        .map(|(k, v)| (t.find_node(|n| n.content_slice().contains(k)).unwrap(), k, v))
        .sorted_by_key(|(a, _, _)| t.topo_depth(*a))
        .map(|(a, k, v)| (a, std::iter::once(k).chain(v.iter()).cloned().collect::<Vec<_>>()))
        .collect::<Vec<(usize, Vec<usize>)>>();
    let species = todos.iter().map(|ts| r.species[ts.1[0]]).collect::<Vec<_>>();
    for todo in &todos {
        trace!(
            "Clearing {:?}",
            t[todo.0].content_iter().map(|g| r.genes[*g].as_str()).collect_vec()
        );
        t[todo.0].content_clear();
    }

    trace!("Tandems found in {:?}", species.iter().map(|s| r.species_name(*s)).collect_vec());

    while let Some(todo) = todos.last().and_then(|ts| ts.1.first().map(|x| (ts.0, *x))) {
        trace!("Opening at {:?}", todo);
        let anchor = todo.0;
        let genes = vec![todo.1];

        let mut current_best: (usize, Vec<usize>) = (anchor, genes);
        let mut current_mrca = anchor;

        while let Some(new_mrca) = t[current_mrca].parent {
            let all_descendants = t.full_descendants(new_mrca);
            let new_genes: Vec<GeneID> = todos
                .iter()
                .filter_map(
                    |(anchor, gs)| {
                        if all_descendants.contains(anchor) {
                            Some(gs[0])
                        } else {
                            None
                        }
                    },
                )
                .collect_vec();
            let cs = new_genes.iter().map(|g| r.species[*g]).collect::<IntSet<_>>().len() as f32
                / r.span(t[new_mrca].tag).len() as f32;
            if cs >= MIN_CS && new_genes.len() > current_best.1.len() {
                current_best = (new_mrca, new_genes)
            }
            trace!("{}: cs = {}", r.species_name(t[new_mrca].tag), cs);
            current_mrca = new_mrca;
        }
        trace!("Closing.");
        trace!(
            "Extracting {:?} from {:?}@{}",
            current_best.1.iter().map(|g| r.genes[*g].as_str()).collect::<Vec<_>>(),
            current_best.1.iter().map(|g| r.species_name(r.species[*g])).collect::<HashSet<_>>(),
            r.species_name(t[current_best.0].tag),
        );
        let alpha = t.add_node(&[], t[current_best.0].tag, Some(t[current_best.0].parent.unwrap()));
        t.move_node(current_best.0, alpha);
        let _ = reconcile(t, &current_best.1, r, Some(alpha), None, false);
        todos.iter_mut().for_each(|ts| ts.1.retain(|g| !current_best.1.contains(g)));

        todos.retain(|t| !t.1.is_empty());
    }

    Ok(())
}

/// Represents the crierion used to sort grafting candidates
struct GraftCriteria<F> {
    elc: i64,
    synteny: OrderedFloat<F>,
    divergence: OrderedFloat<F>,
}

/// A set of metadata relative to a putative anchor point for grafting a sub-tree
struct CandidateParent {
    // the node id of this candidate in the tree
    node: usize,
    // the node in the tree that will be used to compute th grafting criterion
    context_node: usize,
    // the species covered by the candidate
    species: IntSet<SpeciesID>,
    // the genes descending from the candidate
    leaves: Vec<usize>,
}

fn make_final_tree(t: &mut PolytomicGeneTree, register: &Register) -> usize {
    fn show_graft<F: ordered_float::Float + std::fmt::Display>(
        parent: &CandidateParent,
        criterion: &GraftCriteria<F>,
        t: &PolytomicGeneTree,
        register: &Register,
    ) {
        info!(
            "{:4}/{:4} ({}) {:4} leaves DELC: {:2} SYN: {:2.2} DV: {:2.2} T: {} {}",
            parent.node,
            parent.context_node,
            register.species_name(register.species[t[parent.context_node].tag]),
            parent.leaves.len(),
            criterion.elc,
            criterion.synteny,
            criterion.divergence,
            register.species_name(t[parent.node].tag),
            parent.leaves.first().map(|g| register.genes[*g].as_str()).unwrap_or("EMPTY")
        );
    }

    let to_keep = t[1].children.iter().map(|&k| t[k].tag).collect::<IntSet<_>>();
    debug!(
        "Partially pruning tree, keeping {:?}",
        to_keep.iter().map(|&s| register.species_name(s)).collect_vec()
    );
    inner_prune_tree(t, 1, &to_keep);
    let mut todo = t[1]
        .children
        .iter()
        .copied()
        .sorted_by_cached_key(|x| {
            (
                -register.species_tree.node_topological_depth(t[*x].tag).unwrap(),
                (t.descendant_leaves(*x).len() as i64),
            )
        }) // .pop(), so implicit reverse ordering
        .collect::<Vec<_>>();

    let aa = todo.pop().unwrap();
    let new_root = reconcile_upstream(t, aa, register, None);

    while let Some(to_graft) = todo.pop() {
        let candidate_parents_id = t
            .descendants(new_root)
            .iter()
            .cloned()
            .filter(|&b| b != to_graft && t[b].tag == t[to_graft].tag)
            .filter(|&b| !t.descendants(to_graft).contains(&b))
            .collect::<Vec<_>>();

        let candidate_parents: Vec<CandidateParent> = candidate_parents_id
            .iter()
            .map(|i| {
                // the context node is chosen to be the oldest parent of the
                // grafting node such that it encompasses no other grafting
                // node.
                let mut context_node = *i;
                while let Some(maybe_context_node) = t[context_node].parent {
                    if !candidate_parents_id
                        .iter()
                        .filter(|c| *c != i)
                        .any(|c| t.descendants(maybe_context_node).contains(c))
                    {
                        context_node = maybe_context_node;
                    } else {
                        break;
                    }
                }

                CandidateParent {
                    node: *i,
                    context_node,
                    species: t
                        .descendant_leaves(context_node)
                        .iter()
                        .map(|g| register.species[*g])
                        .collect::<IntSet<_>>(),
                    leaves: t
                        .descendant_leaves(context_node)
                        .intersection(&register.core_set)
                        .cloned()
                        .collect::<Vec<_>>(),
                }
            })
            .collect();

        let mut leaves = HashMap::new();
        // NOTE: do not intersect with the core set here, as grafted tree may
        // very well not contain any member of it.
        leaves.insert(to_graft, t.descendant_leaves(to_graft).into_iter().cloned().collect_vec());
        let mut speciess = HashMap::new();
        speciess.insert(
            to_graft,
            view_cloned(&register.species, t.descendant_leaves(to_graft))
                .collect::<IntSet<SpeciesID>>(),
        );

        let log = view(&register.genes, t.descendant_leaves(to_graft))
            .any(|s| register.tracked.contains(s));
        if log {
            println!("Injecting {}", register.species_name(t[to_graft].tag));
        }

        let mut candidate_parents = candidate_parents
            .par_iter()
            .map(|b| {
                let mrca = register
                    .species_tree
                    .mrca(speciess[&to_graft].iter().chain(b.species.iter()).cloned())
                    .unwrap();
                let elc = if speciess[&to_graft].is_disjoint(&b.species) {
                    register.elc_from(
                        b.species.iter().cloned().chain(speciess[&to_graft].iter().cloned()),
                        mrca,
                    ) - register.elc_from(b.species.iter().cloned(), mrca)
                } else {
                    register.elc_from(b.species.iter().cloned(), mrca)
                        + register.elc_from(speciess[&to_graft].iter().cloned(), mrca)
                };

                let synteny = leaves[&to_graft]
                    .iter()
                    .map(|l| register.synteny.masked(&[*l], &b.leaves).max())
                    .sum::<f32>()
                    / leaves[&to_graft].len() as f32;
                let divergence = register.divergence.masked(&leaves[&to_graft], &b.leaves).min();

                (
                    b,
                    GraftCriteria {
                        elc,
                        synteny: OrderedFloat(round(synteny, 2)),
                        divergence: OrderedFloat(round(divergence, 2)),
                    },
                )
            })
            .collect::<Vec<_>>();
        candidate_parents.sort_by_key(|c| (c.1.elc, -c.1.synteny, c.1.divergence));

        if log {
            println!("\n\n=== BEFORE ===");
            for b in &candidate_parents {
                let parent = b.0;
                let criterion = &b.1;
                show_graft(parent, criterion, t, register);
            }
        }

        let mut method = "ELC";

        if candidate_parents.iter().map(|c| c.1.synteny).any(|s| s >= RELAXED_SYNTENY_THRESHOLD) {
            // debug!("Filtering by synteny");
            candidate_parents.retain(|c| c.1.synteny >= RELAXED_SYNTENY_THRESHOLD);
            candidate_parents.sort_by_key(|c| {
                (
                    -c.1.synteny,   // Synteny
                    c.1.elc,        // ELC
                    c.1.divergence, // Sequence
                )
            });
            method = "SYN";
        } else if candidate_parents.iter().map(|c| c.1.divergence).any(|d| d <= OrderedFloat(0.5)) {
            // debug!("Filtering by divergence");
            candidate_parents.retain(|c| c.1.divergence <= OrderedFloat(0.5));
            candidate_parents.sort_by_key(|c| {
                (
                    c.1.divergence, // Sequence
                    c.1.elc,        // ELC
                    -c.1.synteny,   // Synteny
                )
            });
            method = "SEQ";
        }

        if log {
            println!("=== AFTER ===");
            for b in &candidate_parents {
                let parent = b.0;
                let criterion = &b.1;
                show_graft(parent, criterion, t, register);
            }
        }
        if let Some(optimum) = candidate_parents.first() {
            let parent = optimum.0;
            let child = to_graft;
            t[child].metadata.insert("METHOD".into(), method.to_owned());

            if log {
                println!(
                    "Parent: {}/{}",
                    t[parent.node].content_len(),
                    t[parent.node].children.len()
                );
            }

            // If the child is plugged on a completely empty subtree,
            // it should replace it to lock its docking points
            if t.descendant_leaves(parent.node).is_empty() {
                t.move_node(child, t[parent.node].parent.unwrap());
                t.delete_node(parent.node);
            } else {
                if !t[parent.node].content_is_empty() {
                    let leaves = t[parent.node].content_slice().to_vec();
                    t.clear_leaves(parent.node);
                    let _content = t.add_node(
                        &leaves,
                        register
                            .species_tree
                            .mrca(clonable_view_cloned(&register.species, &leaves))
                            .unwrap(),
                        Some(parent.node),
                    );
                }

                if t[parent.node].children.len() > 1 {
                    let node_alpha = t.add_node(&[], t[parent.node].tag, None);
                    for c in t[parent.node].children.clone().into_iter() {
                        t.move_node(c, node_alpha);
                    }
                    t.move_node(node_alpha, parent.node);
                }

                t.move_node(child, parent.node);
            }
        }
    }

    let mut tops = t[1]
        .children
        .iter()
        .copied()
        .sorted_by_cached_key(|c| {
            t[1].children
                .iter()
                .filter(|o| *o != c)
                .map(|o| {
                    OrderedFloat(-jaccard(
                        &view(&register.species, t.descendant_leaves(*c))
                            .cloned()
                            .collect::<HashSet<_>>(),
                        &view(&register.species, t.descendant_leaves(*o))
                            .cloned()
                            .collect::<HashSet<_>>(),
                    ))
                })
                .max()
        })
        .collect::<Vec<_>>();
    let mut root = new_root;
    while let Some(to_plug) = tops.pop() {
        if !tops.is_empty() {
            let node_alpha = t.add_node(&[], t[to_plug].tag, Some(root));
            t.move_node(to_plug, node_alpha);
            root = t.add_node(&[], t[to_plug].tag, Some(node_alpha));
            t[node_alpha].tag = register
                .species_tree
                .mrca(view_cloned(&register.species, t.descendant_leaves(node_alpha)))
                .unwrap();
        } else {
            t.move_node(to_plug, root);
        }

        if !t.descendant_leaves(root).is_empty() {
            t[root].tag = register
                .species_tree
                .mrca(view_cloned(&register.species, t.descendant_leaves(root)))
                .unwrap();
        }
    }

    assert!(t.full_descendant_leaves(1).is_empty());
    assert!(t[1].children.is_empty());
    assert!(t[1].content_is_empty());
    t.delete_node(1);

    new_root
}

fn inner_prune_tree(tree: &mut PolytomicGeneTree, root: usize, tag_to_keep: &IntSet<SpeciesID>) {
    tree.disable_leave_cache();
    info!("Pruning empty leaf nodes -- from {}", tree.nodes().count());
    loop {
        let todos = tree
            .nodes()
            .copied()
            .filter(|&k| tree[k].children.is_empty() && tree[k].content_is_empty())
            .filter(|&k| !tag_to_keep.contains(&tree[k].tag))
            .collect::<Vec<_>>();
        if todos.is_empty() {
            break;
        } else {
            tree.delete_nodes(&todos);
        }
    }
    info!("...to {}", tree.nodes().count());

    info!("Pruning scaffolding -- from {}", tree.nodes().count());
    loop {
        let todo = tree
            .nodes()
            .copied()
            .filter(|&k| !tag_to_keep.contains(&tree[k].tag))
            .find(|&k| k != root && tree[k].children.len() == 1 && tree[k].content_is_empty());

        if let Some(k) = todo {
            let parent = tree[k].parent.unwrap();

            tree.move_node(tree[k].children[0], parent);
            tree.delete_node(k);
        } else {
            break;
        }
    }
    info!("...to {}", tree.nodes().count());
    tree.enable_leave_cache();
}

fn prune_tree(tree: &mut PolytomicGeneTree, root: usize) {
    info!("Pruning empty leaf nodes -- from {}", tree.nodes().count());
    loop {
        let todos = tree
            .nodes()
            .copied()
            .filter(|&k| tree[k].children.is_empty() && tree[k].content_is_empty())
            .collect::<Vec<_>>();
        if todos.is_empty() {
            break;
        } else {
            tree.delete_nodes(&todos);
        }
    }
    info!("...to {}", tree.nodes().count());

    info!("Pruning transitions -- from {}", tree.nodes().count());
    loop {
        let todo = tree
            .nodes()
            .copied()
            .find(|&k| tree[k].children.is_empty() && tree[k].content_len() == 1 && k != root);

        if let Some(k) = todo {
            let content = tree[k].content_slice()[0];
            let parent = tree[k].parent.unwrap();
            tree.add_leave(parent, content);
            tree.delete_node(k);
        } else {
            break;
        }
    }
    info!("...to {}", tree.nodes().count());

    info!("Pruning scaffolding -- from {}", tree.nodes().count());
    loop {
        let todo = tree
            .nodes()
            .copied()
            .find(|&k| k != root && tree[k].children.len() == 1 && tree[k].content_is_empty());

        if let Some(k) = todo {
            let parent = tree[k].parent.unwrap();

            tree.move_node(tree[k].children[0], parent);
            tree.delete_node(k);
        } else {
            break;
        }
    }
    info!("...to {}", tree.nodes().count());
}

/// Given a list of genes and the register, reconcile them with the species tree
/// and graft them on the tree or on a new node.
///
/// # Arguments
///
/// * `t` - The gene tree to graft on
/// * `ps` - The IDs of the genes to reconcile
/// * `register` - The global register
/// * `graft` - If defined, where to graft the reconciliation on the gene tree;
///             otherwise the new subtree root node ID will be returned
/// * `from` - If defined, the species to start the reconciliation from; otherwise
///            the MRCA of `ps` will be used
/// * `full` - if true, (empty) nodes will be created for all descendants of `from`
fn reconcile(
    t: &mut PolytomicGeneTree,
    ps: &[GeneID],
    register: &Register,
    graft: Option<usize>,
    from: Option<SpeciesID>,
    full: bool,
) -> usize {
    trace!(
        "reconciling {:?}",
        ps.iter()
            .map(|g| register.species_name(register.species[*g]).to_owned())
            .collect::<Vec<_>>()
    );
    fn rec_add(
        t: &mut PolytomicGeneTree,
        parent: usize,
        phylo_node: SpeciesID,
        ps: &[GeneID],
        speciess: &HashSet<SpeciesID>,
        register: &Register,
        full: bool,
    ) {
        if register.species_tree[phylo_node].is_leaf() {
            let mut current_root = parent;
            let mut to_plug =
                ps.iter().filter(|p| register.species[**p] == phylo_node).collect::<Vec<_>>();
            if !to_plug.is_empty() {
                if true {
                    let mut current_leaf =
                        t.add_node(&[*to_plug.pop().unwrap()], phylo_node, Some(current_root));
                    for &new in to_plug.into_iter() {
                        current_leaf = t.add_node(&[new], phylo_node, Some(current_leaf));
                    }
                } else {
                    for &new in to_plug.into_iter() {
                        current_root = t.add_node(&[new], phylo_node, Some(current_root));
                    }
                }
            } else if full {
                t.add_node(&[], phylo_node, Some(current_root));
            }
        } else {
            let actual_children = register
                .species_tree
                .children(phylo_node)
                .unwrap()
                .iter()
                .filter(|child| {
                    full || !speciess.is_disjoint(&HashSet::from_iter(
                        register.species_tree.leaves_of(**child).iter().copied(),
                    ))
                })
                .collect::<Vec<_>>();
            if !actual_children.is_empty() {
                let me = if t[parent].tag == phylo_node
                    && (t[parent].children.len() + actual_children.len()) <= 2
                {
                    parent
                } else {
                    t.add_node(&[], phylo_node, Some(parent))
                };

                for child in actual_children.into_iter() {
                    rec_add(t, me, *child, ps, speciess, register, full);
                }
            }
        }
    }

    let speciess = view(&register.species, ps).copied().collect::<HashSet<_>>();
    let mrca =
        from.unwrap_or_else(|| register.species_tree.mrca(speciess.iter().cloned()).unwrap());
    let root = graft.unwrap_or_else(|| t.add_node(&[], mrca, None));

    if ps.is_empty() && !full {
        warn!("Empty reconciling");
    } else {
        rec_add(t, root, mrca, ps, &speciess, register, full);
    }

    root
}

fn reconcile_upstream(
    t: &mut PolytomicGeneTree,
    existing: NodeID,
    register: &Register,
    graft: Option<usize>,
) -> usize {
    fn rec_add(
        t: &mut PolytomicGeneTree,
        existing: NodeID,
        parent: usize,
        phylo_node: SpeciesID,
        register: &Register,
    ) {
        if register.species_tree[phylo_node].is_leaf() {
            if phylo_node == t[existing].tag {
                t.move_node(existing, parent);
            } else {
                t.add_node(&[], phylo_node, Some(parent));
            }
        } else {
            let me = t.add_node(&[], phylo_node, Some(parent));
            assert!(me != existing);
            for child in register.species_tree.children(phylo_node).unwrap().iter() {
                if *child == t[existing].tag {
                    t.move_node(existing, me);
                } else {
                    rec_add(t, existing, me, *child, register);
                }
            }
        }
    }

    let mrca = register.species_tree.root();
    let root = graft.unwrap_or_else(|| t.add_node(&[], mrca, None));

    rec_add(t, existing, root, mrca, register);

    root
}

fn do_family(id: &str, register: &Register, logs_root: &str) -> Result<PolytomicGeneTree> {
    info!("Setting synteny threshold...");
    let tt = register.settings.synteny_threshold.unwrap_or_else(|| find_threshold(register) - 0.1);
    info!("Using {}", tt);

    let mut tree = PolytomicGeneTree::new();
    let root = tree.add_node(&[], 0, None);
    let clusters = clusterize(&register.synteny.masked(&register.core, &register.core), tt);
    let satellites = if clusters.is_empty() {
        // In the edge case where the clustering failed, use the extended as single-element clusters
        // This typically happens in very small trees, or tandem-only trees
        for c in register.extended.iter() {
            tree.add_node(&[*c], 0, Some(root));
        }

        File::create(&format!("{}/{}_vanilla.nwk", &logs_root, id))?.write_all(
            tree.to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t)).as_bytes(),
        )?;

        vec![]
    } else {
        for c in clusters.into_iter() {
            tree.add_node(&c.iter().map(|x| register.core[*x]).collect::<Vec<_>>(), 0, Some(root));
        }
        File::create(&format!("{}/{}_vanilla.nwk", &logs_root, id))?.write_all(
            tree.to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t)).as_bytes(),
        )?;

        info!(
            "Injecting {} extended in {} clusters",
            register.extended.len(),
            tree[1].children.len()
        );
        inject_extended(&mut tree, register, root);
        remove_solos_clusters(&mut tree)
    };

    info!("Tagging tree...");
    let to_remove = tree
        .nodes()
        .copied()
        .filter(|&k| k != root && tree[k].content_is_empty() && tree[k].children.is_empty())
        .collect::<Vec<_>>();
    tree.delete_nodes(&to_remove);
    for k in tree[root].children.clone().into_iter() {
        tree[k].tag = register
            .species_tree
            .mrca(clonable_view_cloned(&register.species, tree[k].content_slice()))
            .unwrap();
    }
    info!("Done.");

    File::create(&format!("{}/{}_clusters.nwk", &logs_root, id))?.write_all(
        tree.to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t)).as_bytes(),
    )?;

    info!("Resolving duplications...");
    resolve_duplications(&mut tree, register);
    for s in satellites.iter().chain(register.solos.iter()) {
        tree.add_node(&[*s], register.species[*s], Some(root));
    }
    info!("Done.");
    File::create(&format!("{}/{}_prototree.nwk", &logs_root, id))?.write_all(
        tree.to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t)).as_bytes(),
    )?;

    info!("Assembling final tree -- {} subtrees", tree[root].children.len());
    tree.enable_leave_cache();
    let root = make_final_tree(&mut tree, register);
    tree.disable_leave_cache();
    File::create(&format!("{}/{}_reconciled.nwk", &logs_root, id))?.write_all(
        tree.to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t)).as_bytes(),
    )?;
    info!("Done.");

    expand_meta(&mut tree, register)?;

    info!("Pruning tree...");
    prune_tree(&mut tree, root);
    info!("Done.");

    for k in tree.nodes() {
        if tree[*k].children.len() > 2 {
            warn!("Final tree is polytomic at {}", register.species_name(tree[*k].tag));
        }
    }

    File::create(&format!("{}/{}_synteny.nwk", &logs_root, id))?.write_all(
        tree.to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t)).as_bytes(),
    )?;

    if register.size() != tree.full_descendant_leaves(root).len() {
        error!(
            "Gene mismatch: expected {}, found {}",
            register.size(),
            tree.full_descendant_leaves(root).len(),
        );
        let mut mine = tree.full_descendant_leaves(root).into_iter().collect::<Vec<_>>();
        mine.sort();
        for gs in mine.windows(2) {
            if gs[0] == gs[1] {
                warn!("Double: {}", register.genes[gs[0]]);
            }
        }
        let mine = HashSet::<_>::from_iter(
            tree.full_descendant_leaves(root).iter().map(|l| register.genes[*l].to_owned()),
        );
        let others = HashSet::from_iter(register.genes.iter().cloned());
        let missings = others.difference(&mine).collect::<Vec<_>>();
        info!("{:?}", missings);
        panic!()
    }

    Ok(tree)
}

pub fn do_file(
    filename: &str,
    speciestree_file: &str,
    db_file: &str,
    syntenies: &str,
    divergences: &str,
    settings: Settings,
    timings: &mut Option<File>,
    tracked: &[String],
) -> Result<String> {
    let mut species_tree = newick::one_from_filename(speciestree_file)
        .with_context(|| anyhow!("while opening {}", speciestree_file.yellow().bold()))?;
    species_tree.cache_leaves();

    let gene_families =
        read_genefile(filename).with_context(|| format!("while parsing {}", filename))?;
    if gene_families.is_empty() || gene_families.len() > 1 {
        bail!("{} should contain a single family", filename)
    }
    let family = &gene_families[0];
    let id = &std::path::Path::new(filename)
        .file_stem()
        .with_context(|| "asdfasdf")?
        .to_str()
        .ok_or_else(|| errors::FileError::InvalidFilename(format!("{:?}", filename)))?;

    info!("===== Family {} -- {} proteins =====", id, family.len());
    let logs_root = format!("{}/{}/", &settings.logs, id);
    std::fs::create_dir_all(&logs_root)?;
    let now = Instant::now();
    let register = make_register(
        id,
        family,
        &GeneBook::cached(db_file, settings.window, "id", family)?,
        &species_tree,
        syntenies,
        divergences,
        tracked,
        settings,
    )?;
    let out_tree = do_family(id, &register, &logs_root)?;
    info!("Done in {:.2}s.", now.elapsed().as_secs_f32());
    if let Some(ref mut timings) = timings {
        timings.write_all(
            format!("{},{},{}\n", id, family.len(), now.elapsed().as_secs_f32()).as_bytes(),
        )?;
    }

    let mut x = newick::from_string(
        out_tree.to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t)),
    )?;
    let out_tree = x.get_mut(0).unwrap();
    chainsaw::annotate_mrcas(out_tree, &species_tree)?;
    chainsaw::annotate_duplications(out_tree, &species_tree, false);

    Ok(Newick::to_newick(out_tree, true))
}
