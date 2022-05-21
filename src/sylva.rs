use crate::utils::*;
use anyhow::*;
use itertools::Itertools;
use log::*;
use newick::*;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::time::Instant;

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

fn sim_matrix(
    t: &PolytomicGeneTree,
    register: &Register,
    ks: &[usize],
    f: &dyn Fn(usize, usize) -> f32,
) {
    for (i, k) in ks.iter().enumerate() {
        print!(
            "{:2}/({:3}) {:20} |  ",
            k,
            t[*k].content.len(),
            register.proteins[t[*k].content[0]]
        );
        for l in ks.iter().skip(i) {
            print!("{:2.2}  ", f(*k, *l))
        }
        println!();
    }
}

struct Register<'a> {
    size: usize,
    landscape_size: Vec<usize>,
    core: Vec<GeneID>,
    proteins: Vec<String>,
    genes: Vec<String>,
    species: Vec<SpeciesID>,
    all_species: HashSet<SpeciesID>,
    synteny: VecMatrix<f32>,
    divergence: VecMatrix<f32>,
    species_tree: &'a NewickTree,
    extended: Vec<GeneID>,
    solos: Vec<GeneID>,
}

#[derive(Debug)]
struct Duplication {
    root: SpeciesID,
    content: Vec<GeneID>,
    species: HashSet<SpeciesID>,
}
impl Duplication {
    pub fn pretty<'a>(&self, register: &'a Register) {
        println!(
            "  {:?}",
            view(&register.proteins, &self.content)
                .cloned()
                .collect::<Vec<String>>()
                .join(" ")
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
            self.species_tree[x]
                .data
                .name
                .as_ref()
                .map(|x| x.to_string())
                .unwrap_or("UNKNWN".to_string())
        } else {
            "UNKNWN".into()
        }
    }
    pub fn span(&self, mrca: SpeciesID) -> &[SpeciesID] {
        self.species_tree.cached_leaves_of(mrca)
    }

    pub fn mrca_span<'a>(&self, species: impl IntoIterator<Item = &'a usize>) -> &[SpeciesID] {
        self.span(self.species_tree.mrca(species).unwrap())
    }

    fn rec_elc(&self, me: SpeciesID, missings: &mut HashSet<SpeciesID>, only_large: bool) -> i64 {
        let my_span = HashSet::<SpeciesID>::from_iter(self.span(me).into_iter().copied());
        return if my_span.is_subset(&missings) {
            if !only_large || my_span.len() > 1 {
                1
            } else {
                0
            }
        } else {
            if self.species_tree[me].is_leaf() {
                0
            } else {
                self.species_tree
                    .children(me)
                    .iter()
                    .map(|n| self.rec_elc(*n, missings, only_large))
                    .sum()
            }
        };
    }

    pub fn elc<'a>(&self, species: impl IntoIterator<Item = &'a usize>) -> i64 {
        let species: Vec<SpeciesID> = species.into_iter().copied().collect();
        if species.is_empty() {
            return 0;
        }

        let mrca = self.species_tree.mrca(&species).expect(&format!(
            "No MRCA for {}",
            species
                .iter()
                .map(|&g| self.species_name(g))
                .collect::<Vec<_>>()
                .join(" ")
        ));
        let mut missings = self
            .span(mrca)
            .iter()
            .copied()
            .collect::<HashSet<SpeciesID>>()
            .difference(&HashSet::from_iter(species.iter().copied()))
            .copied()
            .collect::<HashSet<_>>();
        if missings.is_empty() {
            0
        } else {
            self.rec_elc(mrca, &mut missings, false)
        }
    }

    pub fn ellc<'a>(&self, species: impl IntoIterator<Item = &'a usize>) -> i64 {
        let species: Vec<SpeciesID> = species.into_iter().copied().collect();
        let mrca = self.species_tree.mrca(&species).unwrap();
        let mut missings = HashSet::<SpeciesID>::from_iter(self.span(mrca).iter().copied())
            .difference(&HashSet::from_iter(species.iter().copied()))
            .copied()
            .collect::<HashSet<_>>();
        if missings.is_empty() {
            0
        } else {
            self.rec_elc(mrca, &mut missings, true)
        }
    }

    pub fn make_label(&self, x: GeneID) -> String {
        format!(
            "{}[&&NHX:S={}]",
            &self.proteins[x],
            self.species_name(self.species[x])
        )
    }
}

fn jaccard<T: Eq + Hash>(a: &HashSet<T>, b: &HashSet<T>) -> f32 {
    a.intersection(b).count() as f32 / a.union(b).count() as f32
}

fn vec_inter<T: Copy + Eq + Hash>(a: &[T], b: &[T]) -> HashSet<T> {
    HashSet::<T>::from_iter(a.iter().copied())
        .intersection(&HashSet::<T>::from_iter(b.iter().copied()))
        .copied()
        .collect::<HashSet<T>>()
}

fn parse_dist_matrix<S: AsRef<str>>(filename: &str, ids: &[S]) -> Result<VecMatrix<f32>> {
    let n = ids.len();
    let mut r = VecMatrix::<f32>::with_elem(n, n, 0.);
    let mut tmp = VecMatrix::<f32>::with_elem(n, n, 0.);
    let mut r2g = HashMap::<String, usize>::new();

    let mut lines = BufReader::new(File::open(filename)?).lines();
    let nn = lines
        .next()
        .ok_or(anyhow!("No size found in distmat"))??
        .parse::<usize>()?;
    ensure!(
        nn == n,
        format!(
            "{} does not have the expected size (expected {}, found {})",
            filename, n, nn
        )
    );
    for (i, l) in lines.enumerate() {
        let l = l?;
        let mut s = l.split('\t');
        let gene = s.next().ok_or(anyhow!("fsdf"))?.to_owned();
        r2g.insert(gene, i);
        for (j, x) in s.enumerate() {
            tmp[(i, j)] = x.parse::<f32>()?;
        }
    }

    for i in 0..n {
        let gi = ids[i].as_ref();
        for j in 0..n {
            let gj = ids[j].as_ref();
            r[(i, j)] = tmp[(
                *r2g.get(gi)
                    .with_context(|| format!("`{}` not found in {}", gi, filename))?,
                *r2g.get(gj)
                    .with_context(|| format!("`{}` not found in {}", gj, filename))?,
            )];
        }
    }

    Ok(r)
}

fn underscore2point(s: &str) -> String {
    let s = s.replace("_", ".");
    s[0..1].to_uppercase() + &s[1..]
}
fn make_register<'a>(
    id: &str,
    proteins: &[String],
    book: &GeneBook,
    species_tree: &'a NewickTree,
    syntenies: &str,
    divergences: &str,
) -> Result<Register<'a>> {
    info!("Building register");
    let genes = proteins
        .iter()
        .map(|p| {
            book.get(p)
                .map(|p| p.gene.to_owned())
                .with_context(|| format!("while looking up {}", p))
        })
        .collect::<Result<Vec<_>>>()?;

    info!("Parsing synteny matrix");
    let synteny_matrix = &format!("{}/{}.dist", syntenies, id);
    let synteny = parse_dist_matrix(&synteny_matrix, &proteins)
        .with_context(|| format!("while reading synteny matrix `{}`", synteny_matrix))?;

    info!("Parsing divergence matrix");
    let divergence_matrix = &format!("{}/{}.dist", divergences, id);
    let divergence = parse_dist_matrix(&divergence_matrix, &proteins)
        .with_context(|| format!("failed to read sequence matrix `{}`", divergence_matrix))?;
    let species2id = species_tree
        .leaves()
        .map(|l| {
            (
                species_tree[l]
                    .data
                    .name
                    .as_ref()
                    .expect("Found nameless leaf in species tree")
                    .to_owned(),
                l,
            )
        })
        .collect::<HashMap<_, _>>();

    info!("Storing gene data");
    let landscape_sizes = proteins
        .iter()
        .map(|p| {
            book.get(p)
                .map(|x| x.landscape.len())
                .with_context(|| format!("while looking up {}", p))
        })
        .collect::<Result<Vec<_>>>()
        .map_err(|e| anyhow!(format!("Can not find `{}` in database", e)))?;

    let mut core = proteins
        .iter()
        .enumerate()
        .filter_map(|(i, _)| {
            if landscape_sizes[i] > CORE_THRESHOLD {
                Some(i)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    if core.is_empty() {
        warn!("Core is empty; using all non-solo instead");
        core = (0..proteins.len())
            .filter(|&p| landscape_sizes[p] > 1)
            .collect::<Vec<_>>();
    }
    let species = proteins
        .iter()
        .map(|p| species2id[&underscore2point(&book.get(p).unwrap().species)])
        .collect::<Vec<SpeciesID>>();
    let all_species = HashSet::from_iter(species.iter().cloned());

    let extended = proteins
        .iter()
        .enumerate()
        .filter_map(|(i, _)| {
            if 1 < landscape_sizes[i] && landscape_sizes[i] <= CORE_THRESHOLD {
                Some(i)
            } else {
                None
            }
        })
        .collect::<HashSet<_>>()
        .difference(&HashSet::<GeneID>::from_iter(core.iter().copied()))
        .copied()
        .collect::<Vec<_>>();
    let solos = Vec::from_iter((0..proteins.len()).filter(|&p| landscape_sizes[p] == 1));
    let register = Register {
        size: proteins.len(),
        landscape_size: landscape_sizes,
        core,
        proteins: proteins.to_owned(),
        genes,
        species,
        all_species,
        synteny,
        divergence,
        species_tree,
        extended,
        solos,
    };
    Ok(register)
}

fn mean(x: &[f32]) -> f32 {
    x.iter().sum::<f32>() as f32 / x.len() as f32
}

fn find_threshold(register: &Register) -> f32 {
    const CLUSTER_THRESHOLD: usize = 1;

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
    tts.sort_by(|a, b| a.partial_cmp(&b).unwrap());
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
            clusters
                .iter_mut()
                .for_each(|c| c.iter_mut().for_each(|x| *x = register.core[*x]));

            clusters
        })
        .collect::<Vec<_>>();

    info!("Computing optimization targets");
    for (&t, clusters) in tts.iter().zip(clusterss.into_iter()) {
        let redundancies = clusters.iter().map(|c| c.len()).sum::<usize>() as f32
            / clusters
                .iter()
                .map(|c| {
                    view(&register.species, c)
                        .copied()
                        .collect::<HashSet<SpeciesID>>()
                        .len()
                })
                .sum::<usize>() as f32;

        let compacities = mean(
            &clusters
                .iter()
                .map(|c| {
                    let my_species =
                        HashSet::<SpeciesID>::from_iter(view(&register.species, c).copied());
                    my_species.len() as f32 / register.mrca_span(&my_species).len() as f32
                })
                .collect::<Vec<f32>>(),
        );

        reds.push(round(redundancies, 2));
        counts.push(clusters.len() as i64);
        avg_comp.push(round(compacities, 2));
    }

    info!("Finding optimal threshold");
    let mut cmpct_maxima = (1..tts.len() - 1)
        .filter_map(|i| {
            if avg_comp[i - 1] < avg_comp[i] && avg_comp[i] >= avg_comp[i + 1] {
                Some(i)
            } else {
                None
            }
        })
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
    debug!(
        "THR: {:.3} CNT: {:3} CMPCT: {:.2} RED: {:.2}",
        tts[i], counts[i], avg_comp[i], reds[i]
    );

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
    for i in 0..n {
        clusters.entry(r[i]).or_insert(vec![]).push(i);
    }

    clusters.into_values().collect::<Vec<Vec<usize>>>()
}

fn remove_solos_clusters(t: &mut PolytomicGeneTree) -> Vec<NodeID> {
    let solo_clusters = t[1]
        .children
        .iter()
        .copied()
        .filter(|&n| t[n].children.is_empty() && t[n].content.len() == 1)
        .collect::<Vec<_>>();
    let r = solo_clusters
        .iter()
        .copied()
        .flat_map(|n| t[n].content.iter())
        .copied()
        .collect::<Vec<_>>();
    t.delete_nodes(&solo_clusters);
    r
}

fn inject_satellites(
    t: &mut PolytomicGeneTree,
    satellites: &[NodeID],
    register: &Register,
) -> Vec<NodeID> {
    let mut cached_elcs = t
        .nodes()
        .copied()
        .filter(|c| !t[*c].content.is_empty())
        .map(|c| (c, register.elc(view(&register.species, &t[c].content))))
        .collect::<HashMap<_, _>>();
    let mut true_solos = vec![];
    let mut satellites = satellites.to_vec();

    satellites.sort_by_cached_key(|s| {
        cached_elcs
            .iter()
            .map(|(k, v)| {
                (
                    register.elc(
                        view(&register.species, &t[*k].content)
                            .chain([register.species[*s]].iter()),
                    ) - v,
                    -register
                        .species_tree
                        .node_topological_depth(register.species[*s]),
                    &register.proteins[*s],
                )
            })
            .min()
            .unwrap()
    });

    for &id in satellites.iter() {
        let log = register.proteins[id] == "ENSLLEP00000023490";
        let mut candidate_clusters = t
            .nodes()
            .copied()
            .filter(|c| !t[*c].content.is_empty())
            .sorted()
            .par_bridge()
            .filter(|c| {
                register
                    .mrca_span(view(&register.species, &t[*c].content))
                    .contains(&register.species[id])
            })
            .map(|c| {
                let c_content = &t[c].content;
                let synteny = register.synteny.masked(c_content, &[id]).max();
                let sequence = register.divergence.masked(c_content, &[id]).min();
                let delta_elc = register
                    .elc(view(&register.species, c_content).chain([register.species[id]].iter()))
                    - cached_elcs[&c];

                (
                    c,
                    (delta_elc, OrderedFloat(sequence), OrderedFloat(synteny)),
                )
            })
            .collect::<Vec<_>>();

        let syntenies = candidate_clusters
            .iter()
            .map(|c| c.1 .2)
            .collect::<Vec<_>>();
        let divergences = candidate_clusters
            .iter()
            .map(|c| c.1 .1)
            .collect::<Vec<_>>();

        if syntenies.iter().any(|&s| s >= RELAXED_SYNTENY_THRESHOLD) {
            candidate_clusters.sort_by_key(|c| {
                (
                    -c.1 .2, // Synteny
                    c.1 .0,  // ELC
                )
            });
        } else if divergences.iter().any(|d| *d <= OrderedFloat(0.5)) {
            candidate_clusters.sort_by_key(|c| {
                (
                    c.1 .0.max(0), // ELC
                    c.1 .1,        // Divergence
                )
            });
        } else {
            candidate_clusters.sort_by_key(|c| {
                (
                    c.1 .0,  // ELC
                    c.1 .1,  // Divergence
                    -c.1 .2, // Synteny
                )
            });
        }

        if log {
            info!("SATELLITE {}", register.proteins[id]);
            for cc in &candidate_clusters {
                let c = cc.1;
                info!(
                    "    {:20} --> #{:3} -- SYN: {:.2} ΔELC: {:3} DV: {:2.2}",
                    view(&register.proteins, &t[cc.0].content)
                        .sorted()
                        .next()
                        .unwrap(),
                    cc.0,
                    c.2,
                    c.0,
                    c.1
                );
            }
        }
        if let Some(cluster) = candidate_clusters.first() {
            let parent = cluster.0;
            t[parent].content.push(id);
            t[parent].tag = register
                .species_tree
                .mrca(view(&register.species, &t[parent].content))
                .unwrap();
            cached_elcs.insert(
                parent,
                register.elc(view(&register.species, &t[parent].content)),
            );
        } else {
            info!("Creating a new cluster for {}", register.proteins[id]);
            let new_cluster = t.add_node(&[id], register.species[id], Some(1));
            cached_elcs.insert(
                new_cluster,
                register.elc(view(&register.species, &t[new_cluster].content)),
            );
        }
    }

    true_solos
}

fn inject_solos(t: &mut PolytomicGeneTree, register: &Register) {
    for &id in register.solos.iter() {
        let log = ["ENSTMTP00000017759"].contains(&register.proteins[id].as_str());

        let mut candidate_clusters = t
            .nodes()
            .copied()
            .par_bridge()
            .filter(|c| !t[*c].content.is_empty())
            .map(|c| {
                let c_content = &t[c].content;
                let delta_elc = register
                    .elc(view(&register.species, c_content).chain([register.species[id]].iter()))
                    - register.elc(view(&register.species, c_content));
                let not_already_in =
                    view(&register.species, c_content).any(|s| *s == register.species[id]) as i8;
                let in_span = register
                    .mrca_span(view(&register.species, c_content))
                    .contains(&register.species[id]) as i8;
                let sequence = register.divergence.masked(&[id], c_content).min();
                (
                    c,
                    (delta_elc, in_span, OrderedFloat(sequence), not_already_in),
                )
            })
            .collect::<Vec<_>>();
        candidate_clusters.sort_by_key(|c| c.1);

        if log {
            info!("SOLO {}", register.proteins[id]);
            for cc in &candidate_clusters {
                let c = cc.1;
                info!(
                    "    {:20} -- ΔELC: {:3} IN:{} NARY: {} DV: {:2.2}",
                    view(&register.proteins, &t[cc.0].content)
                        .sorted()
                        .next()
                        .unwrap(),
                    c.0,
                    c.1,
                    c.3,
                    c.2
                );
            }
        }
        if let Some(cluster) = candidate_clusters.first() {
            let parent = cluster.0;
            t[parent].content.push(id);
            t[parent].tag = register
                .species_tree
                .mrca(view(&register.species, &t[parent].content))
                .unwrap();
        } else {
            debug!("{} is missing a home", &register.proteins[id]);
            t.add_node(&[id], register.species[id], Some(1));
        }
    }
}

fn inject_extended(t: &mut PolytomicGeneTree, register: &Register) {
    let mut extended = register.extended.to_owned();

    let mut local_synteny = register.synteny.clone();
    for i in extended.iter() {
        for j in 0..register.size {
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
                HashSet::<usize>::from_iter(t[c].content.iter().copied())
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
        .filter(|c| !t[*c].content.is_empty())
        .map(|c| (c, register.elc(view(&register.species, &t[c].content))))
        .collect::<HashMap<_, _>>();

    let mut touched = true;
    while touched {
        let mut new_solos = vec![];
        touched = false;

        for id in extended.iter() {
            let log = ["ENSNSUP00000004150"].contains(&register.proteins[*id].as_str());

            let mut candidate_clusters = t
                .nodes()
                .copied()
                .par_bridge()
                .filter(|c| !t[*c].content.is_empty())
                .map(|c| {
                    let c_content = &t[c].content;
                    let delta_elc = register.elc(
                        view(&register.species, c_content).chain([register.species[*id]].iter()),
                    ) - cached_elcs[&c];
                    let divergence = round(register.divergence.masked(&[*id], c_content).min(), 2);
                    let synteny = -round(local_synteny.masked(&[*id], &core_content[&c]).max(), 2);
                    (
                        c,
                        (delta_elc, OrderedFloat(divergence), OrderedFloat(synteny)),
                    )
                })
                .collect::<Vec<_>>();

            let syntenies = candidate_clusters
                .iter()
                .map(|c| c.1 .2)
                .collect::<Vec<_>>();
            let elcs = candidate_clusters
                .iter()
                .map(|c| c.1 .0)
                .collect::<Vec<_>>();

            if syntenies.iter().any(|&s| s <= -RELAXED_SYNTENY_THRESHOLD) {
                // Synteny is negated
                if elcs.iter().any(|&e| e < 3) {
                    candidate_clusters.retain(|c| c.1 .0 <= 3);
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
                info!("EXTENDED {}", register.proteins[*id],);
                for cc in &candidate_clusters {
                    let c = cc.1;
                    info!(
                        "    {:20} --> SYN: {:.2} DV: {:2.2} ΔELC: {} CRT: {}",
                        register.proteins[t[cc.0].content[0]],
                        c.2,
                        c.1,
                        c.0,
                        register.species_tree.node_topological_depth(
                            register
                                .species_tree
                                .mrca(view(&register.species, &t[cc.0].content))
                                .unwrap()
                        ) - register.species_tree.node_topological_depth(
                            register
                                .species_tree
                                .mrca(view(
                                    &register.species,
                                    t[cc.0].content.iter().chain([register.species[*id]].iter())
                                ))
                                .unwrap()
                        )
                    );
                }
            }

            if let Some(cluster) = candidate_clusters.first() {
                let parent = cluster.0;
                t[parent].content.push(*id);
                t[parent].tag = register
                    .species_tree
                    .mrca(view(&register.species, &t[parent].content))
                    .unwrap();
                cached_elcs.insert(
                    parent,
                    register.elc(view(&register.species, &t[parent].content)),
                );
            } else {
                new_solos.push(*id);
            }
        }
        extended = new_solos;
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
            species: HashSet::from_iter([seed_species].into_iter()),
        })
        .collect::<Vec<_>>();
    sources.get_mut(&seed_species).unwrap().clear();
    let log = view(&register.proteins, ds.iter().flat_map(|d| d.content.iter())).any(|s| {
        [
            "ENSDCDP00000006390",
            "ENSHHUP00000059113",
            "ENSOKIP00005017062",
            "ENSOMYP00000060796",
            "ENSOTSP00005012508",
            "ENSPKIP00000011431",
            "ENSSARP00000000739",
            "ENSSGRP00000084880",
            "ENSSSAP00000087432",
            "ENSSTUP00000100855",
        ]
        .contains(&s.as_str())
    });
    if log {
        dbg!(seed_species);
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
            synteny: OrderedFloat<f32>,
            species: SpeciesID,
        }
        let mut filled = HashMap::<(usize, SpeciesID), bool>::new();
        let mut links = Vec::new();
        let mut putatives = ds.iter().map(|_| Vec::new()).collect::<Vec<_>>();

        for (i, d) in ds.iter().enumerate() {
            let new_span = total_span
                .iter()
                .copied()
                .filter(|x| !d.species.contains(&x) && sources.contains_key(&x))
                .collect::<Vec<_>>();

            for species in new_span {
                let candidates = sources.entry(species).or_default();
                candidates.sort_by_cached_key(|&c| {
                    OrderedFloat(-register.synteny.masked(&[c], &d.content).max())
                });

                for c in candidates {
                    let synteny = OrderedFloat(register.synteny.masked(&[*c], &d.content).max());
                    links.push(Link {
                        arm: i,
                        gene: *c,
                        synteny,
                        species,
                    });
                }
            }
        }
        links.sort_by_key(|l| (-i64::try_from(ds[l.arm].content.len()).unwrap(), -l.synteny));

        for l in links {
            if !filled.get(&(l.arm, l.species)).unwrap_or(&false)
                && sources[&l.species].contains(&l.gene)
            {
                putatives[l.arm].push((l.gene, l.species, *l.synteny));
                sources
                    .get_mut(&l.species)
                    .unwrap()
                    .retain(|x| *x != l.gene);
                filled.insert((l.arm, l.species), true);
            }
        }

        for (i, d) in ds.iter_mut().enumerate() {
            for (new, species, _) in putatives[i].iter() {
                d.content.push(*new);
                d.species.insert(*species);
            }
        }

        // if log {
        //     for (i, pp) in putatives.iter().enumerate() {
        //         for p in pp.iter() {
        //             println!(
        //                 "{} gets - {}/{} {}",
        //                 i,
        //                 register.species_name(p.1),
        //                 register.proteins[p.0],
        //                 p.2
        //             );
        //         }
        //     }
        // }

        let mut dcss = VecMatrix::<f32>::with_elem(ds.len(), ds.len(), 0.);
        for i in 0..ds.len() {
            for j in i..ds.len() {
                dcss[(i, j)] = ds[i].species.len() as f32 / ds[j].species.len() as f32;
                dcss[(j, i)] = ds[j].species.len() as f32 / ds[i].species.len() as f32;
                // dcss[(i, j)] = ds[i].species.len() as f32 / ds[j].species.intersection(&ds[i].species).count() as f32;
                // dcss[(j, i)] = ds[j].species.len() as f32 / ds[i].species.intersection(&ds[j].species).count() as f32;
            }
        }

        for (i, d) in ds.iter_mut().enumerate() {
            let putative_dcss = (0..dcss.ncols())
                .map(|j| dcss[(i, j)])
                .filter(|&x| 0.6 <= x && x <= 1.5)
                .map(|x| OrderedFloat(x))
                .collect::<Vec<_>>();
            let best_dcs = if putative_dcss.is_empty() {
                OrderedFloat(-1.0)
            } else {
                *putative_dcss.iter().max().unwrap()
            };
            let filling = d.species.len() as f32
                / register
                    .mrca_span(&d.species)
                    .iter()
                    .copied()
                    .collect::<HashSet<_>>()
                    .intersection(&register.all_species)
                    .count() as f32;
            let filling_threshold = if total_span.len() <= 4 { 0.5 } else { 0.5 };

            if log {
                // dbg!(i, filling_threshold, filling);
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
        // if log { ds.iter().for_each(|d| d.pretty(register))}
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
        d.root = register.species_tree.mrca(&d.species).unwrap();
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
    for i in &tree[reference].content {
        sources.entry(register.species[*i]).or_default().push(*i);
    }

    let mut seed_speciess = sources
        .keys()
        .copied()
        .filter(|s| sources[s].len() > 1)
        .sorted_by_cached_key(|s| {
            (
                -(sources[s].len() as i64),
                -register.species_tree.node_topological_depth(*s),
            )
        })
        .collect::<Vec<_>>();
    let duplicated_species = sources
        .keys()
        .copied()
        .filter(|s| sources[s].len() > 1)
        .collect::<HashSet<_>>();

    while let Some(seed_species) = seed_speciess.get(0) {
        let new_family = grow_duplication(&mut sources, *seed_species, register);
        dups.push(new_family);
        sources.remove(seed_species);
        seed_speciess.retain(|s| sources.get(s).map(|v| v.len() > 1).unwrap_or(false));
    }

    assert!(sources.values().all(|l| l.len() <= 1));

    // Remaining singletons from previously duplicated species should be treated independently
    let remaining_species = HashSet::<SpeciesID>::from_iter(
        sources
            .iter()
            .filter_map(|(k, v)| if !v.is_empty() { Some(k) } else { None })
            .copied(),
    );
    if !remaining_species.is_empty() {
        for s in remaining_species.intersection(&duplicated_species) {
            dups.push(vec![Duplication {
                root: register.species_tree.mrca(&[*s]).unwrap(),
                content: sources[s].clone(),
                species: HashSet::from_iter([*s].into_iter()),
            }]);
            sources.remove(s);
        }
    }

    // If the remaining species make for a subset of the existing duplications, they should go under
    let remaining_species = HashSet::from_iter(
        sources
            .iter()
            .filter_map(|(k, v)| if !v.is_empty() { Some(k) } else { None })
            .copied(),
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
        let remaining_genes = sources
            .values()
            .filter(|v| !v.is_empty())
            .map(|v| v[0])
            .collect::<Vec<_>>();
        dups.push(vec![Duplication {
            root: register.species_tree.mrca(&remaining_species).unwrap(),
            content: remaining_genes,
            species: remaining_species.into_iter().collect(),
        }]);
    }

    dups
}

fn resolve_duplications(t: &mut PolytomicGeneTree, register: &Register) {
    let roots = t[1].children.clone();

    for &i in &roots {
        let log = false; // view(&register.proteins, &t[i].content).any(|x| x == "ENSMSIP00000026938");
        let mut dups = create_duplications(register, t, i);
        dups.sort_by_cached_key(|f| {
            let all_species = f
                .iter()
                .flat_map(|d| d.content.iter().map(|&x| register.species[x]));
            let unique_species = HashSet::<usize>::from_iter(all_species).len() as i64;
            -unique_species
        });

        let cluster_mrca = register
            .species_tree
            .mrca(view(&register.species, &t[i].content))
            .unwrap();
        let my_owns = t[i]
            .content
            .iter()
            .copied()
            .filter(|g| {
                dups.iter()
                    .flat_map(|ds| ds.iter())
                    .all(|d| !d.content.contains(g))
            })
            .collect::<Vec<_>>();
        t[i].content.clear();
        reconcile(t, &my_owns, register, Some(i), Some(cluster_mrca), true);

        for f in dups.iter_mut() {
            // We use pop, i.e. implicit reverse
            f.sort_by_cached_key(|d| {
                d.content
                    .iter()
                    .unique_by(|g| register.species[**g])
                    .count()
            });

            while let Some(d) = f.pop() {
                let log = view(&register.proteins, &d.content).any(|x| x == "ENSMSIP00000026938");
                if log {
                    println!(
                        "{} {}",
                        register.proteins[d.content[0]],
                        register.species_name(d.root)
                    );
                }
                let mut root_candidates = t
                    .descendants(i)
                    .iter()
                    .copied()
                    .filter(|n| t[*n].tag == d.root)
                    .sorted_by_cached_key(|n| {
                        let leaves = t.descendant_leaves(*n);
                        let synteny = if leaves.is_empty() {
                            OrderedFloat(0.0)
                        } else {
                            OrderedFloat(register.synteny.masked(&d.content, &leaves).max())
                        };
                        let t_depth = t.topo_depth(*n) as i64;
                        (-t_depth, -synteny)
                    })
                    .collect::<Vec<_>>();

                if log {
                    for &c in root_candidates.iter() {
                        println!(
                            "TOPODPTH: {:3} SYN: {:2.2}",
                            t.topo_depth(c),
                            register
                                .synteny
                                .masked(&d.content, &t.descendant_leaves(c))
                                .max()
                        );
                    }
                }

                let current_root = if let Some(root) = root_candidates.first() {
                    *root
                } else {
                    i
                };

                // If the putative root is a leaf node, make it a non-leaf node
                if !t[current_root].content.is_empty() {
                    let content = t.add_node(
                        &t[current_root].content.clone(),
                        register
                            .species_tree
                            .mrca(
                                t[current_root]
                                    .content
                                    .iter()
                                    .map(|g| &register.species[*g]),
                            )
                            .unwrap(),
                        Some(current_root),
                    );
                    t[current_root].content.clear();
                }

                // If the putative root has a single child, we can reconcile directly in it
                // Otherwise, we need to create an intermediate node
                if t[current_root].children.len() > 1 {
                    let node_alpha = t.add_node(&vec![], t[current_root].tag, None);
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
                        false,
                    );
                } else {
                    reconcile(t, &d.content, register, Some(current_root), None, false);
                }
            }
        }
    }
}

fn make_final_tree(t: &mut PolytomicGeneTree, register: &Register) -> usize {
    t.cache_descendants(1);
    let mut todo = t[1]
        .children
        .iter()
        .copied()
        .sorted_by_cached_key(|x| {
            (
                (t.descendant_leaves(*x).len() as i64),
                // -register.species_tree.node_topological_depth(t[*x].tag),
            )
        }) // .pop(), so implicit reverse ordering
        .collect::<Vec<_>>();

    let aa = todo.pop().unwrap();
    let new_root = reconcile_upstream(t, aa, register, None);
    t.cache_descendants(new_root);

    while let Some(a) = todo.pop() {
        let candidate_parents = t
            .descendants(new_root)
            .into_iter()
            .filter(|&b| b != new_root && b != a && t[b].tag == t[a].tag)
            .filter(|&b| !t.cached_descendants(a).unwrap().contains(&b))
            .collect::<Vec<_>>();

        let mut leaves = HashMap::new();
        leaves.insert(
            a,
            vec_inter(&t.descendant_leaves(a), &register.core)
                .into_iter()
                .collect::<Vec<_>>(),
        );
        let mut speciess = HashMap::new();
        speciess.insert(
            a,
            view(&register.species, &leaves[&a])
                .copied()
                .collect::<HashSet<SpeciesID>>(),
        );
        let mut context_nodes = HashMap::new();
        for i in &candidate_parents {
            let mut context_node = *i;
            'search_context: while t.descendant_leaves(context_node).len() as f32
                <= leaves[&a].len() as f32 / 2.0 + 1.
                && t[context_node].parent.is_some()
            {
                context_node = t[context_node].parent.unwrap();
                if candidate_parents
                    .iter()
                    .filter(|&c| c != i)
                    .any(|c| t.cached_descendants(context_node).unwrap().contains(c))
                {
                    break 'search_context;
                }
            }
            context_nodes.insert(*i, context_node);
            leaves.insert(
                *i,
                vec_inter(&t.descendant_leaves(context_node), &register.core)
                    .into_iter()
                    .collect::<Vec<_>>(),
            );
            speciess.insert(
                *i,
                t.descendant_leaves(context_node)
                    .iter()
                    .map(|g| register.species[*g])
                    .collect::<HashSet<_>>(),
            );
        }

        let log = view(&register.proteins, &t.descendant_leaves(a))
            .any(|s| ["ENSACCP00020023244"].contains(&s.as_str()));

        let mut candidate_parents = candidate_parents
            .par_iter()
            .filter(|b| !leaves[b].is_empty())
            .map(|b| {
                let missing_left = speciess[&a].difference(&speciess[&b]);
                let missing_right = speciess[&b].difference(&speciess[&a]);
                let elc = register.elc(missing_left) + register.elc(missing_right);
                let synteny = register.synteny.masked(&leaves[&a], &leaves[b]).max();
                let divergence = register.divergence.masked(&leaves[&a], &leaves[b]).min();

                (
                    b,
                    (
                        elc,
                        OrderedFloat(round(synteny, 2)),
                        OrderedFloat(round(divergence, 2)),
                    ),
                )
            })
            .collect::<Vec<_>>();
        candidate_parents.sort_by_key(|c| (c.1 .0, -c.1 .1, c.1 .2));

        if log {
            println!("\n\n=== BEFORE ===");
            for b in &candidate_parents {
                let bb = b.1;
                let b = b.0;
                info!(
                    "{:20} {:4}/{:4} ({:3}/{:3} leaves) DELC: {:2} SYN: {:2.2} DV: {:2.2} T: {} {}",
                    register.species_name(t[*b].tag),
                    b,
                    context_nodes[b],
                    t.descendant_leaves(*b).len(),
                    leaves[b].len(),
                    bb.0,
                    bb.1,
                    bb.2,
                    register.species_name(t[*b].tag),
                    leaves[b]
                        .get(0)
                        .map(|g| register.proteins[*g].as_str())
                        .unwrap_or("EMPTY")
                );
            }
        }

        if candidate_parents
            .iter()
            .map(|c| c.1 .1)
            .any(|s| s >= RELAXED_SYNTENY_THRESHOLD)
        {
            // debug!("Filtering by synteny");
            candidate_parents.retain(|c| c.1 .1 >= RELAXED_SYNTENY_THRESHOLD);
            candidate_parents.sort_by_key(|c| {
                (
                    -c.1 .1, // Synteny
                    c.1 .0,  // ELC
                    c.1 .2,  // Sequence
                )
            })
        }
        if candidate_parents
            .iter()
            .map(|c| c.1 .2)
            .any(|d| d <= OrderedFloat(0.5))
        {
            // debug!("Filtering by divergence");
            candidate_parents.retain(|c| c.1 .2 <= OrderedFloat(0.5));
            candidate_parents.sort_by_key(|c| {
                (
                    c.1 .2,  // Sequence
                    c.1 .0,  // ELC
                    -c.1 .1, // Synteny
                )
            })
        }

        if log {
            println!("=== AFTER ===");
            for b in &candidate_parents {
                let bb = b.1;
                let b = b.0;
                info!(
                    "{:4}/{:4} ({:4} leaves) DELC: {:2} SYN: {:2.2} DV: {:2.2} T: {} {}",
                    b,
                    context_nodes[b],
                    leaves[b].len(),
                    bb.0,
                    bb.1,
                    bb.2,
                    register.species_name(t[*b].tag),
                    leaves[b]
                        .get(0)
                        .map(|g| register.proteins[*g].as_str())
                        .unwrap_or("EMPTY")
                );
            }
        }
        if let Some(cluster) = candidate_parents.first() {
            let child = a;
            let parent = *cluster.0;

            if log {
                println!(
                    "Parent: {}/{}",
                    t[parent].content.len(),
                    t[parent].children.len()
                );
            }

            // If the child is plugged on a completely empty subtree, it should replace it to lock its docking points
            if t.descendant_leaves(parent).is_empty() {
                t.move_node(child, t[parent].parent.unwrap());
                t.delete_node(parent);
                t.cache_descendants(t[child].parent.unwrap());
            } else {
                if !t[parent].content.is_empty() {
                    let _content = t.add_node(
                        &t[parent].content.clone(),
                        register
                            .species_tree
                            .mrca(view(&register.species, &t[parent].content))
                            .unwrap(),
                        Some(parent),
                    );
                    t[parent].content.clear();
                }

                if t[parent].children.len() > 1 {
                    let node_alpha = t.add_node(&[], t[parent].tag, None);
                    for c in t[parent].children.clone().into_iter() {
                        t.move_node(c, node_alpha);
                    }
                    t.move_node(node_alpha, parent);
                }

                t.move_node(child, parent);
                t.cache_descendants(parent);
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
                        &view(&register.species, t.descendant_leaves(*c).iter())
                            .copied()
                            .collect::<HashSet<_>>(),
                        &view(&register.species, t.descendant_leaves(*o).iter())
                            .copied()
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
                .mrca(view(&register.species, &t.descendant_leaves(node_alpha)))
                .unwrap();
        } else {
            t.move_node(to_plug, root);
        }

        if !t.descendant_leaves(root).is_empty() {
            t[root].tag = register
                .species_tree
                .mrca(view(&register.species, &t.descendant_leaves(root)))
                .unwrap();
        }
    }

    assert!(t.descendant_leaves(1).len() == 0);
    assert!(t[1].children.len() == 0);
    assert!(t[1].content.len() == 0);
    t.delete_node(1);

    for k in t.nodes().copied().collect::<Vec<_>>().iter() {
        if !t.descendant_leaves(*k).is_empty() {
            t[*k].tag = register
                .species_tree
                .mrca(view(&register.species, &t.descendant_leaves(*k)))
                .unwrap();
        }
    }
    new_root
}

fn prune_tree(tree: &mut PolytomicGeneTree, root: usize) {
    info!("Pruning empty leaf nodes -- from {}", tree.nodes().count());
    loop {
        let todos = tree
            .nodes()
            .copied()
            .filter(|&k| tree[k].children.is_empty() && tree[k].content.is_empty())
            .collect::<Vec<_>>();
        if todos.is_empty() {
            break;
        } else {
            tree.delete_nodes(&todos);
        }
    }
    info!("...to {}", tree.nodes().count());

    info!("Pruning empty leaf nodes -- from {}", tree.nodes().count());
    loop {
        let todos = tree
            .nodes()
            .copied()
            .filter(|&k| tree[k].children.is_empty() && tree[k].content.len() == 1 && k != root)
            .collect::<Vec<_>>();

        if let Some(k) = todos.get(0) {
            let content = tree[*k].content[0];
            let parent = tree[*k].parent.unwrap();
            tree[parent].content.push(content);
            tree.delete_node(*k);
        } else {
            break;
        }
    }
    info!("...to {}", tree.nodes().count());

    info!("Pruning scaffolding -- from {}", tree.nodes().count());
    // Non used speciation nodes
    loop {
        let todos = tree
            .nodes()
            .copied()
            .filter(|&k| k != root && tree[k].children.len() == 1 && tree[k].content.is_empty())
            .collect::<Vec<_>>();
        if let Some(k) = todos.get(0) {
            tree.move_node(tree[*k].children[0], tree[*k].parent.unwrap());
            tree.delete_node(*k);
        } else {
            break;
        }
    }
    info!("...to {}", tree.nodes().count());
}

fn reconcile(
    t: &mut PolytomicGeneTree,
    ps: &[GeneID],
    register: &Register,
    graft: Option<usize>,
    from: Option<SpeciesID>,
    full: bool,
) -> usize {
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
            let to_plug = ps
                .iter()
                .filter(|p| register.species[**p] == phylo_node)
                .collect::<Vec<_>>();
            if !to_plug.is_empty() {
                for &new in to_plug.into_iter() {
                    current_root = t.add_node(&[new], phylo_node, Some(current_root));
                }
            } else if full {
                t.add_node(&[], phylo_node, Some(current_root));
            }
        } else {
            let actual_children = register
                .species_tree
                .children(phylo_node)
                .into_iter()
                .filter(|child| {
                    full || !speciess.is_disjoint(&HashSet::from_iter(
                        register.species_tree.leaves_of(**child).iter().copied(),
                    ))
                })
                .collect::<Vec<_>>();
            if actual_children.is_empty() {
                return;
            } else {
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
    let mrca = from.unwrap_or_else(|| register.species_tree.mrca(&speciess).unwrap());
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
            for child in register.species_tree.children(phylo_node).into_iter() {
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

fn do_family(
    family: &[String],
    id: &str,
    register: &Register,
    logs_root: &str,
) -> Result<PolytomicGeneTree> {
    info!("===== Family {} -- {} proteins =====", id, family.len());

    info!("Optimizing threshold");
    let tt = find_threshold(&register);

    let mut tree = PolytomicGeneTree::new();
    let root = tree.add_node(&[], 0, None);
    let clusters = clusterize(&register.synteny.masked(&register.core, &register.core), tt);
    for c in clusters.into_iter() {
        tree.add_node(
            &c.iter().map(|x| register.core[*x]).collect::<Vec<_>>(),
            0,
            Some(root),
        );
    }
    File::create(&format!("{}/{}_vanilla.nwk", &logs_root, id))?.write_all(
        &tree
            .to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t))
            .as_bytes(),
    )?;

    info!(
        "Injecting {} extended in {} clusters",
        register.extended.len(),
        tree[1].children.len()
    );
    inject_extended(&mut tree, &register);
    let satellites = remove_solos_clusters(&mut tree);
    File::create(&format!("{}/{}_withextended.nwk", &logs_root, id))?.write_all(
        &tree
            .to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t))
            .as_bytes(),
    )?;

    let mut clusters = tree[root]
        .children
        .iter()
        .copied()
        .sorted_by_cached_key(|c| {
            let topo_depth = register.species_tree.node_topological_depth(
                register
                    .species_tree
                    .mrca(view(&register.species, &tree[*c].content))
                    .unwrap(),
            );
            let size = tree[*c].content.len();
            (-topo_depth, size)
        })
        .collect::<Vec<_>>();

    assert!(tree[1]
        .children
        .iter()
        .all(|k| !tree[*k].content.is_empty()));
    info!("Packing clusters");
    info!("From {}...", tree[1].children.len());
    // sim_matrix(&tree, &register, &clusters, &|i, j| {
    //     jaccard(
    //         &view(&register.species, &tree[i].content)
    //             .copied()
    //             .collect::<HashSet<_>>(),
    //         &view(&register.species, &tree[j].content)
    //             .copied()
    //             .collect::<HashSet<_>>(),
    //     )
    // });
    struct Bin {
        members: Vec<usize>,
        species: HashSet<SpeciesID>,
        content: Vec<GeneID>,
    }
    if let Some(c) = clusters.pop() {
        let mut bins = vec![Bin {
            members: vec![c],
            species: HashSet::from_iter(view(&register.species, &tree[c].content).copied()),
            content: tree[c].content.clone(),
        }];

        while let Some(k) = clusters.pop() {
            let mut candidates = bins
                .iter_mut()
                .filter(|b| {
                    jaccard(
                        &b.species,
                        &HashSet::from_iter(view(&register.species, &tree[k].content).copied()),
                    ) <= 0.01
                })
                .sorted_by_cached_key(|b| {
                    (
                        OrderedFloat(-round(
                            register.synteny.masked(&b.content, &tree[k].content).max(),
                            2,
                        )),
                        OrderedFloat(round(
                            register
                                .divergence
                                .masked(&b.content, &tree[k].content)
                                .min(),
                            2,
                        )),
                        -(b.species.len() as i64),
                    )
                });

            if let Some(b) = candidates.next() {
                b.members.push(k);
                b.species.extend(view(&register.species, &tree[k].content));
                b.content.extend_from_slice(&tree[k].content);
            } else {
                bins.push(Bin {
                    members: vec![k],
                    species: HashSet::from_iter(view(&register.species, &tree[k].content).copied()),
                    content: tree[k].content.clone(),
                })
            }
        }

        for mut bin in bins.into_iter() {
            let e = register.elc(&bin.species) as f32;
            let ee = register.ellc(&bin.species);
            let merged_compactness = bin.species.len() as f32
                / (HashSet::<SpeciesID>::from_iter(
                    register.species_tree.mrca(&bin.species).iter().copied(),
                )
                .intersection(&HashSet::from_iter(register.all_species.iter().copied()))
                .count() as f32);
            // let all_compactnesses = todo!();
            let all_elcs = bin
                .members
                .iter()
                .map(|k| {
                    register.elc(&HashSet::<SpeciesID>::from_iter(
                        view(&register.species, &tree[*k].content).copied(),
                    ))
                })
                .collect::<Vec<_>>();

            if e <= (1.1 * all_elcs.iter().sum::<i64>() as f32).ceil() {
                if let Some(merger) = bin.members.pop() {
                    while let Some(merged) = bin.members.pop() {
                        tree.merge_nodes(merger, merged, &|a: &mut Vec<GeneID>, b: &[GeneID]| {
                            a.extend_from_slice(b)
                        });
                    }
                }
            }
        }
    }

    // sim_matrix(&tree, &register, &tree[1].children, &|i, j| {
    //     jaccard(
    //         &view(&register.species, &tree[i].content)
    //             .copied()
    //             .collect::<HashSet<_>>(),
    //         &view(&register.species, &tree[j].content)
    //             .copied()
    //             .collect::<HashSet<_>>(),
    //     )
    // });
    info!("...to {}", tree[1].children.len());
    File::create(&format!("{}/{}_merged.nwk", &logs_root, id))?.write_all(
        &tree
            .to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t))
            .as_bytes(),
    )?;

    info!(
        "Injecting {} satellites in {} clusters...",
        satellites.len(),
        tree[1].children.len()
    );
    let satellites = inject_satellites(&mut tree, &satellites, &register);
    info!("Done.");
    File::create(&format!("{}/{}_withsatellites.nwk", &logs_root, id))?.write_all(
        &tree
            .to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t))
            .as_bytes(),
    )?;

    info!("Tagging tree...");
    let to_remove = tree
        .nodes()
        .copied()
        .filter(|&k| k != root && tree[k].content.is_empty() && tree[k].content.is_empty())
        .collect::<Vec<_>>();
    tree.delete_nodes(&to_remove);
    for k in tree[root].children.clone().into_iter() {
        tree[k].tag = register
            .species_tree
            .mrca(view(&register.species, &tree[k].content))
            .unwrap();
    }
    info!("Done.");

    info!(
        "Injecting {} solos in {} clusters",
        register.solos.len(),
        tree[1].children.len()
    );
    inject_solos(&mut tree, &register);
    info!("Done.");
    File::create(&format!("{}/{}_clusters.nwk", &logs_root, id))?.write_all(
        &tree
            .to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t))
            .as_bytes(),
    )?;

    info!("Resolving duplications...");
    resolve_duplications(&mut tree, &register);
    info!("Done.");
    File::create(&format!("{}/{}_prototree.nwk", &logs_root, id))?.write_all(
        &tree
            .to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t))
            .as_bytes(),
    )?;

    info!(
        "Assembling final tree -- {} subtrees",
        tree[root].children.len()
    );

    let root = make_final_tree(&mut tree, &register);
    info!("Done.");
    File::create(&format!("{}/{}_reconciled.nwk", &logs_root, id))?.write_all(
        &tree
            .to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t))
            .as_bytes(),
    )?;

    info!("Pruning tree...");
    prune_tree(&mut tree, root);
    info!("Done.");

    for k in tree.nodes() {
        if tree[*k].children.len() > 2 {
            warn!(
                "Final tree is polytomic at {}",
                register.species_name(tree[*k].tag)
            );
        }
    }

    File::create(&format!("{}/{}_synteny.nwk", &logs_root, id))?.write_all(
        &tree
            .to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t))
            .as_bytes(),
    )?;

    if register.size != tree.descendant_leaves(root).len() {
        error!(
            "Gene mismatch: expected {}, found {}",
            register.size,
            tree.descendant_leaves(root).len()
        );
        let mut mine = tree.descendant_leaves(root);
        mine.sort();
        for gs in mine.windows(2) {
            if gs[0] == gs[1] {
                warn!("Double: {}", register.proteins[gs[0]]);
            }
        }
        let mine = HashSet::<String>::from_iter(
            tree.descendant_leaves(root)
                .iter()
                .map(|l| register.proteins[*l].to_owned()),
        );
        let others = HashSet::from_iter(register.proteins.iter().cloned());
        let missings = others.difference(&mine).collect::<Vec<_>>();
        info!("{:?}", missings);
        return Err(anyhow!("genes mismatch"));
    }

    Ok(tree)
}

pub fn do_file(
    filename: &str,
    batch: &str,
    book: &GeneBook,
    speciestree_file: &str,
    syntenies: &str,
    divergences: &str,
) -> Result<String> {
    let mut species_tree = newick::one_from_filename(&speciestree_file)
        .with_context(|| format!("can not open `{}`", speciestree_file))?;
    species_tree.cache_leaves();

    let gene_families =
        read_genefile(filename).with_context(|| format!("while parsing {}", filename))?;
    if gene_families.is_empty() || gene_families.len() > 1 {
        bail!("{} should contain a single family", filename)
    }
    let family = &gene_families[0];

    let now = Instant::now();

    let id = &std::path::Path::new(filename)
        .file_stem()
        .with_context(|| "asdfasdf")?
        .to_str()
        .with_context(|| format!("`{}` is not a valid file name", filename))?;
    let logs_root = format!("logs/{}/", batch);
    std::fs::create_dir_all(&logs_root)?;

    let register = make_register(id, &family, &book, &species_tree, syntenies, divergences)?;
    let out_tree = do_family(family, id, &register, &logs_root)?;

    info!("Done in {:.2}s.", now.elapsed().as_secs_f32());
    Ok(out_tree.to_newick(&|l| register.make_label(*l), &|t| register.species_name(*t)))
}
