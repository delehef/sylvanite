use std::{
    collections::{HashMap, HashSet},
    io::BufWriter,
};

use anyhow::*;
use identity_hash::IntMap;
use itertools::Itertools;
use log::*;
use logging_timer::time;
use newick::{NewickTree, NodeID};
use ordered_float::NotNan;
use rayon::prelude::*;
use syntesuite::genebook::GeneBook;

use crate::{
    dede::{Matrix, VecMatrix},
    errors::{self, RuntimeError},
    utils::{parse_dist_matrix_from_sim, read_genefile},
};

/// A phylogeny representing a .newick file.
#[derive(PartialEq, Debug)]
pub struct Phylogeny {
    /// The name of the current node.
    ///
    /// Can be empty for internal nodes.
    name: String,

    /// The children of the current node.
    ///
    /// Empty for leafs, and distances to the parent are optional.
    children: Vec<(Phylogeny, Option<f32>)>,
}

impl Phylogeny {
    /// Create a new phylogeny node from a name and list of children.
    pub fn new(name: &str, children: Vec<(Phylogeny, Option<f32>)>) -> Phylogeny {
        Phylogeny { name: name.to_string(), children }
    }

    /// Create a new leaf node.
    pub fn new_leaf(name: &str) -> Phylogeny {
        Phylogeny::new(name, vec![])
    }

    /// Join two trees into a new named parent node, with given edge lengths.
    pub fn join_with_name(name: &str, l: Phylogeny, dl: f32, r: Phylogeny, dr: f32) -> Phylogeny {
        Phylogeny::new(name, vec![(l, Some(dl)), (r, Some(dr))])
    }

    /// Join two trees into a new anonymous parent node, with given edge lengths.
    pub fn join(l: Phylogeny, dl: f32, r: Phylogeny, dr: f32) -> Phylogeny {
        Phylogeny::join_with_name("", l, dl, r, dr)
    }

    fn to_newick(&self) -> NewickTree {
        fn rec_insert(p: &Phylogeny, t: &mut NewickTree, me: NodeID) {
            for (c, l) in p.children.iter() {
                let c_id = t.add_node(
                    Some(me),
                    newick::Data {
                        name: if c.name.is_empty() { None } else { Some(c.name.clone()) },
                        attrs: Default::default(),
                    },
                );

                t.get_mut(c_id).unwrap().set_branch(l.unwrap_or(1.));
                rec_insert(c, t, c_id);
            }
        }

        let mut t = NewickTree::new();
        let root = t.add_node(
            None,
            newick::Data {
                name: if self.name.is_empty() { None } else { Some(self.name.clone()) },
                attrs: Default::default(),
            },
        );
        rec_insert(self, &mut t, root);
        t
    }
}

fn nj(m: &VecMatrix<f32>, ids: &[String]) -> Phylogeny {
    let mut parts: Vec<Option<Phylogeny>> =
        ids.iter().map(|name| Some(Phylogeny::new_leaf(name))).collect();

    let mut distances = m.clone();
    let mut active: Vec<usize> = (0..distances.nrows()).collect();

    while active.len() > 2 {
        let sum_d = |i: usize| -> f32 { active.iter().map(|&k| distances[(i, k)]).sum::<f32>() };

        // Find i,j for which Q(i,j) is minimal.
        let q = |&(&i, &j): &(&usize, &usize)| -> NotNan<f32> {
            let r =
                NotNan::new((active.len() - 2) as f32 * distances[(i, j)] - sum_d(i) - sum_d(j))
                    .unwrap();
            r
        };

        // Find minimal distance pair.
        let (&i, &j) = active
            .iter()
            .cartesian_product(active.iter())
            .filter(|&(&i, &j)| i != j)
            .min_by_key(q)
            .unwrap();

        let (i, j) = (i.min(j), i.max(j));

        // Compute distance from merged vertex to the nodes being merged.
        let di = distances[(i, j)] / 2. + (sum_d(i) - sum_d(j)) / (2. * (active.len() as f32 - 2.));
        let dj = distances[(i, j)] - di;

        // Remove j from positions considered in later iterations.
        active.remove(active.iter().position(|&x| x == j).unwrap());

        // Compute all other distances.
        active.iter().filter(|&&k| k != i).for_each(|&k| {
            let dk = (distances[(i, k)] + distances[(j, k)] - distances[(i, j)]) / 2.;
            distances[(i, k)] = dk;
            distances[(k, i)] = dk;
        });
        parts[i] =
            Some(Phylogeny::join(parts[i].take().unwrap(), di, parts[j].take().unwrap(), dj));
    }

    // Merge the two remaining vertices with the given distance.
    if let [i, j] = active[..] {
        let d = distances[(i, j)] / 2.;
        parts[i] = Some(Phylogeny::join(parts[i].take().unwrap(), d, parts[j].take().unwrap(), d));
    }

    parts[0].take().unwrap()
}

#[time]
fn make_approximate_gene_tree(f: &str, book: &GeneBook, syntenies: &str) -> Result<NewickTree> {
    let family = read_genefile(f).with_context(|| anyhow!("while parsing {}", f))?;
    let family_species =
        family.iter().map(|g| book.get(g).map(|r| r.species)).collect::<Result<Vec<_>>>()?;
    let id = &std::path::Path::new(f)
        .file_stem()
        .unwrap()
        .to_str()
        .ok_or_else(|| errors::FileError::InvalidFilename(format!("{:?}", f)))?;

    let synteny_matrix = &format!("{}/{}.dist", syntenies, id);
    let synteny = parse_dist_matrix_from_sim(synteny_matrix, &family).map_err(|e| {
        RuntimeError::FailedToReadMatrix { source: e, filename: synteny_matrix.to_owned() }
    })?;

    info!("===== Family {} -- {} proteins =====", id, family.len());
    Ok(nj(&synteny, &family_species).to_newick())
}

const WINDOW_SIZE: usize = 15;
pub(crate) fn build_species_tree(
    db_file: &str,
    bags_files: &[String],
    syntenies: &str,
) -> Result<NewickTree> {
    let book = &GeneBook::in_memory(db_file, WINDOW_SIZE, "id")?;
    let species = book.species();
    let species2id: HashMap<String, usize> =
        species.iter().enumerate().map(|(i, s)| (s.to_owned(), i)).collect();

    // 1. Create one approximate gene tree per bag
    info!("Approximating gene trees");
    let genes_trees = bags_files
        .par_iter()
        .map(|f| {
            make_approximate_gene_tree(f, &book, syntenies)
                .with_context(|| anyhow!("while processing {}", f))
        })
        .collect::<Result<Vec<_>>>()?;

    // 2. Generate an aggregated distance matrix from these trees
    info!("Computing meta-matrix");
    let mut meta_matrix = VecMatrix::<f32>::new_zero(species.len(), species.len());
    fn rec_insert(
        t: &NewickTree,
        n: NodeID,
        m: &mut VecMatrix<f32>,
        species2id: &HashMap<String, usize>,
    ) {
        let children = t.children(n).unwrap();
        // TODO: skip duplicated nodes?
        for c in children {
            let species = t
                .leaves_of(*c)
                .iter()
                .map(|l| t.get(*l).unwrap().data().name.clone().unwrap())
                .collect::<HashSet<_>>();
            for s1 in species.iter() {
                for s2 in species.iter() {
                    if s1 != s2 {
                        let i = species2id[s1];
                        let j = species2id[s2];
                        m[(i, j)] += 1.; // TODO: normalize by tree size?
                    }
                }
            }

            rec_insert(t, *c, m, species2id);
        }
    }

    for tree in genes_trees.iter() {
        rec_insert(tree, tree.root(), &mut meta_matrix, &species2id)
    }
    let max_distance = meta_matrix.iter().fold(0f32, |a, &b| a.max(b)); // TODO: replace with NaNFloat
    for i in 0..meta_matrix.nrows() {
        for j in 0..meta_matrix.ncols() {
            meta_matrix[(i, j)] = 1. - meta_matrix[(i, j)] / max_distance;
        }
    }
    crate::utils::write_matrix(
        &meta_matrix,
        &species,
        BufWriter::new(Box::new(std::fs::File::create("meta.dist").unwrap())),
    )?;

    // 3. Use NJ to generate the final tree
    info!("Computing species tree");
    Ok(nj(&meta_matrix, &species).to_newick())
}
