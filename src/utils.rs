use anyhow::*;
use log::*;
use rusqlite::Connection;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Read;
use std::io::Write;
use std::sync::Mutex;

use crate::errors;

#[allow(dead_code)]
pub enum GeneBook {
    InMemory(HashMap<String, Gene>),
    Cached(HashMap<String, Gene>),
    Inline(Mutex<Connection>, usize),
}

#[derive(Clone, Default)]
pub struct Gene {
    pub species: String,
    pub landscape: Vec<usize>,
}

impl GeneBook {
    fn parse_landscape(landscape: &str) -> Vec<usize> {
        if landscape.is_empty() {
            Vec::new()
        } else {
            landscape.split('.').map(|x| x.parse::<usize>().unwrap()).collect::<Vec<_>>()
        }
    }

    fn get_rows<P: rusqlite::Params>(
        mut query: rusqlite::Statement,
        params: P,
        window: usize,
    ) -> Result<HashMap<String, Gene>> {
        let genes = query
            .query_map(params, |r| {
                std::result::Result::Ok((
                    r.get::<_, String>(0)?, // id
                    r.get::<_, String>(1)?, // left tail
                    r.get::<_, String>(2)?, // right tail
                    r.get::<_, usize>(3)?,  // ancestral id
                    r.get::<_, String>(4)?, // species
                ))
            })?
            .collect::<Result<Vec<_>, _>>()?;

        Ok(genes
            .into_iter()
            .map(|g| {
                let mut left_landscape = Self::parse_landscape(&g.1);
                left_landscape.reverse();
                left_landscape.truncate(window);
                left_landscape.reverse();

                let mut right_landscape = Self::parse_landscape(&g.2);
                right_landscape.truncate(window);

                (
                    g.0.clone(),
                    Gene {
                        species: g.4,
                        landscape: left_landscape
                            .into_iter()
                            .chain([g.3].into_iter())
                            .chain(right_landscape.into_iter())
                            .collect(),
                    },
                )
            })
            .collect())
    }

    pub fn in_memory(filename: &str, window: usize) -> Result<Self> {
        info!("Caching the database...");

        let conn = Connection::open(filename).map_err(|e| errors::DataError::FailedToConnect {
            source: e,
            filename: filename.into(),
        })?;
        let query = conn.prepare(
            "SELECT id, left_tail_ids, right_tail_ids, ancestral_id, species FROM genomes",
        )?;
        let r = Self::get_rows(query, [], window)?;
        info!("Done.");
        Ok(GeneBook::InMemory(r))
    }

    pub fn cached<S: AsRef<str>>(filename: &str, window: usize, ids: &[S]) -> Result<Self> {
        info!("Caching the database...");

        let conn = Connection::open(filename).map_err(|e| errors::DataError::FailedToConnect {
            source: e,
            filename: filename.into(),
        })?;

        let query = conn.prepare(&format!(
            "SELECT id, left_tail_ids, right_tail_ids, ancestral_id, species FROM genomes WHERE id IN ({})",
            std::iter::repeat("?").take(ids.len()).collect::<Vec<_>>().join(", ")
        ))?;
        let r = Self::get_rows(
            query,
            rusqlite::params_from_iter(ids.iter().map(|s| s.as_ref())),
            window,
        )?;
        info!("Done.");
        Ok(GeneBook::Cached(r))
    }

    #[allow(dead_code)]
    pub fn inline(filename: &str, window: usize) -> Result<Self> {
        let conn = Connection::open(filename).map_err(|e| errors::DataError::FailedToConnect {
            source: e,
            filename: filename.into(),
        })?;
        Ok(GeneBook::Inline(Mutex::new(conn), window))
    }

    pub fn get(&self, g: &str) -> Result<Gene> {
        match self {
            GeneBook::InMemory(book) | GeneBook::Cached(book) => book
                .get(g)
                .cloned()
                .ok_or_else(|| errors::RuntimeError::IdNotFound(g.to_owned()).into()),
            GeneBook::Inline(conn_mutex, window) => {
                let conn = conn_mutex.lock().expect("MUTEX POISONING");
                let mut query = conn.prepare(
                    "SELECT left_tail_ids, right_tail_ids, ancestral_id, species FROM genomes WHERE id=?",
                )?;
                query
                    .query_row(&[g], |r| {
                        let species = r.get::<_, String>(3)?;

                        let mut left_landscape = Self::parse_landscape(&r.get::<_, String>(0)?);
                        left_landscape.reverse();
                        left_landscape.truncate(*window);
                        left_landscape.reverse();

                        let mut right_landscape = Self::parse_landscape(&r.get::<_, String>(1)?);
                        right_landscape.truncate(*window);

                        rusqlite::Result::Ok(Gene {
                            species,
                            landscape: left_landscape
                                .into_iter()
                                .chain([r.get::<usize, _>(2)?].into_iter())
                                .chain(right_landscape.into_iter())
                                .collect(),
                        })
                    })
                    .with_context(|| "while accessing DB")
            }
        }
    }
}

pub fn write_dist_matrix<'a, T: std::fmt::Display, L: std::fmt::Display>(
    m: &[T],
    ids: &[L],
    mut out: BufWriter<Box<dyn std::io::Write + 'a>>,
) -> Result<()> {
    let n = ids.len();
    assert!(m.len() == n * (n - 1) / 2);

    // 1. Write the number of elements
    writeln!(out, "{}", ids.len())?;

    // 2. Write the matrix itself
    for (i, id) in ids.iter().enumerate() {
        write!(out, "{}", id)?;
        for j in 0..i {
            write!(out, "\t{:.5}", m[i * (i - 1) / 2 + j])?;
        }
        write!(out, "\t{:.5}", 0.)?;
        for j in i + 1..ids.len() {
            write!(out, "\t{:.5}", m[j * (j - 1) / 2 + i])?;
        }
        writeln!(out)?;
    }

    // 3. Flush the output
    Ok(out.flush()?)
}

pub fn read_genefile(filename: &str) -> Result<Vec<Vec<String>>> {
    let mut filecontent = String::new();
    File::open(filename)
        .map_err(|e| errors::FileError::CannotOpen { source: e, filename: filename.into() })?
        .read_to_string(&mut filecontent)?;
    let filecontent = filecontent.trim();
    if filecontent.starts_with('(') && filecontent.ends_with(';') {
        let trees = newick::from_string(&filecontent)?;
        trees
            .into_iter()
            .map(|tree| {
                tree.leaves()
                    .map(|l| {
                        tree[l]
                            .data
                            .name
                            .as_ref()
                            .map(|s| s.to_owned())
                            .with_context(|| format!("nameless leaf found in {:?}", filename))
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect()
    } else {
        Ok(vec![filecontent.split('\n').map(|s| s.to_owned()).collect()])
    }
}
