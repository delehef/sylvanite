use anyhow::*;
use log::*;
use rusqlite::Connection;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Read;
use std::io::Write;
use std::sync::{Arc, Mutex};

const REQUEST: &str = "SELECT gene, protein, left_tail_ids, right_tail_ids, ancestral_id, species FROM genomes WHERE protein=?";

pub enum GeneBook {
    Cached(HashMap<String, Gene>),
    Inline(Mutex<Connection>, usize),
}

#[derive(Clone, Default)]
pub struct Gene {
    pub gene: String,
    pub species: String,
    pub landscape: Vec<usize>,
}

impl GeneBook {
    fn parse_landscape(landscape: &str) -> Vec<usize> {
        if landscape.is_empty() {
            Vec::new()
        } else {
            landscape
                .split('.')
                .map(|x| x.parse::<usize>().unwrap())
                .collect::<Vec<_>>()
        }
    }

    pub fn cached(filename: &str, window: usize) -> Result<Self> {
        info!("Parsing the database...");

        let conn = Connection::open(filename)
            .with_context(|| format!("while connecting to {}", filename))?;
        let mut query = conn.prepare(
        "SELECT gene, protein, left_tail_ids, right_tail_ids, ancestral_id, species FROM genomes",
    )?;
        let genes = query
            .query_map([], |r| {
                std::result::Result::Ok((
                    r.get::<_, String>(0)?,
                    r.get::<_, String>(1)?,
                    r.get::<_, String>(2)?,
                    r.get::<_, String>(3)?,
                    r.get::<_, usize>(4)?,
                    r.get::<_, String>(5)?,
                ))
            })?
            .collect::<Result<Vec<_>, _>>()?;

        let r = genes
            .into_iter()
            .map(|g| {
                let mut left_landscape = Self::parse_landscape(&g.2);
                left_landscape.reverse();
                left_landscape.truncate(window);
                left_landscape.reverse();

                let mut right_landscape = Self::parse_landscape(&g.3);
                right_landscape.truncate(window);

                (
                    g.1.clone(),
                    Gene {
                        gene: g.0,
                        species: g.5,
                        landscape: left_landscape
                            .into_iter()
                            .chain([g.4].into_iter())
                            .chain(right_landscape.into_iter())
                            .collect(),
                    },
                )
            })
            .collect();
        info!("Done.");
        Ok(GeneBook::Cached(r))
    }

    pub fn inline(filename: &str, window: usize) -> Result<Self> {
        let conn = Connection::open(filename)
            .with_context(|| format!("while connecting to {}", filename))?;
        Ok(GeneBook::Inline(Mutex::new(conn), window))
    }

    pub fn get(&self, g: &str) -> Result<Gene> {
        match self {
            GeneBook::Cached(book) => book
                .get(g)
                .map(|g| g.clone())
                .ok_or(anyhow!("key not found")),
            GeneBook::Inline(conn_mutex, window) => {
                let mut conn = conn_mutex.lock().expect("MUTEX POISONING");
                let mut query = conn.prepare(
                    "SELECT gene, protein, left_tail_ids, right_tail_ids, ancestral_id, species FROM genomes",
                )?;
                query
                    .query_row(&[g], |r| {
                        let gene = r.get::<_, String>(0)?;
                        let species = r.get::<_, String>(5)?;

                        let mut left_landscape = Self::parse_landscape(&r.get::<_, String>(2)?);
                        left_landscape.reverse();
                        left_landscape.truncate(*window);
                        left_landscape.reverse();

                        let mut right_landscape = Self::parse_landscape(&r.get::<_, String>(3)?);
                        right_landscape.truncate(*window);

                        rusqlite::Result::Ok(Gene {
                            gene: gene,
                            species: species,
                            landscape: left_landscape
                            .into_iter()
                            .chain([r.get::<usize, _>(4)?].into_iter())
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
