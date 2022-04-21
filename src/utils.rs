use anyhow::*;
use log::*;
use rayon::prelude::*;
use rusqlite::Connection;
use std::collections::HashMap;
pub type GeneBook = HashMap<String, Gene>;
pub struct Gene {
    pub gene: String,
    pub species: String,
    pub landscape: Vec<usize>,
}

pub fn read_db(filename: &str, window: usize) -> Result<GeneBook> {
    info!("Parsing the database...");
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

    let conn = Connection::open(filename)?;
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

    info!("Done.");
    Ok(genes
        .into_par_iter()
        .map(|g| {
            let mut left_landscape = parse_landscape(&g.2);
            left_landscape.reverse();
            left_landscape.truncate(window);
            left_landscape.reverse();

            let mut right_landscape = parse_landscape(&g.3);
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
        .collect())
}
