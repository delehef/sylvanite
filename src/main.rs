use rusqlite::Connection;
use std::fs::File;
use std::io;
use std::io::BufWriter;
use std::io::Write;
use std::sync::Mutex;
use std::{
    collections::HashMap,
    io::BufRead,
    path::{Path, PathBuf},
};

use anyhow::*;
use clap::Parser;

use rayon::prelude::*;
mod align;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Settings {
    #[clap(required = true)]
    infiles: Vec<String>,
    #[clap(short, long)]
    outdir: Option<String>,
    #[clap(short, long)]
    database: String,
    #[clap(short, long, default_value_t = 15)]
    window: usize,
    #[clap(short, long, default_value_t = 0)]
    threads: usize,
}

struct Gene {
    gene: String,
    protein: String,
    landscape: Vec<usize>,
}

struct Register {
    genes: HashMap<String, Gene>,
}

fn write_dist_matrix<'a, T: std::fmt::Display, L: std::fmt::Display>(
    m: &[T],
    ids: &[L],
    mut out: BufWriter<Box<dyn std::io::Write + 'a>>,
) -> Result<()> {
    let n = ids.len();
    assert!(m.len() == n * n);

    // 1. Write the number of elements
    writeln!(out, "{}", ids.len())?;

    // 2. Write the matrix itself
    for (i, id) in ids.iter().enumerate() {
        write!(out, "{}", id)?;
        for j in 0..ids.len() {
            write!(out, "\t{:.5}", m[n * i + j])?;
        }
        writeln!(out)?;
    }

    // 3. Flush the output
    Ok(out.flush()?)
}

fn read_db(filename: &str, window: usize) -> Result<Register> {
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
        "SELECT gene, protein, left_tail_ids, right_tail_ids, ancestral_id FROM genomes",
    )?;
    let genes = query
        .query_map([], |r| {
            std::result::Result::Ok((
                r.get::<_, String>(0)?,
                r.get::<_, String>(1)?,
                r.get::<_, String>(2)?,
                r.get::<_, String>(3)?,
                r.get::<_, usize>(4)?,
            ))
        })?
        .collect::<Result<Vec<_>, _>>()?;

    Ok(Register {
        genes: genes
            .into_par_iter()
            .map(|g| {
                let mut left_landscape = parse_landscape(&g.2);
                left_landscape.reverse();
                left_landscape.truncate(window);
                left_landscape.reverse();

                let mut right_landscape = parse_landscape(&g.3);
                right_landscape.truncate(window);

                (
                    g.0.clone(),
                    Gene {
                        gene: g.0,
                        protein: g.1,
                        landscape: left_landscape
                            .into_iter()
                            .chain([g.4].into_iter())
                            .chain(right_landscape.into_iter())
                            .collect(),
                    },
                )
            })
            .collect(),
    })
}

fn process_file(filename: &str, register: &Register, settings: &Settings) -> Result<PathBuf> {
    let mut outfile = if let Some(ref outdir) = settings.outdir {
        PathBuf::from(outdir)
    } else {
        PathBuf::from(Path::new(filename).parent().unwrap())
    };
    outfile.set_file_name(Path::new(filename).with_extension("mat"));
    let out: BufWriter<Box<dyn std::io::Write>> = BufWriter::with_capacity(30_000_000, Box::new(File::create(&outfile)?));

    let genes: Vec<String> = io::BufReader::new(File::open(filename)?)
        .lines()
        .map(|l| l.with_context(|| "while reading genes"))
        .collect::<Result<Vec<_>>>()?;
    let m = Mutex::new(vec![0f32; genes.len().pow(2)]);
    for (i, g1) in genes.iter().enumerate() {
        genes[0..i]
            .par_iter()
            .enumerate()
            .map(|(j, g2)| {
                let g1 = register
                    .genes
                    .get(g1)
                    .with_context(|| format!("`{}` not found in database", g1))?;
                let g2 = register
                    .genes
                    .get(g2)
                    .with_context(|| format!("`{}` not found in database", g2))?;

                let score =
                    align::score_landscape(&g1.landscape, &g2.landscape, &|x, y| x.max(y) as f32);
                m.lock()
                    .map(|mut m| {
                        m[i * genes.len() + j] = score;
                        m[j * genes.len() + i] = score;
                    })
                    .expect("MUTEX POISONING");
                Ok(())
            })
            .collect::<Result<Vec<_>>>()?;
    }
    write_dist_matrix(&m.into_inner().expect("BROKEN MUTEX"), &genes, out)?;
    Ok(outfile)
}

fn main() -> Result<()> {
    let args = Settings::parse();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();
    println!("Using {} threads", rayon::current_num_threads());

    let register = read_db(&args.database, args.window)?;

    for f in args.infiles.iter() {
        println!("Processing {}", f);
        let out = process_file(&f, &register, &args)?;
        println!("Result written to {:?}\n", out);
    }

    Ok(())
}
