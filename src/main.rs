use indicatif::ProgressBar;
use log::*;
use rusqlite::Connection;
use std::fs::File;
use std::io::BufWriter;
use std::io::Read;
use std::io::Write;
use std::sync::Mutex;
use std::time::Instant;
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

use anyhow::*;
use clap::Parser;

use rayon::prelude::*;
mod align;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Settings {
    #[clap(short, long)]
    verbose: bool,
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

fn read_db(filename: &str, window: usize) -> Result<Register> {
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

    info!("Done.");
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
                    g.1.clone(),
                    Gene {
                        // gene: g.0,
                        // protein: g.1,
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

fn read_genefile(filename: &str) -> Result<Vec<String>> {
    let mut filecontent = String::new();
    File::open(filename)?.read_to_string(&mut &mut filecontent)?;
    let filecontent = filecontent.trim();
    if filecontent.starts_with("(") && filecontent.ends_with(";") {
        let tree = newick::from_string(&filecontent)?;
        tree.leaves()
            .map(|l| {
                tree[l]
                    .data
                    .name
                    .as_ref()
                    .map(|s| s.to_owned())
                    .with_context(|| format!("nameless leaf found in `{}`", filename))
            })
            .collect::<Result<Vec<_>>>()
    } else {
        Ok(filecontent.split("\n").map(|s| s.to_owned()).collect())
    }
}

fn process_file(filename: &str, register: &Register, settings: &Settings) -> Result<PathBuf> {
    let mut outfile = if let Some(ref outdir) = settings.outdir {
        PathBuf::from(outdir)
    } else {
        PathBuf::from(Path::new(filename).parent().unwrap())
    };
    outfile.set_file_name(Path::new(filename).with_extension("dist"));
    let out: BufWriter<Box<dyn std::io::Write>> =
        BufWriter::with_capacity(30_000_000, Box::new(File::create(&outfile)?));
    let genes = read_genefile(filename)?;

    let n = genes.len();
    let m = Mutex::new(vec![0f32; n * (n - 1) / 2]);
    let bar = if settings.verbose && atty::is(atty::Stream::Stdout) {
        Some(ProgressBar::new(genes.len() as u64))
    } else {
        None
    };
    for (i, g1) in genes.iter().enumerate() {
        bar.as_ref().map(|b| b.inc(1));
        genes[0..i].par_iter().enumerate().try_for_each(|(j, g2)| {
            let gg1 = register
                .genes
                .get(g1)
                .with_context(|| format!("`{}` not found in database", g1))?;
            let gg2 = register
                .genes
                .get(g2)
                .with_context(|| format!("`{}` not found in database", g2))?;

            let score =
                align::score_landscape(&gg1.landscape, &gg2.landscape, &|x, y| x.max(y) as f32);
            m.lock()
                .map(|mut m| {
                    m[i * (i - 1) / 2 + j] = score;
                })
                .expect("MUTEX POISONING");
            Ok(())
        })?;
    }
    bar.map(|b| b.finish_and_clear());
    write_dist_matrix(&m.into_inner().expect("BROKEN MUTEX"), &genes, out)?;
    Ok(outfile)
}

fn main() -> Result<()> {
    let args = Settings::parse();
    stderrlog::new()
        .timestamp(stderrlog::Timestamp::Off)
        .verbosity(if args.verbose { 3 } else { 2 })
        .show_level(false)
        .init()
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();
    debug!("Using {} threads", rayon::current_num_threads());

    let register = read_db(&args.database, args.window)?;

    for f in args.infiles.iter() {
        info!("Processing {}", f);
        let now = Instant::now();
        let out = process_file(&f, &register, &args)?;
        debug!(
            "Done in {}s. Result written to {:?}\n",
            now.elapsed().as_secs(),
            out
        );
    }

    Ok(())
}
