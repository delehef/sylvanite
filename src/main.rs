use anyhow::*;
use clap::Parser;
use dede::*;
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
use utils::*;

use rayon::prelude::*;
mod align;
mod dede;
mod polytomic_tree;
mod sylva;
mod utils;

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

fn read_genefile(filename: &Path) -> Result<Vec<String>> {
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
                    .with_context(|| format!("nameless leaf found in {:?}", filename))
            })
            .collect::<Result<Vec<_>>>()
    } else {
        Ok(filecontent.split("\n").map(|s| s.to_owned()).collect())
    }
}

fn process_file(filename: &Path, register: &GeneBook, settings: &Settings) -> Result<PathBuf> {
    let outdir = if let Some(ref outdir) = settings.outdir {
        PathBuf::from(outdir)
    } else {
        PathBuf::from(Path::new(filename).parent().unwrap())
    };
    let outfile = outdir.join(
        Path::new(filename)
            .with_extension("dist")
            .file_name()
            .unwrap(),
    );

    let genes = read_genefile(filename)?;

    let n = genes.len();
    let m = Mutex::new(vec![0f32; n * (n - 1) / 2]);
    let bar = if settings.verbose && n > 5000 && atty::is(atty::Stream::Stdout) {
        Some(ProgressBar::new(genes.len() as u64))
    } else {
        None
    };
    for (i, g1) in genes.iter().enumerate() {
        bar.as_ref().map(|b| b.inc(1));
        genes[0..i].par_iter().enumerate().try_for_each(|(j, g2)| {
            let gg1 = register
                .get(g1)
                .with_context(|| format!("`{}` not found in database", g1))?;
            let gg2 = register
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

    write_dist_matrix(
        &m.into_inner().expect("BROKEN MUTEX"),
        &genes,
        BufWriter::with_capacity(30_000_000, Box::new(File::create(&outfile)?)),
    )?;
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

    // let lines =
    //     std::fs::read_to_string("/users/ldog/delehell/duplications/data/SuperTrees.nhx")
    //         .unwrap();
    // let lines = lines.split('\n').collect::<Vec<_>>();
    let register = read_db(&args.database, args.window).unwrap();
    // sylva::do_family(&lines[720], 720, "pipo", &register);
    // return Ok(());

    for p in args.infiles.iter() {
        let files = if Path::new(&p).is_dir() {
            std::fs::read_dir(p)?
                // .with_context(|| format!("while opening directory `{}`", p))?
                .map(|x| Ok(x?.path()))
                .collect::<Result<Vec<_>>>()?
        } else {
            vec![PathBuf::from(p)]
        };
        for f in files {
            info!("Processing {:?}", f);
            let now = Instant::now();
            let out = process_file(&f, &register, &args)?;
            debug!(
                "Done in {}s. Result written to {:?}",
                now.elapsed().as_secs(),
                out
            );
        }
    }

    Ok(())
}
