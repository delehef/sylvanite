use anyhow::*;
use clap::{Parser, Subcommand};
use log::*;
use std::fs::File;
use std::io::prelude::*;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;
use utils::*;

use rayon::prelude::*;
mod align;
mod dede;
mod polytomic_tree;
mod sylva;
mod synteny;
mod utils;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Cli {
    #[clap(short, long)]
    verbose: bool,
    #[clap(short, long)]
    database: String,
    #[clap(long)]
    cache_db: bool,
    #[clap(short, long, default_value_t = 15)]
    window: usize,
    #[clap(short, long, default_value_t = 0)]
    threads: usize,

    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    Align {
        #[clap(required = true)]
        infiles: Vec<String>,
        #[clap(short, long)]
        outdir: Option<String>,
        #[clap(short, long)]
        bar: bool,
    },
    BuildTrees {
        #[clap(short = 'S', long, required = true)]
        species_tree: String,
        #[clap(short, long, required = true)]
        syntenies: String,
        #[clap(short, long, required = true)]
        divergences: String,
        #[clap(required = true)]
        infiles: Vec<String>,
        #[clap(short, long)]
        outdir: Option<String>,
        #[clap(long)]
        timings: Option<String>,
    },
}

fn main() -> Result<()> {
    let args = Cli::parse();
    stderrlog::new()
        .timestamp(stderrlog::Timestamp::Off)
        .verbosity(if args.verbose { 3 } else { 3 })
        .show_level(false)
        .init()
        .unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();
    debug!("Using {} threads", rayon::current_num_threads());

    let gene_book = if args.cache_db {
        GeneBook::cached(&args.database, args.window)
    } else {
        GeneBook::inline(&args.database, args.window)
    }?;

    match args.command {
        Commands::Align {
            infiles,
            outdir,
            bar,
        } => {
            for f in infiles.iter() {
                info!("Processing {:?}", f);
                let now = Instant::now();
                let out = synteny::process_file(&f, &gene_book, &outdir, bar)?;
                debug!(
                    "Done in {}s. Result written to {:?}",
                    now.elapsed().as_secs(),
                    out
                );
            }
            Ok(())
        }
        Commands::BuildTrees {
            species_tree,
            syntenies,
            divergences,
            infiles,
            outdir,
            timings,
        } => {
            let batch_name = "pipo";
            let mut timings = if let Some(timings) = timings {
                let mut timings = File::create(&timings)
                    .with_context(|| format!("while creating {}", &timings))?;
                timings.write_all("file,time".as_bytes())?;
                Some(timings)
            } else {
                None
            };

            let todos = infiles
                .into_iter()
                .flat_map(|f| {
                    if std::fs::metadata(&f).unwrap().is_dir() {
                        std::fs::read_dir(&f)
                            .unwrap()
                            .map(|entry| entry.unwrap().path())
                            .filter(|p| std::fs::metadata(p).unwrap().is_file())
                            .map(|p| p.to_str().unwrap().to_owned())
                            .collect::<Vec<_>>()
                            .into_iter()
                    } else {
                        vec![f].into_iter()
                    }
                })
                .collect::<Vec<String>>();

            for f in todos.iter() {
                let mut input_filename = std::path::PathBuf::from(&f);
                input_filename.set_file_name(format!(
                    "sylvanite_{}",
                    input_filename
                        .file_name()
                        .ok_or(anyhow!("invalid filename found"))?
                        .to_str()
                        .ok_or(anyhow!("invalid filename found"))?
                ));

                let out_file = std::path::PathBuf::from(if let Some(ref outdir) = outdir {
                    format!(
                        "{}/{}",
                        outdir,
                        input_filename.file_name().unwrap().to_str().unwrap()
                    )
                } else {
                    input_filename.to_str().unwrap().to_owned()
                });

                if out_file.exists() {
                    println!("{} already exists; skipping", out_file.display());
                } else {
                    println!("Creating {}...", out_file.display());
                    let now = Instant::now();
                    let tree = sylva::do_file(
                        &f,
                        &batch_name,
                        &gene_book,
                        &species_tree,
                        &syntenies,
                        &divergences,
                    )?;
                    info!("Done in {:.2}s.", now.elapsed().as_secs_f32());
                    if let Some(ref mut timings) = timings {
                        timings.write_all(
                            format!(
                                "{:#?},{}",
                                Path::new(f).file_name().unwrap(),
                                now.elapsed().as_secs_f32()
                            )
                            .as_bytes(),
                        )?;
                    }
                    std::fs::File::create(out_file)?.write_all(tree.as_bytes())?
                }
            }

            Ok(())
        }
    }
}
