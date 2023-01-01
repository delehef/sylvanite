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
        #[clap(long)]
        no_overwrite: bool,
    },
}

fn paths2files<S: AsRef<str>>(fs: &[S]) -> Vec<String> {
    fs.into_iter()
        .flat_map(|f| {
            if std::fs::metadata(f.as_ref()).unwrap().is_dir() {
                std::fs::read_dir(f.as_ref())
                    .unwrap()
                    .map(|entry| entry.unwrap().path())
                    .filter(|p| std::fs::metadata(p).unwrap().is_file())
                    .map(|p| p.to_str().unwrap().to_owned())
                    .collect::<Vec<_>>()
                    .into_iter()
            } else {
                vec![f.as_ref().to_owned()].into_iter()
            }
        })
        .collect::<Vec<String>>()
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

    let register = if args.cache_db {
        Some(GeneBook::in_memory(&args.database, args.window)?)
    } else {
        None
    };

    match args.command {
        Commands::Align {
            infiles,
            outdir,
            bar,
        } => {
            for f in paths2files(&infiles).into_iter() {
                info!("Processing {:?}", f);
                let now = Instant::now();
                let out = synteny::process_file(&f, &args.database, args.window, &outdir, bar)?;
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
            no_overwrite,
        } => {
            let batch_name = "pipo";
            let mut timings = if let Some(timings) = timings {
                let mut timings = File::create(&timings)
                    .with_context(|| format!("while creating {}", &timings))?;
                timings.write_all("file,size,time\n".as_bytes())?;
                Some(timings)
            } else {
                None
            };

            for f in paths2files(&infiles).into_iter() {
                let mut input_filename = std::path::PathBuf::from(&f);
                input_filename.set_file_name(format!(
                    "sylvanite_{}",
                    input_filename
                        .file_name()
                        .ok_or_else(|| anyhow!("invalid filename found"))?
                        .to_str()
                        .ok_or_else(|| anyhow!("invalid filename found"))?
                ));

                let out_file = std::path::PathBuf::from(if let Some(ref outdir) = outdir {
                    if !std::path::Path::new(outdir).exists() {
                        std::fs::create_dir(outdir)
                            .with_context(|| anyhow!("while creating `{}`", outdir))?;
                    }
                    format!(
                        "{}/{}",
                        outdir,
                        input_filename.file_name().unwrap().to_str().unwrap()
                    )
                } else {
                    input_filename.to_str().unwrap().to_owned()
                });

                if out_file.exists() && no_overwrite {
                    info!("{} already exists; skipping", out_file.display());
                } else {
                    let tree = sylva::do_file(
                        &f,
                        batch_name,
                        register.as_ref(),
                        &species_tree,
                        &args.database,
                        args.window,
                        &syntenies,
                        &divergences,
                        &mut timings,
                    )?;
                    std::fs::File::create(out_file)?.write_all(tree.as_bytes())?
                }
            }

            Ok(())
        }
    }
}
