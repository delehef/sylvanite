use anyhow::*;
use clap::{Parser, Subcommand};
use log::*;
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
        out: Option<String>,
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
    }  else {
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
            out,
        } => {
            let batch_name = "pipo";
            for f in infiles {
                let mut out_file = std::path::PathBuf::from(&f);
                out_file.set_file_name(format!("sylvanite_{:?}", out_file.file_name()));
                let tree = sylva::do_file(
                    &f,
                    &batch_name,
                    &gene_book,
                    &species_tree,
                    &syntenies,
                    &divergences,
                )?;
                std::fs::File::create(out_file)?.write_all(tree.as_bytes())?
            }

            Ok(())
        }
    }
}
