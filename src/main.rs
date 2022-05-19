use anyhow::*;
use clap::{Parser, Subcommand};
use log::*;
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
    Build {
        #[clap(required = true)]
        infiles: Vec<String>,
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

    let register = read_db(&args.database, args.window).unwrap();

    match args.command {
        Commands::Align {
            infiles,
            outdir,
            bar,
        } => {
            for f in infiles.iter() {
                info!("Processing {:?}", f);
                let now = Instant::now();
                let out = synteny::process_file(&f, &register, &outdir, bar)?;
                debug!(
                    "Done in {}s. Result written to {:?}",
                    now.elapsed().as_secs(),
                    out
                );
            }
            Ok(())
        }
        Commands::Build { infiles } => {
            for f in infiles {
                sylva::do_file(&f, &register);
            }

            Ok(())
        }
    }
}
