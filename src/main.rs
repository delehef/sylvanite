use anyhow::*;
use clap::{ArgGroup, Parser, Subcommand};
use log::*;
use std::fs::File;
use std::io::Write;
use std::time::Instant;

use crate::errors::FileError;

mod align;
mod dede;
mod errors;
mod polytomic_tree;
mod sylva;
mod synteny;
mod utils;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Cli {
    #[clap(short, long)]
    verbose: bool,
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
        #[clap(short = 'D', long)]
        database: String,
        #[clap(required = true)]
        infiles: Vec<String>,
        #[clap(short, long)]
        outdir: Option<String>,
        #[clap(short, long)]
        bar: bool,
    },
    BuildTrees {
        #[clap(short = 'D', long)]
        database: String,
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

fn paths2files<S: AsRef<str>>(fs: &[S]) -> Result<Vec<String>> {
    let mut r = Vec::new();
    for f in fs {
        if std::fs::metadata(f.as_ref())
            .map_err(|_| errors::FileError::NotFound(f.as_ref().into()))?
            .is_dir()
        {
            let it = std::fs::read_dir(f.as_ref())?
                .map(|entry| entry.unwrap().path())
                .filter(|p| std::fs::metadata(p).unwrap().is_file())
                .map(|p| p.to_str().unwrap().to_owned())
                .collect::<Vec<_>>()
                .into_iter();
            r.extend(it);
        } else {
            r.push(f.as_ref().to_string())
        }
    }
    Ok(r)
}

fn main() -> Result<()> {
    let args = Cli::parse();
    buche::new()
        .timestamp(buche::Timestamp::Off)
        .verbosity(if args.verbose { 4 } else { 2 })
        .init()
        .unwrap();

    rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    debug!("Using {} threads", rayon::current_num_threads());

    match args.command {
        Commands::Align { infiles, outdir, bar, database } => {
            for f in paths2files(&infiles)?.into_iter() {
                info!("Processing {:?}", f);
                let now = Instant::now();
                let out = synteny::process_file(&f, &database, args.window, &outdir, bar)?;
                debug!("Done in {}s. Result written to {:?}", now.elapsed().as_secs(), out);
            }
            Ok(())
        }
        Commands::BuildDatabase { gffs, emf } => {
            dbg!(&emf);
            dbg!(&gffs);
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
            database,
        } => {
            let batch_name = "pipo";
            let mut timings = if let Some(timings) = timings {
                let mut timings = File::create(&timings).map_err(|source| {
                    FileError::WhileCreating { source, filename: timings.into() }
                })?;
                timings.write_all("file,size,time\n".as_bytes())?;
                Some(timings)
            } else {
                None
            };

            for f in paths2files(&infiles)?.into_iter() {
                let mut input_filename = std::path::PathBuf::from(&f);
                input_filename.set_file_name(
                    input_filename
                        .with_extension("nhx")
                        .file_name()
                        .ok_or_else(|| {
                            errors::FileError::InvalidFilename(format!("{:?}", input_filename))
                        })?
                        .to_str()
                        .ok_or_else(|| {
                            errors::FileError::InvalidFilename(format!("{:?}", input_filename))
                        })?,
                );

                let out_file = std::path::PathBuf::from(if let Some(ref outdir) = outdir {
                    if !std::path::Path::new(outdir).exists() {
                        std::fs::create_dir(outdir).map_err(|source| FileError::WhileCreating {
                            source,
                            filename: outdir.into(),
                        })?;
                    }
                    format!(
                        "{}/{}",
                        outdir,
                        input_filename
                            .file_name()
                            .ok_or_else(|| errors::FileError::InvalidFilename(format!(
                                "{:?}",
                                input_filename
                            )))?
                            .to_str()
                            .ok_or_else(|| errors::FileError::InvalidFilename(format!(
                                "{:?}",
                                input_filename
                            )))?
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
                        &species_tree,
                        &database,
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
