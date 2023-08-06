use crate::align;
use crate::errors::{FileError, RuntimeError};
use crate::utils::*;
use anyhow::*;
use indicatif::ProgressBar;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use std::path::*;
use std::sync::Mutex;
use syntesuite::genebook::GeneBook;

pub fn process_file(
    filename: &str,
    db_file: &str,
    window: usize,
    outdir: &Option<String>,
    do_bar: bool,
) -> Result<PathBuf> {
    let outdir = if let Some(ref outdir) = outdir {
        PathBuf::from(outdir)
    } else {
        PathBuf::from(Path::new(filename).parent().unwrap())
    };
    if !outdir.exists() {
        std::fs::create_dir_all(&outdir)
            .map_err(|source| FileError::WhileCreating { source, filename: String::new() })?;
    }
    let outfile = outdir.join(Path::new(filename).with_extension("dist").file_name().unwrap());

    let ids = crate::utils::read_genefile(filename)?;
    if ids.is_empty() {
        bail!("{} is empty", filename)
    }
    let book = GeneBook::cached(db_file, window, "id", &ids)?;

    let n = ids.len();
    let m = Mutex::new(vec![0f32; n * (n - 1) / 2]);
    let bar = if do_bar && n > 5000 && atty::is(atty::Stream::Stdout) {
        Some(ProgressBar::new(ids.len() as u64))
    } else {
        None
    };
    for (i, g1) in ids.iter().enumerate() {
        if let Some(b) = bar.as_ref() {
            b.inc(1)
        }
        ids[0..i].par_iter().enumerate().try_for_each(|(j, g2)| {
            let gg1 = book.get(g1).map_err(|_| RuntimeError::IdNotFound(g1.into()))?;
            let gg2 = book.get(g2).map_err(|_| RuntimeError::IdNotFound(g2.into()))?;

            let score = align::score_landscape(
                &gg1.landscape().collect::<Vec<_>>(),
                &gg2.landscape().collect::<Vec<_>>(),
                &|x, y| x.max(y) as f32,
            );
            m.lock()
                .map(|mut m| {
                    m[i * (i - 1) / 2 + j] = score;
                })
                .expect("MUTEX POISONING");
            Ok(())
        })?;
    }
    if let Some(b) = bar {
        b.finish_and_clear()
    }

    write_dist_matrix(
        &m.into_inner().expect("BROKEN MUTEX"),
        &ids,
        BufWriter::with_capacity(30_000_000, Box::new(File::create(&outfile)?)),
    )?;
    Ok(outfile)
}
