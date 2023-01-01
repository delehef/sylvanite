use crate::align;
use crate::utils::*;
use crate::Commands;
use anyhow::*;
use indicatif::ProgressBar;
use log::*;
use rayon::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use std::path::*;
use std::sync::Mutex;

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
        std::fs::create_dir(&outdir).with_context(|| anyhow!("while creating `{:?}`", outdir))?;
    }
    let outfile = outdir.join(
        Path::new(filename)
            .with_extension("dist")
            .file_name()
            .unwrap(),
    );

    let gene_families = crate::utils::read_genefile(filename)?;
    if gene_families.is_empty() || gene_families.len() > 1 {
        bail!("{} should contain a single family", filename)
    }
    let ids = &gene_families[0];
    let book = GeneBook::cached(db_file, window, ids)?;

    let n = ids.len();
    let m = Mutex::new(vec![0f32; n * (n - 1) / 2]);
    let bar = if do_bar && n > 5000 && atty::is(atty::Stream::Stdout) {
        Some(ProgressBar::new(ids.len() as u64))
    } else {
        None
    };
    for (i, g1) in ids.iter().enumerate() {
        bar.as_ref().map(|b| b.inc(1));
        ids[0..i].par_iter().enumerate().try_for_each(|(j, g2)| {
            let gg1 = book
                .get(g1)
                .with_context(|| format!("`{}` not found in database", g1))?;
            let gg2 = book
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
        &ids,
        BufWriter::with_capacity(30_000_000, Box::new(File::create(&outfile)?)),
    )?;
    Ok(outfile)
}
