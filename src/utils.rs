use anyhow::*;
use newick::Newick;
use std::fs::File;
use std::io::BufWriter;
use std::io::Read;
use std::io::Write;

use crate::errors;

pub fn write_dist_matrix<'a, T: std::fmt::Display, L: std::fmt::Display>(
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

pub fn read_genefile(filename: &str) -> Result<Vec<Vec<String>>> {
    let mut filecontent = String::new();
    File::open(filename)
        .map_err(|e| errors::FileError::CannotOpen { source: e, filename: filename.into() })?
        .read_to_string(&mut filecontent)?;
    let filecontent = filecontent.trim();
    if filecontent.starts_with('(') && filecontent.ends_with(';') {
        let trees = newick::from_string(filecontent)?;
        trees
            .into_iter()
            .map(|tree| {
                tree.leaves()
                    .map(|l| {
                        tree.name(l)
                            .map(|s| s.to_owned())
                            .with_context(|| format!("nameless leaf found in {:?}", filename))
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .collect()
    } else {
        Ok(vec![filecontent.split('\n').map(|s| s.to_owned()).collect()])
    }
}
