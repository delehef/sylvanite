use anyhow::*;
use newick::Newick;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::hash::Hash;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Read;
use std::io::Write;

use crate::dede::Matrix;
use crate::dede::VecMatrix;
use crate::errors;
use crate::errors::MatrixParseError;

pub(crate) fn parse_dist_matrix<S: AsRef<str>>(
    filename: &str,
    ids: &[S],
) -> Result<VecMatrix<f32>> {
    let n = ids.len();
    let mut r = VecMatrix::<f32>::with_elem(n, n, 0.);
    let mut tmp = VecMatrix::<f32>::with_elem(n, n, 0.);
    let mut r2g = HashMap::<String, usize>::new();

    let mut lines = BufReader::new(File::open(filename)?).lines();
    let nn = lines.next().ok_or(MatrixParseError::SizeMissing)??.parse::<usize>()?;
    ensure!(
        nn == n,
        format!("{} does not have the expected size (expected {}, found {})", filename, n, nn)
    );
    for (i, l) in lines.enumerate() {
        let l = l?;
        let mut s = l.split('\t');
        let gene = s.next().ok_or(MatrixParseError::ErroneousLine)?;
        r2g.insert(gene.into(), i);
        for (j, x) in s.enumerate() {
            tmp[(i, j)] = x.parse::<f32>()?;
        }
    }

    for i in 0..n {
        let gi = ids[i].as_ref();
        for j in 0..n {
            let gj = ids[j].as_ref();
            r[(i, j)] = tmp[(
                *r2g.get(gi).with_context(|| format!("`{}` not found in {}", gi, filename))?,
                *r2g.get(gj).with_context(|| format!("`{}` not found in {}", gj, filename))?,
            )];
        }
    }

    Ok(r)
}

pub(crate) fn parse_dist_matrix_from_sim<S: AsRef<str>>(
    filename: &str,
    ids: &[S],
) -> Result<VecMatrix<f32>> {
    let n = ids.len();
    let mut r = VecMatrix::<f32>::with_elem(n, n, 0.);
    let mut tmp = VecMatrix::<f32>::with_elem(n, n, 0.);
    let mut r2g = HashMap::<String, usize>::new();

    let mut lines = BufReader::new(File::open(filename)?).lines();
    let nn = lines.next().ok_or(MatrixParseError::SizeMissing)??.parse::<usize>()?;
    ensure!(nn == n, "{} does not have the expected size (expected {}, found {})", filename, n, nn);
    for (i, l) in lines.enumerate() {
        let l = l?;
        let mut s = l.split('\t');
        let gene = s.next().ok_or(MatrixParseError::ErroneousLine)?;
        r2g.insert(gene.into(), i);
        for (j, x) in s.enumerate() {
            let similarity = x.parse::<f32>()?;
            ensure!(similarity <= 1., "{} is not a similarity score (0 <= x <= 1)", similarity);
            tmp[(i, j)] = 1. - similarity;
        }
    }

    for i in 0..n {
        let gi = ids[i].as_ref();
        for j in 0..n {
            let gj = ids[j].as_ref();
            r[(i, j)] = tmp[(
                *r2g.get(gi).with_context(|| format!("`{}` not found in {}", gi, filename))?,
                *r2g.get(gj).with_context(|| format!("`{}` not found in {}", gj, filename))?,
            )];
        }
    }

    Ok(r)
}

pub(crate) fn jaccard<T: Eq + Hash>(a: &HashSet<T>, b: &HashSet<T>) -> f32 {
    a.intersection(b).count() as f32 / a.union(b).count() as f32
}

pub(crate) fn write_dist_matrix<'a, T: std::fmt::Display, L: std::fmt::Display>(
    m: &[T],
    ids: &[L],
    mut out: BufWriter<Box<dyn std::io::Write + 'a>>,
) -> Result<()> {
    let n = ids.len();
    assert!(m.len() == n * (n - 1) / 2, "{} vs {}", m.len(), n * (n - 1) / 2);

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

pub(crate) fn write_matrix<'a, T: std::fmt::Display, L: std::fmt::Display>(
    m: &impl Matrix<T>,
    ids: &[L],
    mut out: BufWriter<Box<dyn std::io::Write + 'a>>,
) -> Result<()> {
    assert!(m.nrows() == m.ncols(), "m is not squared");
    assert!(m.nrows() == ids.len(), "{} rows vs. {} ids", m.nrows(), ids.len());

    // 1. Write the number of elements
    writeln!(out, "{}", ids.len())?;

    // 2. Write the matrix itself
    for (i, id) in ids.iter().enumerate() {
        write!(out, "{}", id)?;
        for j in 0..ids.len() {
            write!(out, "\t{:.5}", m[(i, j)])?;
        }
        writeln!(out)?;
    }

    // 3. Flush the output
    Ok(out.flush()?)
}

pub(crate) fn read_genefile(filename: &str) -> Result<Vec<String>> {
    let mut filecontent = String::new();
    File::open(filename)
        .map_err(|e| errors::FileError::CannotOpen { source: e, filename: filename.into() })?
        .read_to_string(&mut filecontent)?;
    let filecontent = filecontent.trim();
    if filecontent.starts_with('(') && filecontent.ends_with(';') {
        let tree = newick::one_from_string(&filecontent)?;
        tree.leaves()
            .map(|l| {
                tree.name(l)
                    .map(|s| s.to_owned())
                    .with_context(|| format!("nameless leaf found in {:?}", filename))
            })
            .collect::<Result<Vec<_>>>()
    } else {
        Ok(filecontent.split('\n').map(|s| s.to_owned()).collect())
    }
}
