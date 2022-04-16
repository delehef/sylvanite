use std::cmp::*;
pub const MIN_INFORMATIVE_SYNTENY: usize = 0;

fn align_sw<T1, T2, V>(s1: T1, s2: T2, normalizer: &dyn Fn(usize, usize) -> f32) -> f32
where
    T1: AsRef<[V]>,
    T2: AsRef<[V]>,
    V: PartialEq,
{
    const MATCH: i16 = 1;
    const MISMATCH: i16 = -1;
    const GAP: i16 = 0;

    let s1 = s1.as_ref();
    let s2 = s2.as_ref();

    let w = s1.len() + 1;
    let h = s2.len() + 1;

    let mut m = Vec::<i16>::with_capacity(w * h);
    unsafe { m.set_len(m.capacity()) }
    let mut max_score = 0;
    for i in 0..w {
        m[i] = 0;
    }
    for i in 0..h {
        m[i * w] = 0;
    }

    for i in 1..h {
        for j in 1..w {
            let diag_score = if s2[i - 1] == s1[j - 1] {
                m[(i - 1) * w + (j - 1)] + MATCH
            } else {
                m[(i - 1) * w + (j - 1)] + MISMATCH
            };
            let up_score = m[(i - 1) * w + j] + GAP;
            let left_score = m[i * w + j - 1] + GAP;
            let score = max(max(diag_score, up_score), left_score);

            m[i * w + j] = score;
            if score > max_score {
                max_score = score;
            }
        }
    }

    if max_score >= 1 {
        max_score -= 1;
    }

    max_score as f32 / normalizer(s1.len(), s2.len())
}

type Path = Vec<(usize, Option<usize>)>;
type Weights = Vec<i16>;
fn score_path(path: &Path, w1: &Weights, w2: &Weights) -> f32 {
    const GAP_PENALTY: i16 = 0;
    const REVERSE_PENALTY: f32 = 0.5;

    let mut score = path
        .iter()
        .map(|p| (w1[p.0]).min(p.1.map(|j| w2[j]).unwrap_or(GAP_PENALTY)) as f32)
        .sum::<f32>();

    if path.len() > 2 {
        let mut dir =
            (path[1].1.unwrap_or_default() as i32 - path[0].1.unwrap_or_default() as i32).signum();
        for j in 1..path.len() - 1 {
            if let Some(jj) = path[j + 1].1 {
                if let Some(j) = path[j].1 {
                    let new_dir = (jj as i32 - j as i32).signum();
                    if new_dir != dir {
                        score -= REVERSE_PENALTY;
                        dir = new_dir;
                    }
                }
            }
        }
    }

    if score >= 1. {
        score - 1.
    } else {
        score
    }
}

fn thread<T1, T2, V>(
    i: usize,
    paths: &mut Vec<Path>,
    matches: &[Vec<usize>],
    s1: T1,
    s2: T2,
    w1: &Weights,
    w2: &Weights,
) where
    T1: AsRef<[V]>,
    T2: AsRef<[V]>,
    V: PartialEq,
{
    let s1 = s1.as_ref();
    let s2 = s2.as_ref();
    if i >= s1.len() {
        return;
    }

    if matches[i].is_empty() {
        paths.iter_mut().for_each(|path| path.push((i, None)));
    } else {
        let matches = &matches[i];
        let m = matches[0];
        paths.iter_mut().for_each(|path| {
            if path.iter().any(|p| p.1.is_some() && p.1.unwrap() == m) {
                path.push((i, None))
            } else {
                path.push((i, Some(m)))
            }
        });

        if matches.len() > 1 {
            let old_paths = paths.len();
            for &m in matches[1..].iter() {
                for k in 0..old_paths {
                    if !paths[k].iter().any(|p| p.1.is_some() && p.1.unwrap() == m) {
                        let mut new_path = paths[k].clone();
                        new_path.last_mut().map(|l| *l = (i, Some(m)));
                        paths.push(new_path);
                    }
                }
            }
        }
    }

    thread(i + 1, paths, matches, s1, s2, w1, w2)
}

fn compress<T, V>(s: T) -> (Vec<V>, Weights)
where
    T: AsRef<[V]>,
    V: PartialEq + Copy,
{
    let s = s.as_ref();
    let mut r = Vec::with_capacity(s.len());
    let mut w = Vec::with_capacity(s.len());
    r.push(s[0]);
    w.push(1);

    for g in s.iter().skip(1) {
        if r.last().unwrap() == g {
            w.last_mut().map(|x| *x += 1);
        } else {
            r.push(*g);
            w.push(1);
        }
    }

    (r, w)
}

fn align_threaded<T1, T2, V>(s1_: T1, s2_: T2, normalizer: &dyn Fn(usize, usize) -> f32) -> f32
where
    T1: AsRef<[V]>,
    T2: AsRef<[V]>,
    V: PartialEq + Copy + std::fmt::Debug,
{
    let (s1, w1) = compress(&s1_);
    let (s2, w2) = compress(&s2_);

    let matches = s1
        .iter()
        .map(|g| {
            s2.iter()
                .enumerate()
                .filter_map(|(j, h)| if h == g { Some(j) } else { None })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let pp: usize = matches
        .iter()
        .filter_map(|m| if m.is_empty() { None } else { Some(m.len()) })
        .product();
    if pp > 10000 {
        return score_landscape_nw(s1_, s2_, normalizer);
    }

    let start = (0..s1.len()).find(|&i| !matches[i].is_empty()).unwrap();
    let mut paths = matches[start]
        .iter()
        .map(|&j| vec![(start, Some(j))])
        .collect::<Vec<_>>();
    thread(start + 1, &mut paths, &matches, &s1, &s2, &w1, &w2);
    (paths
        .iter()
        .map(|p| (2. * score_path(p, &w1, &w2)) as i32)
        .max()
        .unwrap() as f32
        / 2.0)
        / normalizer(s1_.as_ref().len(), s2_.as_ref().len())
}

pub fn score_landscape_nw<T1, T2, V>(
    s1: T1,
    s2: T2,
    normalizer: &dyn Fn(usize, usize) -> f32,
) -> f32
where
    T1: AsRef<[V]>,
    T2: AsRef<[V]>,
    V: PartialEq + Copy,
{
    if s1.as_ref().len() <= MIN_INFORMATIVE_SYNTENY || s2.as_ref().len() < MIN_INFORMATIVE_SYNTENY {
        0.0
    } else {
        let dir = align_sw(&s1, &s2, normalizer);
        let rev = align_sw(
            s1,
            s2.as_ref().iter().copied().rev().collect::<Vec<_>>(),
            normalizer,
        );
        dir.max(rev)
    }
}

pub fn score_landscape<T1, T2, V>(s1: T1, s2: T2, normalizer: &dyn Fn(usize, usize) -> f32) -> f32
where
    T1: AsRef<[V]>,
    T2: AsRef<[V]>,
    V: PartialEq + Copy + std::fmt::Debug,
{
    if s1.as_ref().len() <= MIN_INFORMATIVE_SYNTENY || s2.as_ref().len() < MIN_INFORMATIVE_SYNTENY {
        0.0
    } else {
        align_threaded(s1, s2, normalizer)
    }
}
