use std::cmp::*;

pub fn align_sw<T1, T2, V>(s1: T1, s2: T2, normalizer: &dyn Fn(usize, usize) -> f32) -> f32
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

    max_score as f32 / normalizer(s1.len() - 1, s2.len() - 1)
}

type Path = Vec<(usize, usize)>;

fn align_threaded<T1, T2, V>(s1: T1, s2: T2, normalizer: &dyn Fn(f32, f32) -> f32) -> f32
where
    T1: AsRef<[V]>,
    T2: AsRef<[V]>,
    V: PartialEq,
{
    todo!()
}
