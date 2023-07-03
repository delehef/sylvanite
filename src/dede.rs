#![allow(dead_code)]
use std::fmt::{Debug, Display, Formatter};
use std::ops::{Index, IndexMut};

use either::Either;

pub trait Matrix<T>: Index<(usize, usize), Output = T> {
    fn is_square(&self) -> bool;
    fn nrows(&self) -> usize;
    fn ncols(&self) -> usize;
    // fn row(&self, i: usize) -> Vec<T>;
    // fn col(&self, j: usize) -> &dyn Iterator<Item = &T>;
}

pub struct VecMatrix<T> {
    m: Vec<T>,
    r: usize,
    c: usize,
}
impl<T: Debug> Debug for VecMatrix<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(f)?;
        for i in 0..self.r {
            for j in 0..self.c {
                write!(f, "{:10.3?}", self[(i, j)])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl<T> VecMatrix<T> {
    pub fn new(r: usize, c: usize) -> VecMatrix<T> {
        let mut m: Vec<T> = Vec::<T>::with_capacity(r * c);
        unsafe {
            m.set_len(r * c);
        }
        VecMatrix { m, r, c }
    }

    pub fn masked<'a>(&'a self, is: &'a [usize], js: &'a [usize]) -> MaskedMatrix<'a, T> {
        assert!(is.iter().all(|&i| i < self.r));
        assert!(js.iter().all(|&j| j < self.c));
        MaskedMatrix { m: self, is: Either::Left(is), js: Either::Left(js) }
    }

    pub fn masked_from_iter<'a>(
        &'a self,
        is: impl Iterator<Item = usize>,
        js: impl Iterator<Item = usize>,
    ) -> MaskedMatrix<'a, T> {
        MaskedMatrix { m: self, is: Either::Right(is.collect()), js: Either::Right(js.collect()) }
    }

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.m.iter()
    }
}

impl<T> Matrix<T> for VecMatrix<T> {
    fn is_square(&self) -> bool {
        self.r == self.c
    }

    fn nrows(&self) -> usize {
        self.r
    }

    fn ncols(&self) -> usize {
        self.c
    }
    // fn row(&'a self, i: usize) -> &Self::Iter {
    //     &(0..self.ncols()).map(|j| &self[(i, j)])
    // }
    // fn col(&self, j: usize) -> &mut dyn Iterator<Item = &T> {
    //     &mut (0..self.nrows()).map(|i| &self[(i, j)])
    // }
}

impl<T: Clone> VecMatrix<T> {
    pub fn with_elem(r: usize, c: usize, x: T) -> VecMatrix<T> {
        VecMatrix { m: vec![x; r * c], r, c }
    }
}

impl<T: Display> Display for VecMatrix<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.r {
            for j in 0..self.c {
                write!(f, "{:2.3} ", self.m[i * self.c + j])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl<T> Index<usize> for VecMatrix<T> {
    type Output = T;

    fn index(&self, i: usize) -> &Self::Output {
        debug_assert!(i < self.m.len());
        &self.m[i]
    }
}

impl<T> Index<(usize, usize)> for VecMatrix<T> {
    type Output = T;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        debug_assert!(i < self.r);
        debug_assert!(j < self.c);
        &self.m[i * self.r + j]
    }
}

impl<T> IndexMut<(usize, usize)> for VecMatrix<T> {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        debug_assert!(i < self.r);
        debug_assert!(j < self.c);
        &mut self.m[i * self.c + j]
    }
}

impl<T: Clone> std::clone::Clone for VecMatrix<T> {
    fn clone(&self) -> Self {
        VecMatrix { c: self.c, r: self.r, m: self.m.clone() }
    }
}

pub struct MaskedMatrix<'a, T> {
    m: &'a dyn Matrix<T>,
    is: Either<&'a [usize], Vec<usize>>,
    js: Either<&'a [usize], Vec<usize>>,
}

impl<'a, T> MaskedMatrix<'a, T> {
    pub fn is(&self) -> &[usize] {
        match &self.is {
            Either::Left(is) => is,
            Either::Right(is) => is.as_slice(),
        }
    }

    pub fn js(&self) -> &[usize] {
        match &self.js {
            Either::Left(js) => js,
            Either::Right(js) => js.as_slice(),
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = &T> {
        MaskedMatrixIterator { m: self, p: 0 }
    }
}

impl<'a> MaskedMatrix<'a, f32> {
    pub fn max(&self) -> f32 {
        self.iter().fold(0f32, |a, &b| a.max(b))
    }

    pub fn min(&self) -> f32 {
        self.iter().fold(f32::INFINITY, |a, &b| a.min(b))
    }
}

impl<T> Matrix<T> for MaskedMatrix<'_, T> {
    fn is_square(&self) -> bool {
        self.ncols() == self.nrows()
    }
    fn nrows(&self) -> usize {
        self.is.len()
    }
    fn ncols(&self) -> usize {
        self.js.len()
    }
    // fn row(&self, i: usize) -> &mut dyn Iterator<Item = &T> {
    //     &mut (0..self.ncols()).map(|j| &self[(i, j)])
    // }
    // fn col(&self, j: usize) -> &mut dyn Iterator<Item = &T> {
    //     &mut (0..self.nrows()).map(|i| &self[(i, j)])
    // }
}

impl<T> Index<(usize, usize)> for MaskedMatrix<'_, T> {
    type Output = T;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        debug_assert!(i < self.m.nrows());
        debug_assert!(j < self.m.ncols());
        &self.m[(self.is[i], self.js[j])]
    }
}

impl<T: Display> Display for MaskedMatrix<'_, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.nrows() {
            for j in 0..self.ncols() {
                write!(f, "{:2.3} ", self.index((i, j)))?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

pub struct MaskedMatrixIterator<'a, T> {
    m: &'a MaskedMatrix<'a, T>,
    p: usize,
}

impl<'a, T> Iterator for MaskedMatrixIterator<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.p >= self.m.is().len() * self.m.js().len() {
            None
        } else {
            let i = self.p / self.m.ncols();
            let j = self.p % self.m.ncols();
            self.p += 1;
            Some(&self.m[(i, j)])
        }
    }
}

pub fn view<'a, T: Index<usize>>(
    v: &'a T,
    is: impl IntoIterator<Item = &'a usize> + 'a,
) -> impl Iterator<Item = &'a <T as Index<usize>>::Output> + 'a
where
    <T as Index<usize>>::Output: Sized,
{
    is.into_iter().map(|i| &v[*i])
}

pub fn clonable_view_cloned<'a, T, I>(
    v: &'a T,
    is: I,
) -> impl Clone + Iterator<Item = <T as Index<usize>>::Output> + 'a
where
    T: Index<usize> + 'a,
    <T as Index<usize>>::Output: Sized + Clone,

    I: IntoIterator<Item = &'a usize> + 'a,
    <I as IntoIterator>::IntoIter: Clone,
{
    is.into_iter().map(move |i| v[*i].clone())
}

pub fn view_cloned<'a, T, I>(
    v: &'a T,
    is: I,
) -> impl Iterator<Item = <T as Index<usize>>::Output> + 'a
where
    T: Index<usize> + 'a,
    <T as Index<usize>>::Output: Sized + Clone,

    I: IntoIterator<Item = &'a usize> + 'a,
{
    is.into_iter().map(move |i| v[*i].clone())
}
