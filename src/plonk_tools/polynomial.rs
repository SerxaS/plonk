use core::fmt;
use halo2::{arithmetic::Field, halo2curves::ff::PrimeField};
use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

#[derive(Debug, Clone)]
pub struct Polynomial<F> {
    pub(crate) coeff: Vec<F>,
}

impl<F> Polynomial<F>
where
    F: Field + PrimeField,
{
    /// Creates a new struct field.
    pub fn new(coeff: Vec<F>) -> Self {
        Self { coeff }
    }

    /// Returns an iterator over the slice.
    pub fn iter(&self) -> impl Iterator<Item = &F> {
        self.coeff.iter()
    }

    /// Returns the number of elements in the polynomial.
    pub fn len(&self) -> usize {
        self.coeff.len()
    }

    /// Removes the last element from a vector and returns it, or None if it is empty.
    pub fn pop(&mut self) -> Option<F> {
        self.coeff.pop()
    }

    /// Appends an element to the back of a collection.
    pub fn push(&mut self, value: F) {
        self.coeff.push(value)
    }

    /// Computes polynomial addition.
    pub fn add_poly(&self, rhs: &Self) -> Self {
        let mut len = 0;

        if self.len() >= rhs.len() {
            len = self.len();
        } else {
            len = rhs.len();
        }
        let mut sum_poly = Polynomial::new(vec![F::ZERO; len]);

        for (i, _) in self.iter().enumerate() {
            sum_poly[i] += self[i];
        }

        for (i, _) in rhs.iter().enumerate() {
            sum_poly[i] += rhs[i];
        }
        sum_poly
    }

    /// Add Fr element to polynomial.
    pub fn add_fr_poly(&self, rhs: &F) -> Self {
        let len = self.len();
        let mut sum_poly = Polynomial::new(vec![F::ZERO; len]);

        for (i, _) in self.iter().enumerate() {
            sum_poly[i] += self[i];
        }
        sum_poly[0] += rhs;
        sum_poly
    }

    /// Computes polynomial subtraction.
    pub fn sub_poly(&self, rhs: &Self) -> Self {
        let mut len = 0;

        if self.len() >= rhs.len() {
            len = self.len();
        } else {
            len = rhs.len();
        }
        let mut sub_poly = Polynomial::new(vec![F::ZERO; len]);

        for (i, _) in self.iter().enumerate() {
            sub_poly[i] += self[i];
        }

        for (i, _) in rhs.iter().enumerate() {
            sub_poly[i] -= rhs[i];
        }
        sub_poly
    }

    /// Substract Fr element to polynomial.
    pub fn sub_fr_poly(&self, rhs: &F) -> Self {
        let len = self.len();
        let mut sub_poly = Polynomial::new(vec![F::ZERO; len]);

        for (i, _) in self.iter().enumerate() {
            sub_poly[i] += self[i];
        }
        sub_poly[0] -= rhs;
        sub_poly
    }

    /// Computes polynomial multiplication.
    pub fn mul_poly(&self, rhs: &Self) -> Self {
        let poly_len = self.len() + rhs.len() - 1;
        let mut mul_poly = Polynomial::new(vec![F::ZERO; poly_len]);

        for i in 0..self.len() {
            for j in 0..rhs.len() {
                mul_poly[i + j] += self[i] * rhs[j];
            }
        }
        mul_poly
    }

    /// Multiply Fr element to polynomial.
    pub fn mul_fr_poly(&self, rhs: &F) -> Self {
        let poly_len = self.len();
        let mut mul_poly = Polynomial::new(vec![F::ZERO; poly_len]);

        for i in 0..self.len() {
            for j in 0..1 {
                mul_poly[i + j] += self[i] * rhs;
            }
        }
        mul_poly
    }

    /// Computes polynomial long division.
    pub fn div_poly(&mut self, den: &Self) -> (Self, Self) {
        if den.len() > self.len() {
            return (Self::new(vec![F::ZERO]), self.clone());
        }
        let diff = self.len() - den.len();
        let mut quotient = Polynomial::new(vec![F::ZERO; diff + 1]);

        for i in (0..quotient.len()).rev() {
            let n_idx = self.len() - 1 - diff + i;
            let inv_d = den[den.len() - 1].invert().unwrap();
            quotient[i] = self[n_idx].mul(&inv_d);

            for j in 0..den.len() {
                self[n_idx - j] -= quotient[i].mul(&den[den.len() - j - 1]);
            }
        }

        for i in (1..self.len()).rev() {
            if self[i] == F::ZERO {
                self.pop();
            } else {
                break;
            }
        }
        let remainder = self.clone();
        (quotient, remainder)
    }

    /// Computes lagrange interpolation polynomial.
    pub fn lagrange(x_val: &Vec<F>, y_val: &Vec<F>) -> Self {
        let mut int_poly = Polynomial::new(vec![F::ZERO]);

        for i in 0..x_val.len() {
            let mut term = Polynomial::new(vec![y_val[i]]);

            for j in 0..x_val.len() {
                if j != i {
                    let mut num = Polynomial::new(vec![x_val[j].neg(), F::ONE]);
                    let den = Polynomial::new(vec![x_val[i] - x_val[j]]);
                    term = term.mul_poly(&num.div_poly(&den).0);
                }
            }
            int_poly = term.add_poly(&int_poly);
        }
        int_poly
    }

    /// Computes zero(vanishing) polynomial Zh(X).
    pub fn zero_poly(n: usize) -> Self {
        let mut van_poly = Polynomial::new(vec![F::ONE.neg()]);

        for _ in 0..(n - 1) {
            van_poly.push(F::ZERO);
        }
        van_poly.push(F::ONE);
        van_poly
    }

    /// Computes permutation polynomial's accumulator.
    pub fn perm_poly_acc(
        a: &Vec<F>,
        b: &Vec<F>,
        c: &Vec<F>,
        sigma_star_1: &Vec<F>,
        sigma_star_2: &Vec<F>,
        sigma_star_3: &Vec<F>,
        rou_dom: &Vec<F>,
        k1: F,
        k2: F,
        beta: F,
        gamma: F,
    ) -> Self {
        let mut acc = vec![F::ONE; a.len()];

        for i in 0..(a.len() - 1) {
            acc[i + 1] = ((a[i] + (beta * rou_dom[i]) + gamma)
                * (b[i] + (beta * k1 * rou_dom[i]) + gamma)
                * (c[i] + (beta * k2 * rou_dom[i]) + gamma))
                * (((a[i] + (beta * sigma_star_1[i]) + gamma)
                    * (b[i] + (beta * sigma_star_2[i]) + gamma)
                    * (c[i] + (beta * sigma_star_3[i]) + gamma))
                    .invert()
                    .unwrap());
        }
        let mut acc_dom = vec![F::ONE; a.len()];

        for i in 0..1 {
            let mut tmp_dom = vec![F::ONE];

            for j in 0..acc_dom.len() {
                tmp_dom[i] = acc[j] * tmp_dom[i];
                acc_dom[j] = tmp_dom[i];
            }
        }
        let acc_poly = Polynomial::lagrange(&rou_dom, &acc_dom);
        acc_poly
    }

    /// Computes shifted permutation polynomial z(Xω).
    pub fn shift_perm_poly_z(&self, rou_dom: &Vec<F>) -> Self {
        // Extends rou domain for multiply each z_x elements.
        let mut ext_rou = Polynomial::new(vec![F::ONE; self.len()]);

        for i in 0..rou_dom.len() {
            ext_rou[i] = rou_dom[i];
        }

        for i in 0..(rou_dom.len() - 1) {
            ext_rou[i + rou_dom.len()] = rou_dom[i];
        }

        let mut shifted_perm_poly_z = Polynomial::new(vec![F::ONE; self.len()]);

        for i in 0..ext_rou.len() {
            shifted_perm_poly_z[i] = ext_rou[i] * self[i];
        }
        shifted_perm_poly_z
    }

    /// Computes L1(X) polynomial. L1 refers to a Lagrange basis polynomial over our
    /// roots of unity H. Specifically L1(1)=1, but takes the value 0  on each of
    /// the other roots of unity.     
    pub fn l_1_poly(rou_dom: &Vec<F>) -> Self {
        let mut l_1_dom = vec![F::ONE; rou_dom.len()];

        for i in 1..rou_dom.len() {
            l_1_dom[i] = F::ZERO;
        }
        let l1_poly = Polynomial::lagrange(&rou_dom, &l_1_dom);
        l1_poly
    }

    /// Splits t(X) polynomial.
    pub fn split_t_x(t_x: Polynomial<F>) -> (Self, Self, Self) {
        let mut t_x_lo = Polynomial::new(vec![F::ONE; t_x.len() / 3]);
        let mut t_x_mid = Polynomial::new(vec![F::ONE; t_x.len() / 3]);
        let mut t_x_hi = Polynomial::new(vec![F::ONE; t_x.len() / 3]);

        for i in 0..t_x_lo.len() {
            t_x_lo[i] = t_x[i];
        }

        for i in 0..t_x_mid.len() {
            t_x_mid[i] = t_x[i + t_x_mid.len()];
        }

        for i in 0..t_x_hi.len() {
            t_x_hi[i] = t_x[i + (2 * (t_x_hi.len()))];
        }
        (t_x_lo, t_x_mid, t_x_hi)
    }

    /// Computes blind_10 and blind_11 polynomials from random Fr elements.
    pub fn blind_poly(t_x_n: usize, blind_sc_10: F, blind_sc_11: F) -> (Self, Self) {
        let mut b_10_poly = Polynomial::new(vec![F::ZERO; t_x_n]);
        b_10_poly.push(blind_sc_10);
        let mut b_11_poly = Polynomial::new(vec![F::ZERO; t_x_n]);
        b_11_poly.push(blind_sc_11);

        (b_10_poly, b_11_poly)
    }
}

/// Shows polynomial in the equation form.
impl<F: Field> Display for Polynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let result: Vec<String> = self
            .coeff
            .iter()
            .enumerate()
            .map(|(i, c)| format!("({:?})x^{}", c, i))
            .collect();
        write!(f, "{}", result.join(" + "))
    }
}

impl<F: Field> Index<usize> for Polynomial<F> {
    type Output = F;

    fn index(&self, idx: usize) -> &Self::Output {
        self.coeff.index(idx)
    }
}

impl<F: Field> IndexMut<usize> for Polynomial<F> {
    fn index_mut(&mut self, idx: usize) -> &mut F {
        self.coeff.index_mut(idx)
    }
}
