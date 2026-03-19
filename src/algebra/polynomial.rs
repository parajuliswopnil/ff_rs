//! Polynomial implemention

use std::cmp::max;

use crate::algebra::field::Field;

/// Polynomial representing struct
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Polynomial<F>
where
    F: Field,
{
    cofficients: Vec<F>,
}

impl<F> Polynomial<F>
where
    F: Field,
{
    /// creates new instance of polynomials with the list of cofficients
    /// removes the trailing zeros from the cofficient list
    /// the index of the cofficients represents the degree of x
    /// length - 1 represents the degree of the polynomial
    pub fn new(cofficients: Vec<F>) -> Self {
        let mut cofficients = cofficients;
        Self::normalize(&mut cofficients);
        Self { cofficients }
    }

    /// returns degree of the polynomial
    pub fn degree(&self) -> Option<usize> {
        if self.is_zero() {
            return None;
        }
        Some(self.cofficients.len() - 1)
    }

    /// zero polynomial
    pub fn zero() -> Self {
        Self {
            cofficients: vec![],
        }
    }

    /// one
    pub fn one() -> Self {
        Self {
            cofficients: vec![F::one()],
        }
    }

    /// add two polynomials
    pub fn add(&self, other: &Self) -> Self {
        let mut result = vec![];
        let zero = F::zero();
        for i in 0..max(self.cofficients.len(), other.cofficients.len()) {
            let a = self.cofficients.get(i).unwrap_or(&zero);
            let b = other.cofficients.get(i).unwrap_or(&zero);
            let res = a.add(b);
            result.push(res);
        }

        Self::new(result)
    }

    /// subtract two polynomials
    pub fn sub(&self, other: &Self) -> Self {
        let mut result = vec![];
        let zero = F::zero();
        for i in 0..max(self.cofficients.len(), other.cofficients.len()) {
            let a = self.cofficients.get(i).unwrap_or(&zero);
            let b = other.cofficients.get(i).unwrap_or(&zero);
            let res = a.sub(b);
            result.push(res);
        }

        Self::new(result)
    }

    /// multiply two polynomials
    pub fn mul(&self, other: &Self) -> Self {
        let mut intermediate_result = vec![];
        if *self == Self::zero() || *other == Self::zero() {
            return Self::zero();
        }

        for i in 0..self.cofficients.len() {
            let mut in_res = vec![];
            for _ in 0..i {
                in_res.push(F::zero());
            }
            let mut mul_res = other.scalar_mul(&self.cofficients[i]);
            in_res.append(&mut mul_res.cofficients);
            intermediate_result.push(Polynomial {
                cofficients: in_res,
            });
        }
        while intermediate_result.len() > 1 {
            let a = intermediate_result.pop().unwrap();
            let b = intermediate_result.pop().unwrap();
            let aplusb = a.add(&b);
            intermediate_result.push(aplusb);
        }
        intermediate_result[0].add(&Self::zero())
    }

    /// scalar multiply
    pub fn scalar_mul(&self, scalar: &F) -> Self {
        let mut result = vec![];
        for i in 0..self.cofficients.len() {
            let res = self.cofficients[i].mul(scalar);
            result.push(res);
        }
        Self::normalize(&mut result);
        Self {
            cofficients: result,
        }
    }

    /// polynomial evaluation at x
    pub fn eval(&self, x: &F) -> F {
        let mut result = F::zero();
        () = self
            .cofficients
            .iter()
            .rev()
            .map(|coeff| {
                result = result.mul(x);
                result = result.add(coeff);
            })
            .collect();
        result
    }

    /// lead cofficient
    pub fn lead_coeff(&self) -> Option<F> {
        self.cofficients.last().cloned()
    }

    /// is zero
    pub fn is_zero(&self) -> bool {
        Self::zero() == *self
    }

    /// polynomial division
    pub fn div(&self, other: &Self) -> Option<(Self, Self)> {
        if other.is_zero() {
            return None;
        }
        let mut remainder = self.clone();
        let mut quotient = Self::zero();

        while remainder.degree().is_some() && remainder.degree().unwrap() >= other.degree().unwrap()
        {
            let k = remainder.degree().unwrap() - other.degree().unwrap();

            let c = remainder
                .lead_coeff()
                .unwrap()
                .mul(&other.lead_coeff().unwrap().inv().unwrap());

            let mut tx_coeff = vec![F::zero(); k + 1];
            tx_coeff[k] = F::one();

            let tx = Self::new(tx_coeff).scalar_mul(&c);

            quotient = quotient.add(&tx);
            let rem = tx.mul(other);
            remainder = remainder.sub(&rem);
        }
        Some((quotient, remainder))
    }

    fn normalize(cofficients: &mut Vec<F>) {
        while let Some(last) = cofficients.last() {
            if *last == F::zero() {
                cofficients.pop();
                continue;
            }
            break;
        }
    }
}

mod tests {
    #![allow(unused_imports)]
    use super::*;
    use crate::algebra::fp::Fp;

    #[test]
    fn test_polynomial_division_exact() {
        // (x^2 + 3x + 2) / (x + 1) = x + 2
        // quotient = x + 2
        // remainder = 0

        let f = Polynomial::<Fp<7>>::new(vec![Fp::new(2), Fp::new(3), Fp::new(1)]);

        let g = Polynomial::<Fp<7>>::new(vec![Fp::new(1), Fp::new(1)]);

        let (q, r) = f.div(&g).unwrap();

        let expected_q = Polynomial::<Fp<7>>::new(vec![Fp::new(2), Fp::new(1)]);

        let expected_r = Polynomial::<Fp<7>>::zero();

        assert_eq!(q, expected_q);
        assert_eq!(r, expected_r);

        let f_reconstructed = q.mul(&g).add(&r);
        assert_eq!(f, f_reconstructed);
    }

    #[test]
    fn test_polynomial_division_with_remainder() {
        // (x^2 + 1) / (x + 1)
        // quotient = x - 1
        // remainder = 2
        // in F7: -1 ≡ 6

        let f = Polynomial::<Fp<7>>::new(vec![Fp::new(1), Fp::new(0), Fp::new(1)]);

        let g = Polynomial::<Fp<7>>::new(vec![Fp::new(1), Fp::new(1)]);

        let (q, r) = f.div(&g).unwrap();

        let expected_q = Polynomial::<Fp<7>>::new(vec![Fp::new(6), Fp::new(1)]);

        let expected_r = Polynomial::<Fp<7>>::new(vec![Fp::new(2)]);

        assert_eq!(q, expected_q);
        assert_eq!(r, expected_r);

        let f_reconstructed = q.mul(&g).add(&r);
        assert_eq!(f, f_reconstructed);
    }

    #[test]
    fn test_divisor_higher_degree() {
        // (x + 1) / (x^2 + 1)
        // quotient = 0
        // remainder = x + 1

        let f = Polynomial::<Fp<7>>::new(vec![Fp::new(1), Fp::new(1)]);

        let g = Polynomial::<Fp<7>>::new(vec![Fp::new(1), Fp::new(0), Fp::new(1)]);

        let (q, r) = f.div(&g).unwrap();

        assert_eq!(q, Polynomial::<Fp<7>>::zero());
        assert_eq!(r, f);

        let f_reconstructed = q.mul(&g).add(&r);
        assert_eq!(f, f_reconstructed);
    }

    #[test]
    fn test_division_by_zero_polynomial() {
        let f = Polynomial::<Fp<7>>::new(vec![Fp::new(1), Fp::new(2)]);

        let g = Polynomial::<Fp<7>>::zero();

        assert!(f.div(&g).is_none());
    }
}
