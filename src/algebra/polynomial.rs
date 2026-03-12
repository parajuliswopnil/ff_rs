use std::cmp::max;

use crate::algebra::field::Field;

/// Polynomial representing struct
#[derive(Debug, PartialEq, Eq)]
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
    pub fn degree(&self) -> u64 {
        max((self.cofficients.len() - 1) as u64, 0)
    }

    /// zero polynomial
    pub fn zero() -> Self {
        Self {
            cofficients: vec![],
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

    /// polynomial division
    pub fn div(&self, other: &Self) -> (Self, Self) {
        if self.degree() < other.degree() {
            return (Self::zero(), Self::zero().add(self));
        }
        todo!()
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
