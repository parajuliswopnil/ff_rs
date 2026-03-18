//! Finite field extensions

use std::marker::PhantomData;

use crate::algebra::{field::Field, polynomial::Polynomial};

/// trait for the Modulus type
pub trait Modulus<F: Field>
where
    Self: Sized,
{
    /// returns modulus polynomial
    fn inner() -> Polynomial<F>;
}

/// Struct representing the extension field element
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtensionField<F, M>
where
    F: Field,
    M: Modulus<F>,
{
    polynomial: Polynomial<F>,
    phantom_data: PhantomData<M>,
}

impl<F, M> ExtensionField<F, M>
where
    F: Field,
    M: Modulus<F>,
{
    /// creates new instance of extension field
    pub fn new(polynomial: Polynomial<F>) -> Self {
        Self {
            polynomial,
            phantom_data: PhantomData,
        }
    }

    /// add
    pub fn add(&self, other: &Self) -> Self {
        Self::new(self.polynomial.add(&other.polynomial))
    }

    /// multiply
    pub fn mul(&self, other: &Self) -> Option<Self> {
        if let Some(polynomial) = self.polynomial.mul(&other.polynomial).div(&M::inner()) {
            return Some(Self::new(polynomial.1));
        }
        None
    }

    /// return inverse  
    pub fn inv(&self) -> Option<Self> {
        if self.polynomial.is_zero() {
            return None;
        }
        self.extended_euclidean_polynomial()
    }

    /// extended euclidean polynomial to calculate inverse of a polynomial
    fn extended_euclidean_polynomial(&self) -> Option<Self> {
        let mut r1 = M::inner();
        let mut r2 = self.polynomial.clone();

        let mut t1 = Polynomial::<F>::zero();
        let mut t2 = Polynomial::<F>::new(vec![F::new(1)]); // polynomial 1

        while r2 != Polynomial::<F>::zero() {
            let qnr = r1.div(&r2).unwrap(); // can safely unwrap because r2 is not zero

            let t = t1.sub(&t2.mul(&qnr.0));

            r1 = r2;
            r2 = qnr.1;

            t1 = t2;
            t2 = t;
        }

        if r1 != Polynomial::<F>::new(vec![F::new(1)]) {
            return None;
        }

        let inv_lead = r1.lead_coeff()?.inv()?;
        let t1 = t1.scalar_mul(&inv_lead);

        // reduce mod modulus
        let (_, t1) = t1.div(&M::inner()).unwrap();
        Some(Self::new(t1.scalar_mul(&inv_lead)))
    }
}
