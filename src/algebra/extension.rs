//! Finite field extensions

use std::marker::PhantomData;

use crate::algebra::{field::Field, polynomial::Polynomial};

/// trait for the Modulus type
pub trait Modulus<F: Field>
where
    Self: Sized,
{
    /// returns modulus polynomial
    /// instead of doing this, later change this to returns &'static Polynomial...
    /// so that the modulus is not created new every time we call this.
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

    /// return polynomial associated to self
    pub fn polynomial(&self) -> Polynomial<F> {
        self.polynomial.clone()
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

    /// one
    pub fn one() -> Self {
        Self::new(Polynomial::<F>::one())
    }

    /// zero
    pub fn zero() -> Self {
        Self::new(Polynomial::<F>::zero())
    }

    /// extended euclidean polynomial to calculate inverse of a polynomial
    fn extended_euclidean_polynomial(&self) -> Option<Self> {
        let mut r1 = M::inner();
        let mut r2 = self.polynomial.clone();

        let mut t1 = Polynomial::<F>::zero();
        let mut t2 = Polynomial::<F>::one(); // polynomial 1

        while r2 != Polynomial::<F>::zero() {
            let qnr = r1.div(&r2).unwrap(); // can safely unwrap because r2 is not zero

            let t = t1.sub(&t2.mul(&qnr.0));

            r1 = r2;
            r2 = qnr.1;

            t1 = t2;
            t2 = t;
        }

        if r1 != Self::one().polynomial {
            return None;
        }

        let inv_lead = r1.lead_coeff()?.inv()?;
        let t1 = t1.scalar_mul(&inv_lead);

        // reduce mod modulus
        let (_, t1) = t1.div(&M::inner()).unwrap();
        Some(Self::new(t1.scalar_mul(&inv_lead)))
    }
}

mod tests {
    #![allow(unused_imports, dead_code)]
    use std::vec;

    use crate::algebra::fp::Fp;

    use super::*;

    type SomeModulus = Polynomial<Fp<2>>;

    impl Modulus<Fp<2>> for SomeModulus {
        fn inner() -> Polynomial<Fp<2>> {
            Polynomial::new(vec![
                Fp::<2>::new(1),
                Fp::<2>::new(1),
                Fp::<2>::new(0),
                Fp::<2>::new(1),
            ])
        }
    }

    type GF8 = ExtensionField<Fp<2>, SomeModulus>;

    #[test]
    fn test_inverse_basic() {
        let a = GF8::new(Polynomial::new(vec![
            Fp::<2>::new(1),
            Fp::<2>::new(0),
            Fp::<2>::new(1),
        ]));

        let inv = a.inv().unwrap();
        let one = a.mul(&inv).unwrap();

        assert_eq!(one, GF8::one());
    }

    #[test]
    fn test_inverse_of_one() {
        let one = GF8::one();
        let inv = one.inv();

        assert_eq!(one, inv.unwrap());
    }

    #[test]
    fn test_inverse_of_zero() {
        let zero = GF8::zero();
        let inv = zero.inv();

        assert!(inv.is_none());
    }

    #[test]
    fn test_multiplicative_identity() {
        let a = GF8::new(Polynomial::new(vec![Fp::<2>::new(1), Fp::<2>::new(1)]));
        let result = a.mul(&GF8::one());

        assert_eq!(a, result.unwrap());
    }

    #[test]
    fn test_additive_identity() {
        let a = GF8::new(Polynomial::new(vec![Fp::<2>::new(1), Fp::<2>::new(1)]));
        let result = a.add(&GF8::zero());

        assert_eq!(a, result);
    }

    #[test]
    fn test_distributivity() {
        let a = GF8::new(Polynomial::new(vec![
            Fp::<2>::new(1),
            Fp::<2>::new(0),
            Fp::<2>::new(1),
        ]));
        let b = GF8::new(Polynomial::new(vec![Fp::<2>::new(1), Fp::<2>::new(1)]));
        let c = GF8::new(Polynomial::new(vec![Fp::<2>::new(0), Fp::<2>::new(1)]));

        let left = a.mul(&b.add(&c)).unwrap();
        let right = a.mul(&b).unwrap().add(&a.mul(&c).unwrap());

        assert_eq!(left, right);
    }

    #[test]
    fn test_degree_reduction() {
        let a = GF8::new(Polynomial::new(vec![
            Fp::<2>::new(0),
            Fp::<2>::new(0),
            Fp::<2>::new(1),
        ])); // x^2
        let b = GF8::new(Polynomial::new(vec![
            Fp::<2>::new(0),
            Fp::<2>::new(0),
            Fp::<2>::new(1),
        ])); // x^2

        let result = a.mul(&b).unwrap();

        // degree must be < modulus degree (3 in GF(2^3))
        assert!(result.polynomial.degree().unwrap() < 3);
    }

    #[test]
    fn test_inverse_symmetry() {
        let a = GF8::new(Polynomial::new(vec![
            Fp::<2>::new(1),
            Fp::<2>::new(0),
            Fp::<2>::new(1),
        ]));

        let inv1 = a.inv().unwrap();
        let inv2 = inv1.inv().unwrap();

        assert_eq!(inv2, a);
    }

    #[test]
    fn test_all_nonzero_have_inverse() {
        let elements = vec![
            Polynomial::new(vec![Fp::<2>::new(1)]),
            Polynomial::new(vec![Fp::<2>::new(0), Fp::<2>::new(1)]),
            Polynomial::new(vec![Fp::<2>::new(1), Fp::<2>::new(1)]),
            Polynomial::new(vec![Fp::<2>::new(0), Fp::<2>::new(0), Fp::<2>::new(1)]),
            Polynomial::new(vec![Fp::<2>::new(1), Fp::<2>::new(0), Fp::<2>::new(1)]),
            Polynomial::new(vec![Fp::<2>::new(0), Fp::<2>::new(1), Fp::<2>::new(1)]),
            Polynomial::new(vec![Fp::<2>::new(1), Fp::<2>::new(1), Fp::<2>::new(1)]),
        ];

        for coeffs in elements {
            let a = GF8::new(coeffs);
            let inv = a.inv().expect("must have inverse");

            let prod = a.mul(&inv).unwrap();
            assert_eq!(prod, GF8::one());
        }
    }
}
