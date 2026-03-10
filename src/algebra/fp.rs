//! Implementation of a prime field struct

use crate::algebra::field::Field;

/// Field struct defination
#[derive(Debug)]
pub struct Fp<const PR: u64> {
    value: u64,
}

impl<const PR: u64> Field for Fp<PR> {
    fn new(a: u64) -> Self {
        let value = a % PR;
        Self { value }
    }
    // additive operations
    fn add(&self, other: &Self) -> Self {
        let result = (self.value + other.value) % PR;

        Self { value: result }
    }

    fn sub(&self, other: &Self) -> Self {
        let result = (self.value + PR - other.value) % PR;

        Self { value: result }
    }

    // additive identity
    fn zero() -> Self {
        Self { value: 0 }
    }

    // additive inverse
    fn neg(&self) -> Self {
        let value = (PR - self.value) % PR;
        Self { value }
    }

    // multiplicative operations
    fn mul(&self, other: &Self) -> Self {
        let result = (self.value * other.value) % PR;

        Self { value: result }
    }

    // multiplicative identity
    fn one() -> Self {
        Self { value: 1 }
    }

    // multiplicative inverse
    fn inv(&self) -> Option<Self>
    where
        Self: Sized,
    {
        if self.value == 0 {
            return None;
        }
        Self::extended_euclidean(self.value)
    }

    // unity operations
    fn pow(&self, x: u64) -> Self {
        let mut result = 1;
        let mut a = self.value;
        let mut b = x;

        while b > 0 {
            if b % 2 == 1 {
                result = (result * a) % PR;
            }

            a = (a * a) % PR;
            b /= 2;
        }

        Self { value: result }
    }
}

impl<const PR: u64> Fp<PR> {
    fn extended_euclidean(value: u64) -> Option<Self> {
        let mut r1 = PR as i64;
        let mut r2 = value as i64;
        let mut t1 = 0;
        let mut t2 = 1;

        while r2 != 0 {
            let q = r1 / r2;
            let r = r1 % r2;
            let t = t1 - t2 * q;

            r1 = r2;
            r2 = r;

            t1 = t2;
            t2 = t;
        }

        if r1 != 1 {
            return None;
        }

        t1 %= PR as i64;

        if t1 < 0 {
            return Some(Self {
                value: (PR as i64 + t1) as u64,
            });
        }

        Some(Self { value: t1 as u64 })
    }

    #[allow(dead_code)]
    fn fermats_little_theorem(&self) -> Option<Self> {
        // for prime field only
        let value = Fp::<PR>::new(self.value);

        Some(value.pow(PR - 2))
    }
}

#[test]
fn test_extended() {
    let field_elem = Fp::<26>::new(11);
    let result = field_elem.inv();

    assert_eq!(19, result.unwrap().value);

    let field_elem = Fp::<7>::new(3);
    let result = field_elem.inv();

    assert_eq!(5, result.unwrap().value);

    let field_elem = Fp::<7>::new(1);
    let result = field_elem.inv();

    assert_eq!(1, result.unwrap().value);

    let field_elem = Fp::<7>::new(6);
    let result = field_elem.inv();

    assert_eq!(6, result.unwrap().value);
}

#[test]
fn test_fermats() {
    let field_elem = Fp::<7>::new(3);
    let result = field_elem.fermats_little_theorem();

    assert_eq!(5, result.unwrap().value);

    let field_elem = Fp::<7>::new(1);
    let result = field_elem.fermats_little_theorem();

    assert_eq!(1, result.unwrap().value);

    let field_elem = Fp::<7>::new(6);
    let result = field_elem.fermats_little_theorem();

    assert_eq!(6, result.unwrap().value);
}

#[test]
fn test_pow() {
    let field_elem = Fp::<7>::new(3);
    let result = field_elem.pow(2);

    assert_eq!(2, result.value);

    let field_elem = Fp::<7>::new(4);
    let result = field_elem.pow(2);

    assert_eq!(2, result.value);

    let field_elem = Fp::<7>::new(5);
    let result = field_elem.pow(2);

    assert_eq!(4, result.value);

    let field_elem = Fp::<7>::new(6);
    let result = field_elem.pow(2);

    assert_eq!(1, result.value);
}
