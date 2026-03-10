//! Defines the trait to be implemented by a field

/// Field trait defination
pub trait Field {
    /// Creates new instance of the trait
    /// ## Arguments
    /// - `a`: Field element
    fn new(a: u64) -> Self;
    // additive operations

    /// adds self and other field elements and returns the result
    fn add(&self, other: &Self) -> Self;

    /// subtracts self and other field elements and returns the result
    fn sub(&self, other: &Self) -> Self;

    /// additive identity: zero of the field
    fn zero() -> Self;

    /// additive inverse: -self of the field such that self + additive_inverse = additive identity
    fn neg(&self) -> Self;

    /// multipliply self with other field elements and returns the result
    fn mul(&self, other: &Self) -> Self;

    /// multiplicative identity: one of the group such that self * multiplicative identity = self
    fn one() -> Self;

    /// multiplicative inverse: self ^ -1 of the group such that self * mulitplicative inverse =
    /// multiplicative identity
    fn inv(&self) -> Option<Self>
    where
        Self: Sized;

    // unity operations

    /// raises self to the power of x and returns the result
    fn pow(&self, x: u64) -> Self;
}
