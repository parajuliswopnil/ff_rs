//! AES implementation
use crate::algebra::{
    extension::{ExtensionField, Modulus},
    field::Field,
    fp::Fp,
    polynomial::Polynomial,
};

/// AES state
#[derive(Debug, Clone)]
pub struct AESState {
    /// dat
    pub data: [u8; 16],
}

impl AESState {
    /// creates new instance of AES state
    pub fn new(data: [u8; 16]) -> Self {
        Self { data }
    }

    /// shift rows
    pub fn shift_rows(&mut self) {
        let mut temp = [0u8; 16];

        let s = &self.data;

        // Row 0 (no shift)
        temp[0] = s[0];
        temp[4] = s[4];
        temp[8] = s[8];
        temp[12] = s[12];

        // Row 1 (shift left by 1)
        temp[1] = s[5];
        temp[5] = s[9];
        temp[9] = s[13];
        temp[13] = s[1];

        // Row 2 (shift left by 2)
        temp[2] = s[10];
        temp[6] = s[14];
        temp[10] = s[2];
        temp[14] = s[6];

        // Row 3 (shift left by 3)
        temp[3] = s[15];
        temp[7] = s[3];
        temp[11] = s[7];
        temp[15] = s[11];

        self.data = temp;
    }

    /// sub bytes
    pub fn sub_bytes(&mut self) {
        for b in &mut self.data {
            let gf = GF256::from_byte(*b);
            *b = gf.sbox();
        }
    }
}

/// AES Modulus
#[derive(Debug)]
pub struct AESModulus;

impl Modulus<Fp<2>> for AESModulus {
    fn inner() -> Polynomial<Fp<2>> {
        // f(x) = x ^ 8 + x ^ 4 + x ^ 3 + x + 1
        Polynomial::new(vec![
            Fp::<2>::new(1),
            Fp::<2>::new(1),
            Fp::<2>::new(0),
            Fp::<2>::new(1),
            Fp::<2>::new(1),
            Fp::<2>::new(0),
            Fp::<2>::new(0),
            Fp::<2>::new(0),
            Fp::<2>::new(1),
        ])
    }
}

/// AES field
pub type GF256 = ExtensionField<Fp<2>, AESModulus>;

impl GF256 {
    /// constructs GF256 from byte
    pub fn from_byte(value: u8) -> Self {
        let bits = byte_to_bits_lsb(value).to_vec();

        let field_elements: Vec<Fp<2>> = bits.iter().map(|elem| to_field_element(*elem)).collect();

        Self::new(Polynomial::new(field_elements))
    }

    /// converts self to bytes
    pub fn to_byte(&self) -> u8 {
        let polynomial = self.polynomial();
        let cofficients = polynomial.cofficients();
        let mut bits = [0; 8];
        for i in 0..cofficients.len() {
            bits[i] = cofficients[i].value() as u8;
        }
        bits_to_byte_lsb(bits.as_slice())
    }

    /// affine transfromation
    pub fn sbox(&self) -> u8 {
        let y = if self.is_zero() {
            0
        } else {
            self.inv().unwrap().to_byte()
        };

        y ^ ((y << 1) | (y >> 7))
            ^ ((y << 2) | (y >> 6))
            ^ ((y << 3) | (y >> 5))
            ^ ((y << 4) | (y >> 4))
            ^ 0x63
    }

    /// mix column
    pub fn mix_column(&self) -> [u8; 4] {
        // self.mul(:)
        todo!()
    }
}

fn byte_to_bits_lsb(byte: u8) -> [u8; 8] {
    let mut bits = [0u8; 8];

    () = bits
        .iter_mut()
        .enumerate()
        .map(|(i, item)| {
            *item = (byte >> i) & 1;
        })
        .collect();
    bits
}

fn bits_to_byte_lsb(bits: &[u8]) -> u8 {
    let mut byte = 0u8;
    for (i, &bit) in bits.iter().enumerate().take(8) {
        // Shift the bit (0 or 1) to its correct position (2^i)
        // and OR it into the result
        byte |= (bit & 1) << i;
    }
    byte
}

fn to_field_element(value: u8) -> Fp<2> {
    Fp::<2>::new(value as u64)
}

mod tests {
    #[allow(unused_imports, dead_code)]
    use super::*;

    #[test]
    fn test_byte_conversion() {
        let b = 0x63;

        assert_eq!(GF256::from_byte(b).to_byte(), b)
    }

    #[test]
    fn test_inversion() {
        for b in 1..=255 {
            let a = GF256::from_byte(b);
            let inv = a.inv().unwrap();
            let prod = a.mul(&inv).unwrap();
            assert_eq!(prod.to_byte(), 1);
        }
    }

    #[test]
    fn test_sbox() {
        let a = GF256::from_byte(0x00);
        assert_eq!(a.sbox(), 0x63);

        let a = GF256::from_byte(0x53);
        assert_eq!(a.sbox(), 0xED);

        let a = GF256::from_byte(0x7C);
        assert_eq!(a.sbox(), 0x10);

        let a = GF256::from_byte(0x01);
        assert_eq!(a.sbox(), 0x7c);

        let a = GF256::from_byte(0xf0);
        assert_eq!(a.sbox(), 0x8c);

        let a = GF256::from_byte(0x9a);
        assert_eq!(a.sbox(), 0xb8);

        let a = GF256::from_byte(0x2b);
        assert_eq!(a.sbox(), 0xF1);

        let a = GF256::from_byte(0xFF);
        assert_eq!(a.sbox(), 0x16);

        let a = GF256::from_byte(0xFE);
        assert_eq!(a.sbox(), 0xBB);
    }

    #[test]
    fn test_shift_rows() {
        let mut state = AESState::new([
            0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d,
            0x0e, 0x0f,
        ]);

        state.shift_rows();

        assert_eq!(
            state.data,
            [
                0x00, 0x05, 0x0a, 0x0f, 0x04, 0x09, 0x0e, 0x03, 0x08, 0x0d, 0x02, 0x07, 0x0c, 0x01,
                0x06, 0x0b,
            ]
        );
    }
}
