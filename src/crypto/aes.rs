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

    /// helper function to multiply
    fn gf_mul(a: u8, b: u8) -> u8 {
        let x = GF256::from_byte(a);
        let y = GF256::from_byte(b);

        x.mul(&y).unwrap().to_byte()
    }

    fn mul_2(a: u8) -> u8 {
        Self::gf_mul(a, 0x02)
    }

    fn mul_3(a: u8) -> u8 {
        Self::gf_mul(a, 0x03)
    }

    fn mix_single_column(col: &mut [u8; 4]) {
        let a0 = col[0];
        let a1 = col[1];
        let a2 = col[2];
        let a3 = col[3];

        col[0] = Self::mul_2(a0) ^ Self::mul_3(a1) ^ a2 ^ a3;
        col[1] = a0 ^ Self::mul_2(a1) ^ Self::mul_3(a2) ^ a3;
        col[2] = a0 ^ a1 ^ Self::mul_2(a2) ^ Self::mul_3(a3);
        col[3] = Self::mul_3(a0) ^ a1 ^ a2 ^ Self::mul_2(a3);
    }

    /// mix columns
    pub fn mix_column(&mut self) {
        for c in 0..4 {
            let i = c * 4;

            let mut col = [
                self.data[i],
                self.data[i + 1],
                self.data[i + 2],
                self.data[i + 3],
            ];

            Self::mix_single_column(&mut col);

            self.data[i] = col[0];
            self.data[i + 1] = col[1];
            self.data[i + 2] = col[2];
            self.data[i + 3] = col[3];
        }
    }

    /// add round key
    pub fn add_round_key(&mut self, round_key: &[u8; 16]) {
        () = (0..16)
            .map(|i| {
                self.data[i] ^= round_key[i];
            })
            .collect();
    }
}

/// key schedule
pub mod key_schedule {
    use crate::crypto::aes::GF256;

    /// wofd
    pub type Word = [u8; 4];

    /// rotate word
    pub fn rot_word(word: Word) -> Word {
        [word[1], word[2], word[3], word[0]]
    }

    /// sub word
    pub fn sub_word(word: Word) -> Word {
        [
            GF256::from_byte(word[0]).sbox(),
            GF256::from_byte(word[1]).sbox(),
            GF256::from_byte(word[2]).sbox(),
            GF256::from_byte(word[3]).sbox(),
        ]
    }

    /// rcon
    pub fn rcon(round: usize) -> Word {
        assert!((1..=10).contains(&round));

        let mut value = GF256::from_byte(1);

        let two = GF256::from_byte(0x02);

        for _ in 1..round {
            value = value.mul(&two).unwrap();
        }

        [value.to_byte(), 0x00, 0x00, 0x00]
    }

    /// full g function
    pub fn g(word: Word, round: usize) -> Word {
        let mut w = rot_word(word);

        w = sub_word(w);

        let rc = rcon(round);

        for i in 0..4 {
            w[i] ^= rc[i]
        }

        w
    }

    /// round key
    pub fn round_key(words: &[Word], round: usize) -> [u8; 16] {
        let mut key = [0u8; 16];

        for i in 0..4 {
            let word = words[round * 4 + i];

            key[4 * i] = word[0];
            key[4 * i + 1] = word[1];
            key[4 * i + 2] = word[2];
            key[4 * i + 3] = word[3];
        }

        key
    }
}

/// key expansion
pub mod key_expansion {
    use crate::crypto::aes::key_schedule::{self, Word};

    fn xor_words(a: Word, b: Word) -> Word {
        let mut res = Word::default();
        for i in 0..4 {
            res[i] = a[i] ^ b[i];
        }

        res
    }

    /// key expansion
    pub fn key_expansion(key: [u8; 16]) -> Vec<Word> {
        const NK: usize = 4;
        const NB: usize = 4;
        const NR: usize = 10;

        const TOTAL_WORDS: usize = NB * (NR + 1);

        let mut w: Vec<Word> = vec![[0u8; 4]; TOTAL_WORDS];

        // First 4 words come directly from key
        for i in 0..NK {
            w[i] = [key[4 * i], key[4 * i + 1], key[4 * i + 2], key[4 * i + 3]];
        }

        // Generate remaining words
        for i in NK..TOTAL_WORDS {
            let mut temp = w[i - 1];

            if i % NK == 0 {
                temp = key_schedule::g(temp, i / NK);
            }

            w[i] = xor_words(w[i - NK], temp);
        }

        w
    }
}

/// encryption logic
pub mod aes_encrypt {
    use crate::crypto::aes::{AESState, key_expansion, key_schedule};

    /// encrypt block
    pub fn encrypt_block(plaintext: [u8; 16], key: [u8; 16]) -> [u8; 16] {
        let words = key_expansion::key_expansion(key);

        let mut state = AESState::new(plaintext);

        // Initial round
        let rk0 = key_schedule::round_key(&words, 0);
        state.add_round_key(&rk0);

        // Rounds 1–9
        for round in 1..10 {
            state.sub_bytes();

            state.shift_rows();

            state.mix_column();

            let rk = key_schedule::round_key(&words, round);

            state.add_round_key(&rk);
        }

        // Final round
        state.sub_bytes();

        state.shift_rows();

        let rk_last = key_schedule::round_key(&words, 10);

        state.add_round_key(&rk_last);

        state.data
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
    #![allow(unused_imports, dead_code)]
    use super::*;
    use crate::crypto::aes::key_schedule::Word;

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
    #[test]
    fn test_mix_single_column() {
        let mut column = [0xdb, 0x13, 0x53, 0x45];
        AESState::mix_single_column(&mut column);

        assert_eq!(&column, &[0x8e, 0x4d, 0xa1, 0xbc]);
    }

    #[test]
    fn test_add_round_key_basic() {
        let mut state = AESState::new([
            0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd,
            0xee, 0xff,
        ]);

        let round_key = [
            0xff, 0xee, 0xdd, 0xcc, 0xbb, 0xaa, 0x99, 0x88, 0x77, 0x66, 0x55, 0x44, 0x33, 0x22,
            0x11, 0x00,
        ];

        state.add_round_key(&round_key);

        assert_eq!(
            state.data,
            [
                0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                0xff, 0xff,
            ]
        );
    }

    #[test]
    fn test_add_round_key_zero_key() {
        let original = [
            0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e,
            0x0f, 0x10,
        ];

        let mut state = AESState::new(original);

        let zero_key = [0u8; 16];

        state.add_round_key(&zero_key);

        assert_eq!(state.data, original);
    }

    #[test]
    fn test_rot_word_basic() {
        let input: Word = [0x09, 0xcf, 0x4f, 0x3c];

        let result = key_schedule::rot_word(input);

        assert_eq!(result, [0xcf, 0x4f, 0x3c, 0x09]);
    }

    #[test]
    fn test_sub_word_known_example() {
        let input: Word = [0xcf, 0x4f, 0x3c, 0x09];

        let result = key_schedule::sub_word(input);

        assert_eq!(result, [0x8a, 0x84, 0xeb, 0x01]);
    }

    #[test]
    fn test_sub_word_manual() {
        let input: Word = [0x53, 0x7c, 0x01, 0xff];

        let result = key_schedule::sub_word(input);

        assert_eq!(
            result,
            [
                GF256::from_byte(0x53).sbox(),
                GF256::from_byte(0x7c).sbox(),
                GF256::from_byte(0x01).sbox(),
                GF256::from_byte(0xff).sbox(),
            ]
        );
    }

    #[test]
    fn test_rcon_round_5() {
        assert_eq!(key_schedule::rcon(5), [0x10, 0x00, 0x00, 0x00]);
    }

    #[test]
    fn test_rcon_round_10() {
        assert_eq!(key_schedule::rcon(10), [0x36, 0x00, 0x00, 0x00]);
    }

    #[test]
    fn test_g_round_1() {
        let input: Word = [0x09, 0xcf, 0x4f, 0x3c];

        let result = key_schedule::g(input, 1);

        assert_eq!(result, [0x8b, 0x84, 0xeb, 0x01]);
    }

    #[test]
    fn test_key_expansion_first_generated_word() {
        let key = [
            0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf,
            0x4f, 0x3c,
        ];

        let w = key_expansion::key_expansion(key);

        assert_eq!(w[4], [0xa0, 0xfa, 0xfe, 0x17]);
    }

    #[test]
    fn test_round_1_key_words() {
        let key = [
            0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf,
            0x4f, 0x3c,
        ];

        let w = key_expansion::key_expansion(key);

        assert_eq!(w[4], [0xa0, 0xfa, 0xfe, 0x17]);
        assert_eq!(w[5], [0x88, 0x54, 0x2c, 0xb1]);
        assert_eq!(w[6], [0x23, 0xa3, 0x39, 0x39]);
        assert_eq!(w[7], [0x2a, 0x6c, 0x76, 0x05]);
    }

    #[test]
    fn test_key_expansion_word_count() {
        let key = [0u8; 16];

        let w = key_expansion::key_expansion(key);

        assert_eq!(w.len(), 44);
    }

    #[test]
    fn test_encrypt_block_known_vector() {
        let key = [
            0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf,
            0x4f, 0x3c,
        ];

        let plaintext = [
            0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37,
            0x07, 0x34,
        ];

        let expected = [
            0x39, 0x25, 0x84, 0x1d, 0x02, 0xdc, 0x09, 0xfb, 0xdc, 0x11, 0x85, 0x97, 0x19, 0x6a,
            0x0b, 0x32,
        ];

        let ciphertext = aes_encrypt::encrypt_block(plaintext, key);

        assert_eq!(ciphertext, expected);
    }
}
