#![allow(non_upper_case_globals)]

// word size (in number of bits)
const W: u32 = 32;
// degree of recurrence
const N: u32 = 624;
// middle word, an offset used in the recurrence relation defining the series x, 1 ≤ m < n
const M: u32 = 397;
// separation point of one word, or the number of bits of the lower bitmask, 0 ≤ r ≤ w - 1
const R: u32 = 31;
// coefficients of the rational normal form twist matrix
const A: u64 = 0x9908B0DF;
// b, c: TGFSR(R) tempering bitmasks
const B: u32 = 0x9D2C5680;
const C: u32 = 0xEFC60000;
// s, t: TGFSR(R) tempering bit shifts
const S: u32 = 7;
const T: u32 = 15;
// u, d, l: additional Mersenne Twister tempering bit shifts/masks
const U: u32 = 11;
const D: u32 = 0xFFFFFFFF;
const L: u32 = 18;

const F: u32 = 1812433253;

struct Random {
    mt: [u64; N as usize],
    index: usize,
    lower_mask: u32,
    upper_mask: u32,
}

impl Random {
    fn new() -> Self {
        let lower_mask = (1 << R) - 1;
        Random {
            mt: [0u64; N as usize],
            index: N as usize + 1,
            lower_mask,
            upper_mask: !lower_mask,
        }
    }

    fn rand(&mut self, seed: u32) -> u32 {
        self.seed_mt(seed);
        self.twist();

        self.extract_number()
    }

    fn seed_mt(&mut self, seed: u32) {
        self.index = N as usize;
        self.mt[0] = seed as u64;
        for i in 1..(N as usize - 1) {
            self.mt[i] = ((F as u64 * (self.mt[i - 1] ^ (self.mt[i - 1] >> ((W as u64) - 2))) + (i as u64)) as u32) as u64;
        }
    }

    fn extract_number(&mut self) -> u32 {
        let mut y = self.mt[self.index as usize];
        y ^= (y >> U as u64) & D as u64;
        y ^= (y << S as u64) & B as u64;
        y ^= (y << T as u64) & C as u64;
        y ^= y >> L as u64;

        self.index += 1;

        y as u32
    }

    fn twist(&mut self) {
        for i in 0..(N - 1) {
            let x = (self.mt[i as usize] & self.upper_mask  as u64)
                + (self.mt[((i + 1) % N) as usize] & self.lower_mask as u64);
            let mut x_a = x >> 1;
            if x % 2 != 0 {
                x_a ^= A;
            }
            self.mt[i as usize] = self.mt[((i + M) % N) as usize] ^ x_a;
        }

        self.index = 0;
    }
}

#[cfg(test)]
mod tests {
    use crate::Random;

    #[test]
    fn generates_same_number_with_same_seed() {
        let seed = 10;

        let a_number = Random::new().rand(seed);

        let a_same_number = Random::new().rand(seed);

        assert_eq!(a_number, a_same_number);
    }
    #[test]
    fn generates_different_numbers_for_different_seeds() {
        let a_number = Random::new().rand(10);

        let a_different_number = Random::new().rand(11);

        assert_ne!(a_number, a_different_number);
    }
}
