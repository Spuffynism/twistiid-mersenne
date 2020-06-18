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
const A: u32 = 0x9908B0DF;
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

// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
const DEFAULT_SEED: u32 = 5489;

struct Random {
    mt: [u32; N as usize],
    index: usize,
    lower_mask: u32,
    upper_mask: u32,
}

impl Random {
    fn new() -> Self {
        Self::from_seed(DEFAULT_SEED)
    }

    fn from_seed(seed: u32) -> Self {
        let lower_mask = (1 << R) - 1;
        let mut random = Random {
            mt: [0u32; N as usize],
            index: N as usize + 1,
            lower_mask,
            upper_mask: !lower_mask,
        };

        random.seed_mt(seed);

        random
    }

    fn seed_mt(&mut self, seed: u32) {
        self.index = N as usize;
        self.mt[0] = seed;
        for i in 1..(N as usize - 1) {
            let previous = self.mt[i - 1] as u64;
            let operation =
                (F as u64 * (previous ^ (previous >> ((W - 2) as u64))))
                    + (i as u64);

            self.mt[i] = operation as u32;
        }
    }

    fn next(&mut self) -> u32 {
        if self.index >= N as usize {
            self.twist();
        }

        let mut y = self.mt[self.index];
        y ^= (y >> U) & D;
        y ^= (y << S) & B;
        y ^= (y << T) & C;
        y ^= y >> L;

        self.index += 1;

        y
    }

    fn twist(&mut self) {
        for i in 0..(N as usize - 1) {
            let x = (self.mt[i] & self.upper_mask)
                + (self.mt[(i + 1) % N as usize] & self.lower_mask);
            let mut x_a = x >> 1;
            if x % 2 != 0 {
                x_a ^= A;
            }
            self.mt[i] = self.mt[(i + M as usize) % N  as usize] ^ x_a;
        }

        self.index = 0;
    }
}

#[cfg(test)]
mod tests {
    use crate::{DEFAULT_SEED, N, Random, R};

    const SOME_SEED: u32 = 32;

    fn some_generator() -> Random {
        Random::from_seed(SOME_SEED)
    }

    #[test]
    fn new_initializes_with_default_seed() {
        let generator = Random::new();

        let seed = generator.mt[0];

        assert_eq!(seed, DEFAULT_SEED);
    }

    #[test]
    fn from_seed_seeds_with_given_seed() {
        let generator = Random::from_seed(SOME_SEED);

        assert_eq!(generator.mt[0], SOME_SEED);
    }

    #[test]
    fn from_seed_sets_initial_fields() {
        let generator = Random::from_seed(SOME_SEED);

        assert_ne!(generator.mt[1], 0);
        assert_eq!(generator.index, N as usize);
        assert_eq!(generator.lower_mask, (1 << R) - 1);
        assert_eq!(!generator.lower_mask, generator.upper_mask);
    }

    #[test]
    fn from_seed_seeds() {
        let generator = Random::from_seed(SOME_SEED);

        assert_eq!(generator.index, N as usize);
        assert_eq!(generator.mt[0], SOME_SEED);
    }

    #[test]
    fn next_twists_when_needed() {
        let mut generator = some_generator();

        let index_after_twist = 0;
        let index_after_number_extract = 1;
        let index_with_no_twist = index_after_twist + index_after_number_extract;
        assert_ne!(generator.index, index_with_no_twist);

        let _ = generator.next();

        let index_after_first_twist = index_after_twist + index_after_number_extract;
        assert_eq!(generator.index, index_after_first_twist);

        let _ = generator.next();

        let index_after_second_twist = index_after_first_twist + index_after_number_extract;
        assert_eq!(generator.index, index_after_second_twist);
    }

    /// would panic if `twist()` was not called
    #[test]
    fn next_allows_generation_of_more_than_n_numbers() {
        let mut generator = some_generator();
        for i in 0..(N as usize ) * 2 {
            generator.next();
        }
    }

    #[test]
    fn next_generates_expected_number_with_seed() {
        let actual = some_generator().next();
        let expected = 3688901335;

        assert_eq!(expected, actual);
    }

    #[test]
    fn next_generates_same_number_with_same_seed() {
        let a_number = some_generator().next();

        let a_same_number = some_generator().next();

        assert_eq!(a_number, a_same_number);
    }

    #[test]
    fn next_generates_different_numbers_for_different_seeds() {
        let a_number = some_generator().next();

        let a_different_number = Random::from_seed(11).next();

        assert_ne!(a_number, a_different_number);
    }

    #[test]
    fn next_generates_different_numbers_on_sequential_calls() {
        let mut generator = Random::new();
        let a_number = generator.next();

        let a_different_number = generator.next();

        assert_ne!(a_number, a_different_number);
    }

    #[test]
    fn next_bumps_index() {
        let mut generator = Random::new();

        generator.twist();

        let previous_index = generator.index;

        let _ = generator.next();

        let current_index = generator.index;

        assert_eq!(previous_index + 1, current_index);
    }

    #[test]
    fn twist_resets_index() {
        let mut generator = Random::new();

        generator.twist();

        assert_eq!(generator.index, 0);
    }
}
