extern crate bit_set;

use bit_set::BitSet;

#[derive(Clone)]
pub struct Chromosome {
    pub alleles: BitSet,
    pub inverted: bool,
}

pub type Individual = [Chromosome; 2];

pub type Population = Vec<Individual>;
