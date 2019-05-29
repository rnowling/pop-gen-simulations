extern crate bit_set;
extern crate rand;

use bit_set::BitSet;
use rand::prelude::*;

use super::model::*;
use super::parameters::*;

///
/// Utility function for flipping a bit in a BitSet
///
fn flip_bit(alleles: &mut BitSet, i: usize) -> () {
    if alleles.contains(i) {
        alleles.remove(i);
    } else {
        alleles.insert(i);
    }
}

///
/// Simulates the mutation stage.  Every bit in the chromosome alleles is
/// tested for mutation according to the mutation_rate probabilities given to
/// the SimParameters.  If a mutation occurs, the bit is flipped.
///
pub fn mutate(chrom: &mut Chromosome, params: &SimParameters) -> () {
    let mut rng = rand::thread_rng();
    
    // mutate bits
    for i in 0..params.chromosome_length {
        if rng.gen_bool(params.mutation_rate) {
            flip_bit(&mut chrom.alleles, i);
        }
    }
}

trait MutationStrategy {
    fn new(params: SimParameters) -> Self;

    fn reproduce(self, gamete: &mut Chromosome) -> ();
}

// I should find the proper name for this
struct RandomMutationStrategy {
    params: SimParameters,
}

impl MutationStrategy for RandomMutationStrategy {
    fn new(params: SimParameters) -> RandomMutationStrategy {
        RandomMutationStrategy {
            params: params
        }
    }

    fn reproduce(self, gamete: &mut Chromosome) -> () {
        let mut rng = rand::thread_rng();
    
        // mutate bits
        for i in 0..self.params.chromosome_length {
            if rng.gen_bool(self.params.mutation_rate) {
                flip_bit(&mut gamete.alleles, i);
            }
        }
    }
}

