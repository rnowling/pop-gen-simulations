extern crate bit_set;
extern crate rand;

use bit_set::BitSet;
use rand::prelude::*;

use super::model::*;
use super::parameters::*;

///
/// Perform recombination between a parent's chromosomes.
///
/// 1. If neither chromosome has mutations, return an empty chromosome
/// 3. Flip a coin to determine if recombination is happening
/// 3. If true, find a crossover position.  Take mutants from chromosome 1
///    for positions <= crossover position and chromosome 2 for positions
///    > crossover position.
/// 4. If false, randomly pick one of the chromosomes
pub fn recombine(parent: &Individual, params: &SimParameters) -> Chromosome {
    let mut rng = rand::thread_rng();

    let alleles1 = &parent[0].alleles;
    let alleles2 = &parent[1].alleles;

    // if there are no mutations, then return an empty chromosome
    let alleles = if alleles1.is_empty() && alleles2.is_empty() {
        BitSet::new()
    } else if params.recombination_rate > 0.0 && rng.gen_bool(params.recombination_rate) {
        // otherwise, flip a coin to determine if recombination is happening
        let mut mutant_pos = alleles1.union(alleles2).collect::<Vec<usize>>();
        mutant_pos.sort();

        let mut gamete = BitSet::new();

        let crossover_pos = *mutant_pos.choose(&mut rng).unwrap();
        for pos in mutant_pos {
            if (pos <= crossover_pos && alleles1.contains(pos)) ||
                (pos > crossover_pos && alleles2.contains(pos))
            {
                gamete.insert(pos);
            }
        }

        gamete
    } else {
        parent.choose(&mut rng).unwrap().alleles.clone()
    };

    Chromosome { alleles: alleles, inverted: false }
}
