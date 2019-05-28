extern crate bit_set;
extern crate rand;

use rand::prelude::*;

use super::model::*;
use super::mutation::*;
use super::parameters::*;
use super::recombination::*;

///
/// Produces a (haploid) gamete from a (diploid) parent. The basic steps include:
///
/// 1. Generate a new gamete through recombination
/// 2. Mutate alleles at individual sites
///
fn create_gamete(parent: &Individual, params: &SimParameters) -> Chromosome {
    let mut gamete = recombine(parent, params);

    mutate(&mut gamete, params);

    gamete
}

///
/// Mates two parents to produce an offspring.  The basic steps include:
///
/// 1. A gamete is created for each parent
/// 2. The individual is formed and returned
///
/// This function returns an Option<Individual> as future versions may support
/// evaluating an individual for being non-viable.
///
fn mate(parent1: &Individual, parent2: &Individual, params: &SimParameters) -> Option<Individual> {
    let chrom1 = create_gamete(parent1, params);
    let chrom2 = create_gamete(parent2, params);

    Some([chrom1, chrom2])
}

///
/// Produce a new generation of individuals for a population from an existing
/// generation.  The generations will have the same number of individuals.
///
pub fn reproduce(population: &Population, params: &SimParameters) -> Population {
    let size = population.len();
    let mut next_generation: Vec<Individual> = Vec::with_capacity(size);
    let mut rng = rand::thread_rng();
    
    while next_generation.len() < size {
        let idx1 = rng.gen_range(0, size);
        let idx2 = rng.gen_range(0, size);

        let parent1 = &population[idx1];
        let parent2 = &population[idx2];
        let child = mate(parent1, parent2, params);

        match child {
            Some(c) => next_generation.push(c),
            None => {},
        }
    }

    next_generation
}
