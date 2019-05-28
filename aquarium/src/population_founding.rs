extern crate bit_set;
extern crate rand;

use bit_set::BitSet;
use rand::prelude::*;

// how does this differ from
// use crate::model::*; ?

use super::model::*;
use super::parameters::*;

///
/// Generate a new chromosome by randomly generating an allele for each site.
///
fn randomly_generate_chromosome(params: &SimParameters) -> Chromosome {
    let mut rng = rand::thread_rng();

    let mut alleles = BitSet::new();
    for i in 0..params.chromosome_length {
        if rng.gen_bool(0.5) {
            alleles.insert(i);
        }
    }

    Chromosome {
        alleles: alleles,
        inverted: false
    }
}

///
/// Generate a population of individuals by randomly generating two chromosomes
/// for each individual.
///
pub fn randomly_generate_population(params: &SimParameters) -> Population {
    let mut individuals: Vec<Individual> = Vec::new();
    for _ in 0..params.n_individuals {
        let chrom = randomly_generate_chromosome(params);
        let individual = [chrom.clone(), chrom.clone()];
        individuals.push(individual);
    }

    individuals
}

///
/// Clone a population from a single randomly-generated individual.
///
pub fn clone_population(params: &SimParameters) -> Population {
    let mut individuals: Vec<Individual> = Vec::with_capacity(params.n_individuals);

    let chrom = Chromosome {
        alleles: BitSet::new(),
        inverted: false
    };

    let individual = [chrom.clone(), chrom.clone()];

    for _ in 0..params.n_individuals {
        individuals.push(individual.clone());
    }

    individuals
}
