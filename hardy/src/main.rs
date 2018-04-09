extern crate rand;

use rand::{Rng, thread_rng};
use rand::distributions::{Distribution, Range};
use std::collections::HashMap;
use std::mem;

// Simulation of Hardy model
//
// Assumes a single diploid autosome with a single locus.
// There are two alleles, resulting in three possible genotypes.
// There is no mutation and the sex of the parents doesn't matter.
// Mating is random.

#[derive(Debug, PartialEq, Eq, Hash)]
struct Individual(u8, u8);
impl Individual {
    pub fn new(allele1: u8, allele2: u8) -> Individual {
        Individual(allele1, allele2)
    }
}

struct Population {
    pub n_individuals: usize,
    pub parents: Vec<Individual>,
    pub children: Vec<Individual>
}

impl Population {
    pub fn new<R: Rng>(n: usize, rng: &mut R) -> Population {
        let genotypes = Range::new_inclusive(0usize, 2usize);
        let mut parents: Vec<Individual> = Vec::new();
        let mut children: Vec<Individual> = Vec::new();

        for _ in 0..n {
            let genotype_1 = match genotypes.sample(rng) {
                0 => Individual(2u8, 0u8),
                1 => Individual(1u8, 1u8),
                2 => Individual(0u8, 2u8),
                _ => Individual(0u8, 0u8)
            };
            parents.push(genotype_1);

            let genotype_2 = match genotypes.sample(rng) {
                0 => Individual(2u8, 0u8),
                1 => Individual(1u8, 1u8),
                2 => Individual(0u8, 2u8),
                _ => Individual(0u8, 0u8)
            };
            children.push(genotype_2);
        }

        Population {
            n_individuals: n,
            parents: parents,
            children: children
        }
    }

    pub fn count_genotypes(&mut self) -> HashMap<Individual, u32> {
        let mut counts: HashMap<Individual, u32> = HashMap::new();
    
        for individual in self.parents {
            let count = counts.entry(individual).or_insert(0u32);
            *count += 1;
        }

        counts
    }

    pub fn swap_generations(&mut self) -> () {
        mem::swap(&mut self.parents, &mut self.children);
    }
}

const HOMOZYGOUS: u8 = 2;

fn mate<R: Rng>(population: &mut Population, rng: &mut R) -> () {
    let parent_range = Range::new(0usize, population.n_individuals);
    let allele_range = Range::new_inclusive(0usize, 1usize);

    for child in population.children {
        child.0 = 0;
        child.1 = 0;
        
        let parent1 = population.parents[parent_range.sample(rng)];
        let parent2 = population.parents[parent_range.sample(rng)];

        if parent1.0 == HOMOZYGOUS {
            child.0 = 1;
        } else if parent1.1 == HOMOZYGOUS {
            child.1 = 1;
        } else {
            if allele_range.sample(rng) == 0 {
                child.0 = 1;
            } else {
                child.1 = 1;
            }
        }

        if parent2.0 == HOMOZYGOUS {
            child.0 += 1;
        } else if parent2.1 == HOMOZYGOUS {
            child.1 += 1;
        } else {
            if allele_range.sample(rng) == 0 {
                child.0 += 1;
            } else {
                child.1 += 1;
            }
        }
    }
}


fn main() {
    println!("Hello, world!");

    let population_size = 10000;
    let generations = 10;

    let mut rng = thread_rng();

    let population = Population::new(population_size,
                                     &mut rng);

    for _i in 0..generations {
        mate(&mut population, &mut rng);

       population.swap_generations();
    }    
}
