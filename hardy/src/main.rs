extern crate rand;

use rand::{Rng, thread_rng};
use rand::distributions::{Distribution, Range};
use std::collections::BTreeMap;
use std::cmp::Ordering;
use std::fmt;
use std::mem;

// Simulation of Hardy model
//
// Assumes a single diploid autosome with a single locus.
// There are two alleles, resulting in three possible genotypes.
// There is no mutation and the sex of the parents doesn't matter.
// Mating is random.

const HOMOZYGOUS: u8 = 2;


#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct Individual(u8, u8);

impl fmt::Display for Individual {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({}, {})", self.0, self.1)
    }
}

impl Ord for Individual {
    fn cmp(&self, other: &Individual) -> Ordering {
        if self.0 < other.0 {
            Ordering::Less
        } else if self.0 > other.0 {
            Ordering::Greater
        } else if self.1 < other.1 {
            Ordering::Less
        } else if self.1 > other.1 {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}

impl PartialOrd for Individual {
    fn partial_cmp(&self, other: &Individual) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// impl Individual {
//     pub fn new(allele1: u8, allele2: u8) -> Individual {
//         Individual(allele1, allele2)
//     }
// }

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

    pub fn genotype_frequencies(&self) -> BTreeMap<Individual, f64> {
        let mut counts: BTreeMap<Individual, f64> = BTreeMap::new();
    
        for genotype in &self.parents {
            let count = counts.entry(*genotype).or_insert(0f64);
            *count += 1.0f64;
        }

        let sum: f64 = counts.values().sum();

        for (_, count) in &mut counts {
            *count /= sum;
        }

        counts
    }

    pub fn swap_generations(&mut self) -> () {
        mem::swap(&mut self.parents, &mut self.children);
    }

    pub fn mate<R: Rng>(&mut self, rng: &mut R) -> () {
        let parent_range = Range::new(0usize, self.n_individuals);

        for mut child in self.children.iter_mut() {
            child.0 = 0;
            child.1 = 0;

            let parent1: &Individual = &self.parents[parent_range.sample(rng)];
            let parent2: &Individual = &self.parents[parent_range.sample(rng)];

            match parent1 {
                &Individual(HOMOZYGOUS, _) => child.0 = 1,
                &Individual(_, HOMOZYGOUS) => child.1 = 1,
                _ => if rng.gen::<bool>() {
                    child.0 = 1;
                } else {
                    child.1 = 1;
                }
            }

            match parent2 {
                &Individual(HOMOZYGOUS, _) => child.0 += 1,
                &Individual(_, HOMOZYGOUS) => child.1 += 1,
                _ => if rng.gen::<bool>() {
                    child.0 += 1;
                } else {
                    child.1 += 1;
                }
            }

        }

    }

}

fn main() {
    let population_size = 100000;
    let generations = 10;

    let mut rng = thread_rng();

    let mut population = Population::new(population_size,
                                     &mut rng);
    for i in 0..generations {
        println!("Population: {}", i);
        for (genotype, freq) in &population.genotype_frequencies() {
            println!("Genotype: {}, Frequency: {}", genotype, freq);
        }
        println!();
        
        population.mate(&mut rng);
        population.swap_generations();
    }    
}
