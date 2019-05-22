extern crate bit_set;
extern crate rand;

use bit_set::BitSet;

use rand::prelude::*;
use rand::seq::SliceRandom;

use std::collections::HashMap;

use super::parameters::*;

#[derive(Clone)]
pub struct Chromosome {
    pub alleles: BitSet,
    pub inverted: bool,
}

pub type Individual = [Chromosome; 2];

pub type Population = Vec<Individual>;

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
fn mutate(chrom: &mut Chromosome, params: &SimParameters) -> () {
    let mut rng = rand::thread_rng();
    
    // mutate bits
    for i in 0..params.chromosome_length {
        if rng.gen_bool(params.mutation_rate) {
            flip_bit(&mut chrom.alleles, i);
        }
    }
}

///
/// Perform recombination between a parent's chromosomes.
///
/// 1. If neither chromosome has mutations, return an empty chromosome
/// 3. Flip a coin to determine if recombination is happening
/// 3. If true, find a crossover position.  Take mutants from chromosome 1
///    for positions <= crossover position and chromosome 2 for positions
///    > crossover position.
/// 4. If false, randomly pick one of the chromosomes
fn recombine(parent: &Individual, params: &SimParameters) -> Chromosome {
    let mut rng = rand::thread_rng();

    let alleles1 = &parent[0].alleles;
    let alleles2 = &parent[1].alleles;

    // if there are no mutations, then return an empty chromosome
    let alleles = if alleles1.is_empty() && alleles2.is_empty() {
        BitSet::new()
    } else {
        // otherwise, flip a coin to determine if recombination is happening
        match rng.gen_bool(params.recombination_rate) {
            true => {
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
            },

            false => parent.choose(&mut rng).unwrap().alleles.clone()
        }
    };

    Chromosome { alleles: alleles, inverted: false }
}


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
fn reproduce(population: &Population, params: &SimParameters) -> Population {
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
fn randomly_generate_population(params: &SimParameters) -> Population {
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
fn clone_population(params: &SimParameters) -> Population {
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

///
/// Structure for capturing the simulation state.
///
#[derive(Clone)]
pub struct Simulation {
    pub params: SimParameters,
    pub current_generation: Option<Population>,
}

impl Simulation {
    ///
    /// Create a new simulation
    ///
    pub fn new(params: SimParameters) -> Simulation {
        Simulation {
            params: params,
            current_generation: None
        }
    }
    
    ///
    /// Initialize simulation by generating an initial population
    ///
    pub fn initialize(&mut self) -> () {
        println!("Initializing!");
        let strategy = &self.params.population_initialization_strategy;
        let population = match strategy {
            PopulationInitializationStrategy::ClonedFromSingleIndividual => clone_population(&self.params),
            PopulationInitializationStrategy::AllRandomIndividuals => randomly_generate_population(&self.params)
        };
        
        self.current_generation = Some(population);
    }

    ///
    /// Simulate one generation
    ///
    pub fn step(&mut self) -> () {
        self.current_generation = match self.current_generation {
            None => None,
            Some(ref g) => Some(reproduce(g, &self.params)),
        };
    }

    ///
    /// Identifies and removes positions where all samples are
    /// homozygous.
    ///
    pub fn trim(&mut self) -> () {
        let matrix_option = self.to_matrix();

        if matrix_option.is_none() {
            return;
        }

        let matrix = matrix_option.unwrap();
        let current_generation = self.current_generation.as_ref().unwrap();
        
        let mut fixed = BitSet::new();

        // find positions
        for (pos, genotypes) in matrix {
            let sum = genotypes.iter()
                .fold(0usize, |sum, val| sum + *val as usize);
                    
            // 2 chromosomes per individual
            if sum == 2 * self.params.n_individuals {
                fixed.insert(pos);
            }
        }

        // filter out mutations
        let mut trimmed = Vec::with_capacity(self.params.n_individuals);
        for ref indiv in current_generation {
            let mut chrom1 = indiv[0].clone();
            let mut chrom2 = indiv[1].clone();
            chrom1.alleles.difference_with(&fixed);
            chrom2.alleles.difference_with(&fixed);
            trimmed.push([chrom1, chrom2]);
        }
        
        self.current_generation = Some(trimmed);
    }

    ///
    /// Print out all individuals
    ///
    pub fn print(&mut self) -> () {
        match self.current_generation {
            Some(ref g) =>
                for ref indiv in g {
                    print!("[");
                    for i in indiv[0].alleles.iter() {
                        print!("{}, ", i);
                    }
                    print!("], ");
                    
                    print!("[");
                    for i in indiv[1].alleles.iter() {
                        print!("{}, ", i);
                    }
                    println!("] ");
                },
            None => println!("Uninitialized!")
        }
    }

    ///
    /// Convert mutations to a sparse matrix of genotypes.
    /// The keys of the resulting map are the positions, while
    /// the values are vectors of the individuals' genotypes.
    ///
    pub fn to_matrix(&self) -> Option<HashMap<usize, Vec<u8>>> {
        match self.current_generation {
            None => None,
            Some(ref g) => {
                let n_individuals = self.params.n_individuals;

                let mut matrix : HashMap<usize, Vec<u8>> = HashMap::new();

                for (idx, individual) in g.iter().enumerate() {
                    for pos in individual[0].alleles.iter() {
                        match matrix.get_mut(&pos) {
                            Some(ref mut genotype_counts) => genotype_counts[idx] += 1,
                            None => {
                                let mut genotype_counts = vec![0u8; n_individuals];
                                genotype_counts[idx] += 1;
                                matrix.insert(pos, genotype_counts);
                            }
                        }
                    }

                    for pos in individual[1].alleles.iter() {
                        match matrix.get_mut(&pos) {
                            Some(ref mut genotype_counts) => genotype_counts[idx] += 1,
                            None => {
                                let mut genotype_counts = vec![0u8; n_individuals];
                                genotype_counts[idx] += 1;
                                matrix.insert(pos, genotype_counts);
                            }
                        }
                    }
                }

                Some(matrix)
            }
        }
    }
}
