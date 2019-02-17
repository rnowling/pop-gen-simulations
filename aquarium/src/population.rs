extern crate bit_set;
extern crate rand;

use bit_set::BitSet;

use rand::prelude::*;

#[derive(Clone)]
pub struct Chromosome {
    pub alleles: BitSet,
    pub inverted: bool,
}

pub type Individual = (Chromosome, Chromosome);

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
fn run_mutation(chrom1: &mut Chromosome, chrom2: &mut Chromosome, params: &SimParameters) -> () {
    let mut rng = rand::thread_rng();
    
    // mutate bits
    for i in 0..params.chromosome_length {
        if rng.gen_bool(params.mutation_rate) {
            flip_bit(&mut chrom1.alleles, i);
        }

        if rng.gen_bool(params.mutation_rate) {
            flip_bit(&mut chrom2.alleles, i);
        }
    }
}


///
/// Produces a new individual from two parents. The basic steps include:
///
/// 1. Produce two gametes from each parent through recombination
/// 2. Select one gamate from each parent
/// 3. Mutate alleles at individual sites
///
/// This function returns an Option<Individual> as future versions may support
/// evaluating an individual for being non-viable.
///
fn mate(parent1: &Individual, parent2: &Individual, params: &SimParameters) -> Option<Individual> {
    let mut rng = rand::thread_rng();

    // TODO: recombination

    // choose chromosome pairs
    let mut chrom1 = match rng.gen_bool(0.5) {
        true => parent1.0.clone(),
        false => parent1.1.clone(),
    };
        
    let mut chrom2 = match rng.gen_bool(0.5) {
        true => parent2.0.clone(),
        false => parent2.1.clone(),
    };

    run_mutation(&mut chrom1, &mut chrom2, params);

    Some((chrom1, chrom2))
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
        let individual = (randomly_generate_chromosome(params),
                          randomly_generate_chromosome(params));
        individuals.push(individual);
    }

    individuals
}

///
/// Structure of simulation parameters.
///
#[derive(Clone)]
pub struct SimParameters {
    /// number of individuals in a single population
    pub n_individuals: usize,

    /// number of sites per chromosome
    pub chromosome_length: usize,

    /// probability of a mutation occuring given in a rate
    /// of mutations per site / per generation
    pub mutation_rate: f64,
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
            current_generation: None,
        }
    }
    
    ///
    /// Initialize simulation by generating an initial population
    ///
    pub fn initialize(&mut self) -> () {
        println!("Initializing!");
        let population = randomly_generate_population(&self.params);
        
        self.current_generation = Some(population);
    }

    ///
    /// Simulate one generation
    ///
    pub fn step(&mut self) -> () {
        println!("Stepping!");
        
        self.current_generation = match self.current_generation {
            None => None,
            Some(ref g) => Some(reproduce(g, &self.params)),
        };
    }

    ///
    /// Print out all individuals
    ///
    pub fn print(&mut self) -> () {
        match self.current_generation {
            Some(ref g) =>
                for ref indiv in g {
                    print!("[");
                    for i in indiv.0.alleles.iter() {
                        print!("{}, ", i);
                    }
                    println!("]");
                    
                    print!("[");
                    for i in indiv.1.alleles.iter() {
                        print!("{}, ", i);
                    }
                    println!("]");
                    println!();
                },
            None => println!("Uninitialized!")
        }
    }            
}
