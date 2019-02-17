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

pub fn mate(parent1: &Individual, parent2: &Individual, params: &SimParameters) -> Option<Individual> {
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

    // mutate bits
    for i in 0..params.chromosome_length {
        if rng.gen_bool(params.mutation_rate) {
            if chrom1.alleles.contains(i) {
                chrom1.alleles.remove(i);
            } else {
                chrom1.alleles.insert(i);
            }
        }

        if rng.gen_bool(params.mutation_rate) {
            if chrom2.alleles.contains(i) {
                chrom2.alleles.remove(i);
            } else {
                chrom2.alleles.insert(i);
            }
        }
    }

    Some((chrom1, chrom2))
}

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

pub fn randomly_generate_chromosome(params: &SimParameters) -> Chromosome {
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

pub fn randomly_generate_population(params: &SimParameters) -> Population {
    let mut individuals: Vec<Individual> = Vec::new();
    for _ in 0..params.n_individuals {
        let individual = (randomly_generate_chromosome(params),
                          randomly_generate_chromosome(params));
        individuals.push(individual);
    }

    individuals
}

#[derive(Clone)]
pub struct SimParameters {
    pub n_individuals: usize,
    pub chromosome_length: usize,
    pub mutation_rate: f64
}

#[derive(Clone)]
pub struct Simulation {
    pub params: SimParameters,
    pub current_generation: Option<Population>,
}

impl Simulation {
    pub fn new(params: SimParameters) -> Simulation {
        Simulation {
            params: params,
            current_generation: None,
        }
    }

    pub fn initialize(&mut self) -> () {
        println!("Initializing!");
        let population = randomly_generate_population(&self.params);
        
        self.current_generation = Some(population);
    }

    pub fn step(&mut self) -> () {
        println!("Stepping!");
        
        self.current_generation = match self.current_generation {
            None => None,
            Some(ref g) => Some(reproduce(g, &self.params)),
        };
    }

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
