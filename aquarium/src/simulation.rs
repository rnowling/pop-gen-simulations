extern crate bit_set;
extern crate rand;

use bit_set::BitSet;
use std::collections::HashMap;

// how does this differ from
// use crate::model::*; ?

use super::mating::*;
use super::model::*;
use super::parameters::*;
use super::population_founding::*;

///
/// Structure for capturing the simulation state.
///
#[derive(Clone)]
pub struct Simulation {
    params: SimParameters,
    current_generation: Option<Population>,
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
