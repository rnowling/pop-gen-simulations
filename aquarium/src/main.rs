mod config;
mod fileio;
mod population;

use config::*;
use fileio::*;
use population::*;

fn main() {
    let params = SimParameters {
        n_individuals: 100,
        chromosome_length: 20,
        mutation_rate: 1e-3,
        population_initialization_strategy: PopulationInitializationStrategy::ClonedFromSingleIndividual
    };
    
    let mut sim = Simulation::new(params);
    
    sim.initialize();
    sim.print();
    println!();
    for _ in 1..1000 {
        sim.step();
        //sim.print();
        //println!();
    }

    sim.print();
    println!();

    write_genotypes("genotypes.tsv", sim.to_matrix().unwrap());
}
