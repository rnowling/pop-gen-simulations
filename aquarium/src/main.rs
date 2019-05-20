mod config;
mod population;

use config::*;
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
    for _ in 1..100 {
        sim.step();
        sim.print();
        println!();
    }
}
