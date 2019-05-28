extern crate config;

// all modules need to be defined here
// to be seen by each other
mod fileio;
mod model;
mod mutation;
mod parameters;
mod population;
mod population_founding;
mod recombination;


use std::env;
use std::process;

use config::{Config, File};

use fileio::*;
use parameters::*;
use population::*;

fn print_usage() {

}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() == 1 {
        print_usage();
        process::exit(0);
    }

    let config_filename = &args[1];
    
    let mut config = Config::new();
    config.merge(File::with_name(config_filename).required(true)).unwrap();

    let n_steps = config.get::<usize>("n_steps").unwrap();
    let output_file = config.get::<String>("output_file").unwrap();
    
    let params = SimParameters {
        n_individuals: config.get::<usize>("n_individuals").unwrap(),
        chromosome_length: config.get::<usize>("chromosome_length").unwrap(),
        mutation_rate: config.get::<f64>("mutation_rate").unwrap(),
        recombination_rate: config.get::<f64>("recombination_rate").unwrap(),
        population_initialization_strategy: PopulationInitializationStrategy::ClonedFromSingleIndividual
    };
    
    let mut sim = Simulation::new(params);
    
    sim.initialize();
    sim.print();
    println!();
    for i in 1..n_steps {
        sim.step();
        if i % 100 == 0 {
            println!("Trimming!");
            sim.print();
            println!();
            sim.trim();
            sim.print();
            println!();
        }
    }

    sim.trim();
    sim.print();
    println!();

    write_genotypes(&output_file, sim.to_matrix().unwrap());
}
