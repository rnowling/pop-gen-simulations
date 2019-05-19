mod population;

fn main() {
    println!("Hello, world!");

    let params = population::SimParameters {
        n_individuals: 10,
        chromosome_length: 20,
        mutation_rate: 1e-3,
        population_initialization_strategy: population::PopulationInitializationStrategy::ClonedFromSingleIndividual
    };
    
    let mut sim = population::Simulation::new(params);
    
    sim.initialize();
    sim.print();
    for _ in 1..40 {
        sim.step();
        sim.print();
    }
    
}
