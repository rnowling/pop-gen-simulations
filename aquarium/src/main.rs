mod population;

fn main() {
    let params = population::SimParameters {
        n_individuals: 100,
        chromosome_length: 20,
        mutation_rate: 1e-3,
        population_initialization_strategy: population::PopulationInitializationStrategy::ClonedFromSingleIndividual
    };
    
    let mut sim = population::Simulation::new(params);
    
    sim.initialize();
    sim.print();
    println!();
    for _ in 1..10 {
        sim.step();
        sim.print();
        println!();
    }
}
