///
/// Strategy for initializing the population
///
/// ClonedFromSingleIndividual: Founding individual has two empty chromosomes.  All
/// individuals are clones of the founder.
///
/// AllRandomIndividuals: For each individual, randomly generate a chromosome and
/// create two copies.
///
#[derive(Clone)]
pub enum PopulationInitializationStrategy {
    ClonedFromSingleIndividual,
    AllRandomIndividuals
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
    /// of mutations per site per generation
    pub mutation_rate: f64,

    /// strategy for creating initial individuals in the population
    pub population_initialization_strategy: PopulationInitializationStrategy,

    /// probability of a recombination given in a rate
    /// of recombination events per gamete per generation
    pub recombination_rate: f64,
}
