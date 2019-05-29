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
/// Strategy for mating to produce a new generation
///
/// DeterministicCloning: Create one clone of each individual.  (Only
/// mutation and recombination will perturb individuals.)
///
/// RandomSamplingCloning: Create a new generation by sampling with replacement.
/// Multiple copies will be created from some individuals, while no other individuals
/// will be lost.
///
/// RandomMating: New individuals will be created by randomly choosing two parents.
/// Parents can be the same individual so the child may be a clone.
///
#[derive(Clone)]
pub enum MatingStrategy {
    DeterministicCloning,
    RandomSamplingCloning,
    RandomMating
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

    /// strategy for choosing mates
    pub mating_strategy: MatingStrategy
}
