extern crate bit_set;
extern crate rand;

use rand::prelude::*;
use std::cmp::max;

use super::model::*;
use super::mutation::*;
use super::parameters;
use super::parameters::SimParameters;
use super::recombination::*;



struct DeterministicCloningMatingStrategy {
}

impl DeterministicCloningMatingStrategy {
    fn new() -> DeterministicCloningMatingStrategy {
        DeterministicCloningMatingStrategy {}
    }
}


struct RandomSamplingCloningMatingStrategy {
}

impl RandomSamplingCloningMatingStrategy {
    fn new() -> RandomSamplingCloningMatingStrategy {
        RandomSamplingCloningMatingStrategy {}
    }
}

struct RandomMatingStrategy {
}

impl RandomMatingStrategy {
    fn new() -> RandomMatingStrategy {
        RandomMatingStrategy {
        }
    }
}

/// Mating Strategy logic

trait MatingStrategy {
    fn choose_mates(self, population: &Population) -> Iterator<Item = [Individual; 2]>;
}

struct DeterministicCloningIterator<'a> {
    idx: usize,
    population: &'a Population
}

impl<'a> Iterator for DeterministicCloningIterator<'a> {
    type Item = [Individual; 2];

    fn next(&mut self) -> Option<Self::Item> {
        let parent = self.population[self.idx];
        
        self.idx = max(self.idx + 1, self.population.len() - 1);

        Some([parent, parent])
    }
}

impl MatingStrategy for DeterministicCloningMatingStrategy {
    fn choose_mates(self, population: &Population) -> impl Iterator<Item = [Individual; 2]> {
        DeterministicCloningIterator {
            idx: 0,
            population: population
        }
    }
}

struct RandomSamplingCloningIterator<'a> {
    rng: &'a mut ThreadRng,
    population: &'a Population
}

impl<'a> Iterator for RandomSamplingCloningIterator<'a> {
    type Item = [Individual; 2];

    fn next(&mut self) -> Option<Self::Item> {
        let size = self.population.len();
        let idx = self.rng.gen_range(0, size);
        let parent = self.population[idx];
        
        Some([parent, parent])
    }
}

impl MatingStrategy for RandomSamplingCloningMatingStrategy {
    fn choose_mates(self, population: &Population) -> impl Iterator<Item = [Individual; 2]> {
        RandomSamplingCloningIterator {
            rng: &mut rand::thread_rng(),
            population: population
        }
    }
}

struct RandomMatingIterator<'a> {
    rng: &'a mut ThreadRng,
    population: &'a Population
}

impl<'a> Iterator for RandomMatingIterator<'a> {
    type Item = [Individual; 2];

    fn next(&mut self) -> Option<Self::Item> {
        let size = self.population.len();
        let idx1 = self.rng.gen_range(0, size);
        let idx2 = self.rng.gen_range(0, size);
        let parent1 = self.population[idx1];
        let parent2 = self.population[idx2];
        
        Some([parent1, parent2])
    }
}

impl MatingStrategy for RandomMatingStrategy {
    fn choose_mates(self, population: &Population) -> impl Iterator<Item = [Individual; 2]> {
        RandomMatingIterator {
            rng: &mut rand::thread_rng(),
            population: population
        }
    }
}

pub fn mating_strategy_factory(params: SimParameters) -> Box<dyn MatingStrategy> {
    match params.mating_strategy {
        parameters::MatingStrategy::RandomMating =>
            Box::new(RandomMatingStrategy::new()),
        parameters::MatingStrategy::RandomSamplingCloning =>
            Box::new(RandomSamplingCloningMatingStrategy::new()),
        parameters::MatingStrategy::DeterministicCloning =>
            Box::new(DeterministicCloningMatingStrategy::new())
    }
}
