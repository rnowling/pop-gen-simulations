use std::collections::BTreeMap;
use std::mem;

use rand::Rng;
use rand::distributions::{Distribution, Range};

use genetic_architecture::{GeneticArchitecture, Individual};
use samplers::IndividualSampler;

// Describe a population (or deme) of individuals
//
//
pub struct Population {
    pub genetic_architecture: GeneticArchitecture,
    pub n_individuals: usize,
    pub parents: Vec<Individual>,
    pub children: Vec<Individual>
}


impl Population {
    pub fn new<R: Rng, S: IndividualSampler>(genetic_architecture: &GeneticArchitecture, n: usize, sampler: &mut S, rng: &mut R) -> Population {
        let mut parents: Vec<Individual> = Vec::new();
        let mut children: Vec<Individual> = Vec::new();

        for _ in 0..n {
            parents.push(sampler.sample(genetic_architecture, rng));
            children.push(sampler.sample(genetic_architecture, rng));
        }

        Population {
            genetic_architecture: genetic_architecture.clone(),
            n_individuals: n,
            parents: parents,
            children: children
        }
    }

    pub fn genotype_counts(&self) -> Vec<Vec<BTreeMap<Vec<u8>, f64>>> {
        let n_chromosomes = self.genetic_architecture.n_chromosomes;
        let mut chromosome_counts: Vec<Vec<BTreeMap<Vec<u8>, f64>>> = Vec::with_capacity(n_chromosomes);

        let n_gametes = self.genetic_architecture.n_gametes;
        for (c, chrom_arch) in self.genetic_architecture.chromosome_architectures.iter().enumerate() {
            let n_loci = chrom_arch.n_loci;
            let mut loci_genotype_counts: Vec<BTreeMap<Vec<u8>, f64>> = Vec::with_capacity(n_loci);

            for l in 0..n_loci {
                let mut genotype_counts: BTreeMap<Vec<u8>, f64> = BTreeMap::new();
                let n_alleles = chrom_arch.n_alleles_per_locus[l];

                for individual in &self.parents {
                    let mut genotype: Vec<u8> = vec![0u8; n_alleles as usize];
                    
                    for g in 0..n_gametes {
                        let allele = individual[c][g][l] as usize;
                        genotype[allele] += 1u8;
                    }

                    
                    let count = genotype_counts.entry(genotype).or_insert(0f64);
                    *count += 1.0f64;
                }

                loci_genotype_counts.push(genotype_counts);
            }
            
            chromosome_counts.push(loci_genotype_counts);
        }
        
        chromosome_counts
    }

    pub fn swap_generations(&mut self) -> () {
        mem::swap(&mut self.parents, &mut self.children);
    }

    pub fn mate<R: Rng>(&mut self, rng: &mut R) -> () {
        let n_gametes = self.genetic_architecture.n_gametes;

        let parent_range = Range::new(0usize, self.n_individuals);
        let gamete_range = Range::new(0usize, n_gametes);

        for mut child in self.children.iter_mut() {
            let mut parents: Vec<&Individual> = Vec::new();

            for _ in 0..n_gametes {
                parents.push(&self.parents[parent_range.sample(rng)]);
            }
            
            for (i, chromosome) in child.iter_mut().enumerate() {
                chromosome.clear();

                for parent in parents.iter() {
                    let parent_chromosome = &parent[i];

                    // this is probably causing a copy...
                    chromosome.push(parent_chromosome[gamete_range.sample(rng)].clone());
                }
            }
        }
    }
}
