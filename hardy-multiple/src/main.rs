extern crate itertools;
extern crate rand;

use itertools::Itertools;
use rand::{Rng, thread_rng};
use rand::distributions::{Distribution, Range};
use std::collections::BTreeMap;
// use std::cmp::Ordering;
// use std::fmt;
use std::mem;

// Simulation of Hardy model
//
// Assumes a single diploid autosome with a single locus.
// There are two alleles, resulting in three possible genotypes.
// There is no mutation and the sex of the parents doesn't matter.
// Mating is random.

pub type Locus = u8; // each allele had a particular locus is represented as an integer
pub type Gamete = Vec<Locus>;
pub type Chromosome = Vec<Gamete>; // a haploid organism has a single gamete per chromosome, a diploid organism has two gametes per chromosome, and so forth
pub type Individual = Vec<Chromosome>;

fn individual_to_string(individual : &Individual) -> String {
    let mut i_str : Vec<String> = Vec::new();
    
    for c in individual {
        let mut c_str : Vec<String> = Vec::new();
        for g in c {
            let mut s = String::from("(");
            s = s + &g.iter().map( |a| a.to_string() ).join(",");
            s = s + ")";
            c_str.push(s);
        }
        
        let s = vec!["[", &c_str.join(", "), "]"].join("");
        i_str.push(s);
    }
    
    return i_str.join("\n");
}
    
#[derive(Clone)]
pub struct ChromosomeArchitecture {
    pub n_loci: usize,
    pub n_alleles_per_locus: Vec<u8>
}

impl ChromosomeArchitecture {
    pub fn new(n_loci: usize, n_alleles_per_locus: Vec<u8>) -> ChromosomeArchitecture {
        ChromosomeArchitecture {
            n_loci: n_loci,
            n_alleles_per_locus: n_alleles_per_locus
        }
    }
}

// Describe the organism
//
//
#[derive(Clone)]
pub struct GeneticArchitecture {
    pub n_chromosomes: usize,
    pub n_gametes: usize,
    pub chromosome_architectures: Vec<ChromosomeArchitecture>
}

impl GeneticArchitecture {
    pub fn new(n_gametes: usize) -> GeneticArchitecture {
        let chromosomes : Vec<ChromosomeArchitecture> = Vec::new();

        GeneticArchitecture {
            n_chromosomes: 0,
            n_gametes: n_gametes,
            chromosome_architectures: chromosomes
        }
    }

    pub fn add_chromosome(&mut self, n_loci: usize, n_alleles_per_locus: Vec<u8>) -> () {
        let chromosome = ChromosomeArchitecture::new(n_loci, n_alleles_per_locus);
        
        self.n_chromosomes += 1;
        self.chromosome_architectures.push(chromosome);
    }
}

// Describe a population (or deme) of individuals
//
//
struct Population {
    pub genetic_architecture: GeneticArchitecture,
    pub n_individuals: usize,
    pub parents: Vec<Individual>,
    pub children: Vec<Individual>
}

// #[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]

// impl fmt::Display for Locus {
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         write!(f, "({}, {})", self.0, self.1)
//     }
// }

// impl Ord for Locus {
//     fn cmp(&self, other: &Individual) -> Ordering {
//         if self.0 < other.0 {
//             Ordering::Less
//         } else if self.0 > other.0 {
//             Ordering::Greater
//         } else if self.1 < other.1 {
//             Ordering::Less
//         } else if self.1 > other.1 {
//             Ordering::Greater
//         } else {
//             Ordering::Equal
//         }
//     }
// }

// impl PartialOrd for Individual {
//     fn partial_cmp(&self, other: &Individual) -> Option<Ordering> {
//         Some(self.cmp(other))
//     }
// }

// impl Individual {
//     pub fn new(allele1: u8, allele2: u8) -> Individual {
//         Individual(allele1, allele2)
//     }
// }

pub trait IndividualSampler {
    fn sample<R: Rng>(&mut self, genetic_architecture: &GeneticArchitecture, rng: &mut R) -> Individual;
}

impl Population {
    pub fn new<R: Rng>(genetic_architecture: &GeneticArchitecture, n: usize, sampler: &mut UniformSampler, rng: &mut R) -> Population {
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
        
        // let sum: f64 = counts.values().sum();

        // for (_, count) in &mut counts {
        //     *count /= sum;
        // }

        // counts

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

struct UniformSampler {

}

impl UniformSampler {
    pub fn new() -> UniformSampler {
        UniformSampler {

        }
    }
}

impl IndividualSampler for UniformSampler {
    fn sample<R: Rng>(&mut self, gen_arch: &GeneticArchitecture, rng: &mut R) -> Individual {
        let mut individual: Individual = Vec::new();
        let n_gamete = gen_arch.n_gametes;

        for chrom_arch in gen_arch.chromosome_architectures.iter() {
            let mut chromosome: Chromosome = Vec::new();

            for _ in 0..n_gamete {
                let mut gamete: Gamete = Vec::new();

                for alleles_per_locus in &chrom_arch.n_alleles_per_locus {
                    let allele = rng.gen_range(0u8, *alleles_per_locus);
                    gamete.push(allele);
                }

                chromosome.push(gamete);
            }

            individual.push(chromosome);
        }

        individual
    }
}

fn main() {
    let population_size = 100000;
    let generations = 100;
    let n_gametes = 2usize;

    let mut rng = thread_rng();

    let mut genetic_architecture = GeneticArchitecture::new(n_gametes);

    // 2 chromosomes (each chromosome has 2 gametes)
    // each chromosome has 2 loci of 2 alleles
    genetic_architecture.add_chromosome(2usize, vec![2u8, 2u8]);
    genetic_architecture.add_chromosome(2usize, vec![2u8, 2u8]);

    let mut sampler = UniformSampler::new();
    let mut population = Population::new(&genetic_architecture,
                                         population_size,
                                         &mut sampler,
                                         &mut rng);
    for i in 0..generations {
        println!("Generation: {}", i);
        let genotype_counts = population.genotype_counts();
        for (c, chrom_counts) in genotype_counts.iter().enumerate() {
            println!("Chromosome: {}", c);
            for locus_counts in chrom_counts {
                for (genotype, count) in locus_counts {
                    let genotype_string = genotype.iter().map( |gt| gt.to_string() ).join(",");
                    println!("Genotype: ({}), Count: {}", genotype_string, count);
                }
            }

            println!();
        }
        
        // for indiv in &population.parents {
        //     println!("{}", individual_to_string(indiv));
        //     println!();
        // }
        
        population.mate(&mut rng);
        population.swap_generations();
    }    
}
