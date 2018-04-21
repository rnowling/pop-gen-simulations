extern crate itertools;
extern crate rand;

use itertools::Itertools;
use rand::thread_rng;

mod genetic_architecture;
mod population;
mod samplers;

// Simulation of Hardy model
//
// Assumes a single diploid autosome with a single locus.
// There are two alleles, resulting in three possible genotypes.
// There is no mutation and the sex of the parents doesn't matter.
// Mating is random.

fn main() {
    let population_size = 100000;
    let generations = 100;
    let n_gametes = 2usize;

    let mut rng = thread_rng();

    let mut genetic_architecture = genetic_architecture::GeneticArchitecture::new(n_gametes);

    // 2 chromosomes (each chromosome has 2 gametes)
    // each chromosome has 2 loci of 2 alleles
    genetic_architecture.add_chromosome(2usize, vec![2u8, 2u8]);
    genetic_architecture.add_chromosome(2usize, vec![2u8, 2u8]);

    let mut sampler = samplers::UniformSampler::new();
    let mut population = population::Population::new(&genetic_architecture,
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
