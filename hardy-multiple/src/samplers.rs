use rand::Rng;

use genetic_architecture::{Chromosome, Gamete, GeneticArchitecture, Individual};

pub struct UniformSampler {

}

pub trait IndividualSampler {
    fn sample<R: Rng>(&mut self, genetic_architecture: &GeneticArchitecture, rng: &mut R) -> Individual;
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
