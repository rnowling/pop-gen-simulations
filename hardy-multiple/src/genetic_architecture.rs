use itertools::Itertools;

pub type Locus = u8; // each allele had a particular locus is represented as an integer
pub type Gamete = Vec<Locus>;
pub type Chromosome = Vec<Gamete>; // a haploid organism has a single gamete per chromosome, a diploid organism has two gametes per chromosome, and so forth
pub type Individual = Vec<Chromosome>;

pub fn individual_to_string(individual : &Individual) -> String {
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
