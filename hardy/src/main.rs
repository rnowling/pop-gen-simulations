use std::collections::HashMap;

struct Locus {
    alleles: Vec<String>,
    mutation_rate: f64,
    genotype_frequencies: HashMap<String, f64>
}

enum SystemOfMating {
    Random,
}

struct Model {
    population_size: u32,
    mating: SystemOfMating,
    genetic_architecture: Vec<Locus>,
}

impl Model {
    fn empty(population_size: u32, mating: SystemOfMating, genetic_architecture: Vec<Locus>) -> Model {
        Model {
            population_size: population_size,
            mating: mating,
            genetic_architecture : genetic_architecture
        }
    }
}

fn main() {
    println!("Hello, world!");

    let mut loci: Vec<Locus> = Vec::new();
    let mut count = 0;

    while count < 2 {
        let mut genotype_frequencies = HashMap::new();

        genotype_frequencies.insert(String::from("aa"), 0.25);
        genotype_frequencies.insert(String::from("Aa"), 0.50);
        genotype_frequencies.insert(String::from("AA"), 0.25);
        
        let locus = Locus {
            alleles : vec![String::from("a"), String::from("A")],
            mutation_rate : 0.0f64,
            genotype_frequencies : genotype_frequencies,
        };

        loci.push(locus);

        count += 1;
    }
    
    
    let model = Model::empty(100u32,
                             SystemOfMating::Random,
                             loci);

    println!("The model's population has {} members.", model.population_size);
}
