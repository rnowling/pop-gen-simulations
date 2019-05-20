use std::collections::HashMap;

use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;

pub fn write_genotypes(filename: &str, matrix: HashMap<usize, Vec<u8>>) -> () {
    let file = File::create(filename).unwrap();
    let mut file = BufWriter::new(file);

    let mut sorted_keys = matrix.keys().collect::<Vec<&usize>>();
    sorted_keys.sort();
    for pos in sorted_keys {
        write!(file, "{}", pos + 1usize);
        for gt_count in matrix.get(pos).unwrap().iter() {
            write!(file, " {}", gt_count);
        }
        write!(file, "\n");
    }   
}
