// in this script, genes can move around, be duplicated and deleted

use crate::mutation::MutationMap;
use std::collections::HashMap;
use crate::population::Genome;
use rand::rngs::StdRng;
use crate::mutation::Distribution as MutationDistribution;

// TODO need to think of way of rearranging each feature, and taking into account where insertions
// and translocations occur. Also need two different TE compartments, one which copies and inserts
// and one which moves around and inserts. 

// also need to think about how to determine whether a TE inserts into another gene, making it non-functional
// or whether it is upstream or downstream and can augment its function, having a multiplicative effect on its fitness contribution.


// for a given NucElement, store its position in the genome
// which can then be shuffled around by structural mutations, or copied

// hold probability of structural mutation for a given element.
// can then sample from uniform distribution to determine whether a structural mutation occurs, and if so, which one, and where it moves to.
#[derive(Clone)]
pub struct StructureMutationMap {
    pub duplication_rate: f64,
    pub deletion_rate: f64,
    pub inversion_rate: f64,
    pub recombination_rate: f64,
}

// write function which runs through each element and determines whether a structural mutation occurs, and if so, which one, and where it moves to.
pub fn mutate_intra_genome(genome: &mut Genome, mu_dist: &MutationDistribution, pos_dist: &MutationDistribution) {
    let mut thread_rng = rand::thread_rng();

    // For all intra genome comparisons, sample from uniform distribution to determine if variant occurs
    // and poisson distribution to determine where duplication goes

    // store hashmap of positions of genome elements, can store multiple per entry to capture duplications
    let new_positions: HashMap<usize, Vec<usize>> = HashMap::new(); // placeholder for new positions of each element after structural mutations, which will be used to update the mutation maps of each element after all structural mutations have been processed
    let genome_size = genome.seq.len();

    for (current_pos, element) in &mut genome.seq.iter().enumerate() {

        //store element structure positions, with current position first
        let mut element_structure_vec: Vec<i64> = vec![current_pos as i64];
        
        // duplications, can model multiple duplications repeatedly sampling until rand_val is above duplication rate
        let mut rand_val = mu_dist.sample(&mut thread_rng);
        while rand_val < element.structure_mutation_map.duplication_rate {
            // sample from poisson distribution to determine where duplication goes
            let duplication_pos = pos_dist.sample(&mut thread_rng) as i64;

            // determine if position is before or after gene, adjust if too large
            let pos_rand_val = mu_dist.sample(&mut thread_rng);
            let pos_order: i64 = if pos_rand_val < 0.5 { -1 } else { 1 };
            let mut new_pos = current_pos as i64 + (duplication_pos * pos_order);
            if new_pos < 0 {
                new_pos = 0;
            } else if new_pos >= genome_size as i64 {
                new_pos = genome_size as i64 - 1;
            }

            element_structure_vec.push(new_pos);
            rand_val = mu_dist.sample(&mut thread_rng);
        }

        // deletions, only first gene deleted
        rand_val = mu_dist.sample(&mut thread_rng);
        if rand_val < element.structure_mutation_map.deletion_rate {
            let _ = element_structure_vec.remove(0);
        }

        // translocations not explicitely modelled, as captured by simulatenous duplication and deletion event
        
        // inversions, apply to all genes in element structure vec
        for pos in element_structure_vec.iter_mut() {
            rand_val = mu_dist.sample(&mut thread_rng);
            if rand_val < element.structure_mutation_map.inversion_rate {
                *pos *= -1i64;
            }
        }       
    }
}