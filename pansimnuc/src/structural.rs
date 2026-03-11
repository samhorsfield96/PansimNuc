// in this script, genes can move around, be duplicated and deleted

use crate::{mutation::MutationMap, population::NucElement};
use std::collections::HashMap;
use crate::population::{Genome, Population};
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
    /// Cap on the number of duplications per element per generation.
    /// `None` means no cap (behaviour unchanged from before this field was added).
    pub max_duplications: Option<usize>,
}

// write function which runs through each element and determines whether a structural mutation occurs, and if so, which one, and where it moves to.
pub fn mutate_intra_genome(genome: &mut Genome, homology_map: &mut Vec<Vec<Vec<usize>>>, mu_dist: &MutationDistribution, pos_dist: &MutationDistribution) {
    let mut thread_rng = rand::thread_rng();

    // For all intra genome comparisons, sample from uniform distribution to determine if variant occurs
    // and poisson distribution to determine where duplication goes

    // store hashmap of positions of genome elements, can store multiple per entry to capture duplications
    let mut new_positions: HashMap<usize, Vec<i64>> = HashMap::new(); // placeholder for new positions of each element after structural mutations, which will be used to update the mutation maps of each element after all structural mutations have been processed
    let genome_size = genome.seq.len();
    let mut max_position: usize = 0;

    for (current_pos , element) in &mut genome.seq.iter().enumerate() {
        //store element structure positions, with current position first
        let mut element_structure_vec: Vec<i64> = vec![current_pos as i64];
        
        // duplications, can model multiple duplications repeatedly sampling until rand_val is above duplication rate
        let mut rand_val = mu_dist.sample(&mut thread_rng);
        let mut dup_count: usize = 0;
        while rand_val < element.structure_mutation_map.duplication_rate
            && element.structure_mutation_map.max_duplications.map_or(true, |max| dup_count < max)
        {
            dup_count += 1;
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

            // determine maximum position present
            if new_pos as usize > max_position {
                max_position = new_pos as usize;
            }
            if current_pos > max_position {
                max_position = current_pos;
            }
        }

        // deletions, only first gene deleted
        rand_val = mu_dist.sample(&mut thread_rng);
        if rand_val < element.structure_mutation_map.deletion_rate {
            let _ = element_structure_vec.remove(0);
        }

        // translocations not explicitely modelled, as captured by simulatenous duplication and deletion event
        
        // inversions, apply to all genes in element structure vec
        for pos in element_structure_vec.iter_mut() {
            
            // currently zero indexed, need to make 1 indexed so inversion logic works
            *pos += 1i64;

            rand_val = mu_dist.sample(&mut thread_rng);
            if rand_val < element.structure_mutation_map.inversion_rate {
                *pos *= -1i64;
            }
        }
        max_position += 1usize; // account for 1 indexing

        new_positions.insert(current_pos, element_structure_vec); 
    }

    // generate new genome based on all intra-genome variation, current everything is 1-indexed
    let mut new_genome_seq: Vec<i64> = vec![0; max_position];

    for idx in 0..genome_size {
        if let Some(prev_pos) = new_positions.get(&idx) {
            // value exists
            for entry in prev_pos { 
                if *entry == 0 {
                    panic!("Entry for structural varient is 0, for position {} in genome {}", idx, genome.identifier);
                }

                let insertion_pos = entry.abs() as usize - 1; // convert back to 0 indexed
                new_genome_seq[insertion_pos] = *entry;
            }
        } else {
            // value doesn't exist
            panic!("No structural entry for position {} in genome {}", idx, genome.identifier);
        }
    }

    // go through and delete any positions with value 0
    new_genome_seq.retain(|&x| x != 0);

    // generate new genome
    let mut new_genome: Vec<NucElement> = Vec::new();

    for element_idx in new_genome_seq {
        let invert: bool = if element_idx > 0 { false } else { true };
        let mut new_element = genome.seq[element_idx.abs() as usize - 1].clone(); // convert back to 0 indexed

        if invert {
            new_element.strand = !new_element.strand;
        }

        // update homology map for new element
        let element_id = new_element.element_id;
        let mut homology_group = &mut homology_map[element_id][genome.genome_id];
        homology_group.push(element_idx.abs() as usize - 1); // convert back to 0 indexed
        
        new_genome.push(new_element);
    }

    genome.seq = new_genome;

}

pub fn mutate_inter_genome (population: &mut Population) {
    // TODO add recombination between genomes, need to think of way of identifying homologous regions and if each recombination event 
    // has single or double crossovers
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mutation::{Distribution as MutationDistribution, MutationMap};
    use crate::population::{Genome, NucElement};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn make_test_genome() -> Genome {
        let base_map = StructureMutationMap {
            duplication_rate: 0.0,
            deletion_rate: 0.0,
            inversion_rate: 0.0,
            recombination_rate: 0.0,
            max_duplications: None,
        };
        let mut rng = StdRng::seed_from_u64(42);
        let sel_dist = MutationDistribution::new_uniform(0.0, 1.0).unwrap();

        Genome {
            identifier: "test".to_string(),
            genome_id: 0,
            parent: "root".to_string(),
            seq: vec![
                NucElement {
                    seqname: "chr1".to_string(),
                    element_id: 0,
                    feature_id: 0,
                    feature_type: "exon".to_string(),
                    selection_coefficient: 0.0,
                    seq: vec![],
                    mutation_map: MutationMap::new(0, 0, &vec![], &sel_dist, &mut rng),
                    strand: true,
                    structure_mutation_map: base_map.clone(),
                },
                NucElement {
                    seqname: "chr1".to_string(),
                    element_id: 1,
                    feature_id: 1,
                    feature_type: "exon".to_string(),
                    selection_coefficient: 0.0,
                    seq: vec![],
                    mutation_map: MutationMap::new(0, 0, &vec![], &sel_dist, &mut rng),
                    strand: false,
                    structure_mutation_map: base_map,
                },
            ],
        }
    }

    #[test]
    fn inversion_flips_strand() {
        let mut genome = make_test_genome();
        let before_strands: Vec<bool> = genome.seq.iter().map(|e| e.strand).collect();

        for e in &mut genome.seq {
            e.structure_mutation_map.inversion_rate = 1.0;
        }

        let mu = MutationDistribution::new_uniform(0.0, 1.0).unwrap();
        let pos = MutationDistribution::new_uniform(0.0, 1.0).unwrap();
        mutate_intra_genome(&mut genome, &mu, &pos);

        let after_strands: Vec<bool> = genome.seq.iter().map(|e| e.strand).collect();
        assert_ne!(before_strands, after_strands, "strands should flip after forced inversion");
        assert_eq!(genome.seq.len(), before_strands.len(), "inversion should preserve genome length");
    }

    #[test]
    fn deletion_reduces_genome_length() {
        let mut genome = make_test_genome();
        let before_len = genome.seq.len();

        for e in &mut genome.seq {
            e.structure_mutation_map.deletion_rate = 1.0;
        }

        let mu = MutationDistribution::new_uniform(0.0, 1.0).unwrap();
        let pos = MutationDistribution::new_uniform(0.0, 1.0).unwrap();
        mutate_intra_genome(&mut genome, &mu, &pos);

        assert!(genome.seq.len() < before_len, "genome should shrink after forced deletion");
    }

    #[test]
    fn translocation_preserves_length() {
        // A translocation is a simultaneous duplication (copy to new position) and
        // deletion (remove from original position).  The gene moves but the genome
        // length stays the same.
        let mut genome = make_test_genome();
        let before_len = genome.seq.len();
        let before_ids: Vec<usize> = genome.seq.iter().map(|e| e.feature_id).collect();

        for e in &mut genome.seq {
            // Exactly one duplication per element, then delete the original.
            e.structure_mutation_map.duplication_rate = 1.0;
            e.structure_mutation_map.max_duplications = Some(1);
            e.structure_mutation_map.deletion_rate = 1.0;
            e.structure_mutation_map.inversion_rate = 0.0;
        }

        // Use a non-zero offset so duplicates land somewhere other than position 0.
        let mu = MutationDistribution::new_uniform(0.0, 1.0).unwrap();
        let pos = MutationDistribution::new_uniform(1.0, 2.0).unwrap();
        mutate_intra_genome(&mut genome, &mu, &pos);

        assert_eq!(
            genome.seq.len(),
            before_len,
            "translocation should preserve genome length"
        );
        // The set of feature_ids present should be unchanged even if order differs.
        let mut after_ids: Vec<usize> = genome.seq.iter().map(|e| e.feature_id).collect();
        let mut expected = before_ids.clone();
        after_ids.sort();
        expected.sort();
        assert_eq!(after_ids, expected, "translocation should not add or remove genes");
    }
}