// in this script, genes can move around, be duplicated and deleted

use crate::{mutation::MutationMap, population::NucElement};
use std::collections::HashMap;
use crate::population::{Genome, Population};
use rand::seq::SliceRandom;
use rand::Rng;
use crate::mutation::Distribution as MutationDistribution;
extern crate levenshtein;
use levenshtein::levenshtein;

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
    pub duplication_insertion_prob: f64, // probability of duplication inserting upstream of the original position, as opposed to downstream
}

/// Returns a similarity score in [0.0, 1.0] based on normalized edit distance.
/// 1.0 = identical, 0.0 = completely different.
fn calculate_homology(a: &NucElement, b: &NucElement) -> f64 {
    let s = String::from_utf8_lossy(&a.seq);
    let t = String::from_utf8_lossy(&b.seq);

    let m = s.len();
    let n = t.len();

    if m == 0 && n == 0 { return 1.0; }
    if m == 0 || n == 0 { return 0.0; }

    let edit_distance = levenshtein(&s, &t);
    let max_len = m.max(n);
    1.0 - (edit_distance as f64 / max_len as f64)
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
            let pos_order: i64 = if pos_rand_val < element.structure_mutation_map.duplication_insertion_prob { -1 } else { 1 };
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

    let mut prev_contig_id = 0;

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

        // TODO need to determine contig breaks, likely using prev_contig_id, and reassign when new contig is reached
    }

    genome.seq = new_genome;

}

pub fn mutate_inter_genome (population: &mut Population) {
    // TODO use  homology map to identify homologous regions, sample
    // from poisson to determine how many recombination events occur, then size of track

    let mut thread_rng = rand::thread_rng();
    
    // get number of recombination events across whole population
    let n_recombinations = population.recombination_dists[0].sample(&mut thread_rng) as usize;

    let pop_size = population.pop.len();

    // All ordered pairs where donor != recipient
    let all_pairs: Vec<(usize, usize)> = (0..pop_size)
        .flat_map(|i| (0..pop_size).filter(move |&j| j != i).map(move |j| (i, j)))
        .collect();

    let mut recombination_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for _ in 0..n_recombinations {
        let (donor, recipient) = all_pairs.choose(&mut thread_rng).expect("Failed to select a random pair for recombination");
        recombination_map.entry(*donor).or_default().push(*recipient);
    }

    // iterate through each donor and recipient pair
    for (donor, recipients) in recombination_map {
        let donor_genome = &population.pop[donor];
        for recipient in recipients {
            let recipient_genome = &mut population.pop[recipient];

            // look for donor and recipient site, maximum 25 attempts, if not found, skip recombination event
            let mut donor_site_chosen : bool = false;
            let mut attempts = 0;
            let max_attempts = 25;

            let mut start_donor_site: usize;
            let mut start_recipient_site: usize;

            while !donor_site_chosen && attempts < max_attempts {
                let recombination_pos = thread_rng.gen_range(0..donor_genome.seq.len());

                // determine if position in both donor and recipient genome, if not, resample
                let recomb_element = population.homology_map[recombination_pos];

                let donor_has_site = !recomb_element[donor].is_empty();
                let recipient_has_site = recomb_element.len() > recipient && !recomb_element[recipient].is_empty();

                // if both vectors are not empty, then search through each and test to make sure they have sufficient homology
                if donor_has_site && recipient_has_site {
                    for donor_site in &recomb_element[donor] {
                        for recipient_site in &recomb_element[recipient] {
                            // check homology between donor and recipient site, if sufficient, break loop and move to recombination, if not, continue searching
                            let homology = calculate_homology(&donor_genome.seq[*donor_site], &recipient_genome.seq[*recipient_site]);
                            if homology >= population.recombination_threshold {
                                // perform recombination event, break out of loops

                                start_donor_site = *donor_site;
                                start_recipient_site = *recipient_site;

                                donor_site_chosen = true;
                                break;
                            }
                        }
                        if donor_site_chosen {
                            break;
                        }
                    }
                }  
                attempts += 1;
            }

            if donor_site_chosen {
                // now sample from poisson distribution to determine minumum size of recombination track
                let min_recombination_len = population.recombination_dists[1].sample(&mut thread_rng) as usize;

                // determine whether there is a track that can be recombined
                let mut track_found = false;
                let mut end_donor_site = start_donor_site;
                let mut end_recipient_site = start_recipient_site;
                let mut recombination_len = &donor_genome.seq[start_donor_site].seq.len();

                // ensure recombination occurs in single chromosome each
                let donor_contig_id = donor_genome.seq[start_donor_site].contig_id;
                let recipient_contig_id = recipient_genome.seq[start_recipient_site].contig_id;

                while !track_found {
                    // determine length of donor DNA
                    while recombination_len < &min_recombination_len {
                        let new_end_donor_site = end_donor_site + 1;
                        
                        // run off end of contig, assume complete recombination
                        if donor_genome.seq[new_end_donor_site].contig_id != donor_contig_id || new_end_donor_site >= donor_genome.seq.len() {
                            recombination_len += &donor_genome.seq[end_donor_site].seq.len();

                            // find end of recipient track
                            let mut contig_end = false;
                            while !contig_end {
                                let new_end_recipient_site = end_recipient_site + 1;
                                if recipient_genome.seq[new_end_recipient_site].contig_id != recipient_contig_id || new_end_recipient_site >= recipient_genome.seq.len() {
                                    contig_end = true;
                                } else {
                                    end_recipient_site = new_end_recipient_site;
                                }
                            }

                            track_found = true;
                            break;
                        }

                        // else continue going through contig
                        end_donor_site = new_end_donor_site;
                        recombination_len += &donor_genome.seq[end_donor_site].seq.len();


                    }

                    // use to determine homology between sites
                    let donor_site = &donor_genome.seq[end_donor_site];

                    // now iterate through recipient genome until homology found between end and donor site
                    let mut recipient_end_found = false;
                    while !recipient_end_found {
                        let new_end_recipient_site = end_recipient_site + 1;

                        // run off end of contig, assume complete recombination
                        if donor_genome.seq[new_end_donor_site].contig_id != recipient_contig_id || new_end_donor_site >= recipient_genome.seq.len() {
                            recipient_end_found = true;
                            break;
                        }

                        let recipient_site = &recipient_genome.seq[new_end_recipient_site];
                        let homology = calculate_homology(donor_site, recipient_site);
                        if homology >= population.recombination_threshold {
                            track_found = true;
                            recipient_end_found = true;
                        }

                        // continue iterating through recipient contig until homology found, or end of contig reached
                        end_recipient_site = new_end_recipient_site;
                    }
                }

                // perform recombination event, replacing recipient track with donor track
                // clone the donor track first, before any mutable borrow of population.pop
                let mut donor_track: Vec<NucElement> = population.pop[donor].seq[start_donor_site..=end_donor_site].to_vec();

                // update information from recipient track
                for element in &mut donor_track {
                    element.contig_id = recipient_contig_id;
                }

                // update homology map for recipient genome, need to add new positions for each element in donor track, and remove old positions for each element in recipient track
                // remove old positions
                for element_idx in start_recipient_site..=end_recipient_site {
                    let element_id = recipient_genome.seq[element_idx].element_id;
                    let homology_group = &mut population.homology_map[element_id][recipient_genome.genome_id];
                    homology_group.retain(|&pos| pos != element_idx); // remove old position
                }
                
                // now safe to mutably borrow recipient
                recipient_genome.seq.splice(start_recipient_site..=end_recipient_site, donor_track);

                // add new positions
                for element_idx in start_recipient_site..start_recipient_site + donor_track.len() {
                    let element_id = recipient_genome.seq[element_idx].element_id;
                    let homology_group = &mut population.homology_map[element_id][recipient_genome.genome_id];
                    homology_group.push(element_idx); // add new position
                }
            }
        }
    }
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
            duplication_insertion_prob: 0.5,
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
                    structure_mutation_map: base_map.clone(),
                },
                NucElement {
                    seqname: "chr1".to_string(),
                    element_id: 2,
                    feature_id: 2,
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
        
        let mut homology_map: Vec<Vec<Vec<usize>>> = Vec::new();
        for _ in genome.seq.iter() {
            homology_map.push(vec![vec![0]]);
        }

        mutate_intra_genome(&mut genome, &mut homology_map, &mu, &pos);

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

        let mut homology_map: Vec<Vec<Vec<usize>>> = Vec::new();
        for _ in genome.seq.iter() {
            homology_map.push(vec![vec![0]]);
        }

        let mu = MutationDistribution::new_uniform(0.0, 1.0).unwrap();
        let pos = MutationDistribution::new_uniform(0.0, 1.0).unwrap();
        let mut homology_map: Vec<Vec<Vec<usize>>> = vec![vec![vec![]], vec![vec![]]];
        mutate_intra_genome(&mut genome, &mut homology_map, &mu, &pos);

        assert!(genome.seq.len() < before_len, "genome should shrink after forced deletion");
    }

    #[test]
    // run multiple times, as may fail stochastically if duplication and deletion events don't line up as expected
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
                
        let mut homology_map: Vec<Vec<Vec<usize>>> = Vec::new();
        for _ in genome.seq.iter() {
            homology_map.push(vec![vec![0]]);
        }

        // Use a non-zero offset so duplicates land somewhere other than position 0.
        let mu = MutationDistribution::new_uniform(0.0, 0.9).unwrap();
        let pos = MutationDistribution::new_uniform(1.0, 2.0).unwrap();
        mutate_intra_genome(&mut genome, &mut homology_map, &mu, &pos);

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