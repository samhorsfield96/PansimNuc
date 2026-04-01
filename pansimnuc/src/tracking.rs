use crate::demography::MetaPopulation;
use crate::population::Population;

// function to take information about NucElements and write to output file for tracking purposes
pub fn write_tracking_header(out_path: &str) {
    let mut wtr = csv::Writer::from_path(out_path).expect("Could not create output file for tracking.");
    
    // write header
    wtr.write_record(&["element_id", "feature_id", "population_id", "genome_id", "contig_id", "feature_type", "multiplier", "start", "end", "strand", "original_length", "length", "sequence", "log_selection_coefficients"]).expect("Could not write header to tracking output file.");
    
    wtr.flush().expect("Could not flush tracking output file.");
}

pub fn identify_tracked_elements(population: &mut Population, tracking_regions: &Vec<(usize, usize, usize)>) {
    for genome in &mut population.pop {
        // determine position of each element in the genome, and write to output file
        let mut current_start = 0;
        let mut current_contig_id = 0;

        for element in &mut genome.seq {
            if element.contig_id != current_contig_id {
                // reset current start position for new contig
                current_start = 0;
                current_contig_id = element.contig_id;
            }

            let element_end = current_start + element.seq.len(); 

            for (contig_id, start, end) in tracking_regions {
                // make sure there is an overlap between the element and the tracking region, and that they are on the same contig
                if element.contig_id == *contig_id && current_start <= *end && *start <= element_end {
                    element.tracked = true;
                    break; // stop checking other tracking regions once a match is found
                }
            }

            current_start = element_end;
        }
    }
}

pub fn write_tracking_output(out_path: &str, metapopulation: &MetaPopulation) {
    let mut wtr = csv::Writer::from_path(out_path).expect("Could not create output file for tracking.");
    
    // write information for each NucElement in the population
    for population in &metapopulation.populations {
        for genome in &population.pop {
            // determine position of each element in the genome, and write to output file
            let mut current_start = 0;
            let mut current_contig_id = 0;

            for element in &genome.seq {
                if element.contig_id != current_contig_id {
                    // reset current start position for new contig
                    current_start = 0;
                    current_contig_id = element.contig_id;
                }

                let element_end = current_start + element.seq.len();

                if element.tracked {
                    let element_seq = element.seq.iter().map(|&base| Population::decode_base(base)).collect::<Vec<u8>>();
                    let element_selection_coefficients = element.generate_selection_coefficients(); // placeholder function to generate selection coefficients, can be implemented based on specific requirements

                    wtr.write_record(&[
                        element.element_id.to_string(),
                        element.feature_id.to_string(),
                        population.id.to_string(),
                        genome.genome_id.to_string(),
                        element.contig_id.to_string(),
                        element.feature_type.clone(),
                        element.multiplier.to_string(),
                        current_start.to_string(),
                        element_end.to_string(),
                        element.strand.to_string(),
                        element.original_length.to_string(),
                        (element_end - current_start).to_string(),
                        element_seq.iter().map(|&base| base as char).collect::<String>(),
                        element_selection_coefficients.iter().map(|&coeff| coeff.to_string()).collect::<Vec<String>>().join(";")
                    ]).expect("Could not write record to tracking output file.");
                }

                current_start = element_end;

            }
        }
    }
    
    wtr.flush().expect("Could not flush tracking output file.");
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mutation::{Distribution, MutationMap};
    use crate::population::{Genome, NucElement, Population};
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use std::collections::HashMap;

    fn make_element(contig_id: usize, seq_len: usize) -> NucElement {
        let seq = vec![1u8; seq_len];
        let dist = Distribution::new_uniform(0.0, 1.0).unwrap();
        let mut rng: StdRng = StdRng::seed_from_u64(0);
        NucElement {
            contig_id,
            element_id: 0,
            feature_id: 0,
            feature_type: "intergenic".to_string(),
            multiplier: 1.0,
            mutation_map: MutationMap::new(0, 0, &seq, &dist, &mut rng),
            seq,
            strand: true,
            original_length: seq_len,
            frameshift: false,
            tracked: false,
        }
    }

    fn make_population(elements: Vec<NucElement>) -> Population {
        let genome = Genome {
            identifier: "0".to_string(),
            genome_id: 0,
            contig_starts: vec![0],
            parent: "root".to_string(),
            seq_length: elements.iter().map(|e| e.seq.len()).sum(),
            seq: elements,
        };
        Population {
            id: 0,
            generation: 0,
            pop: vec![genome],
            core_vec: vec![],
            selection_dists: vec![],
            mu_dists: vec![],
            indel_dists: vec![],
            structural_mu_dists: vec![],
            recombination_dists: vec![],
            recombination_threshold: 0.9,
            homology_map: vec![],
            feature_map: HashMap::new(),
            max_multiplier_dist: 0,
            n_generations: 1,
            verbose: false,
        }
    }

    #[test]
    fn test_element_inside_region_is_tracked() {
        // two elements at [0,100) and [100,200); region covers [0,500)
        let elements = vec![make_element(0, 100), make_element(0, 100)];
        let mut pop = make_population(elements);

        identify_tracked_elements(&mut pop, &vec![(0, 0, 500)]);

        assert!(pop.pop[0].seq[0].tracked);
        assert!(pop.pop[0].seq[1].tracked);
    }

    #[test]
    fn test_element_outside_region_is_not_tracked() {
        // two elements at [0,100) and [100,200); region covers only [300,500)
        let elements = vec![make_element(0, 100), make_element(0, 100)];
        let mut pop = make_population(elements);

        identify_tracked_elements(&mut pop, &vec![(0, 300, 500)]);

        assert!(!pop.pop[0].seq[0].tracked);
        assert!(!pop.pop[0].seq[1].tracked);
    }

    #[test]
    fn test_wrong_contig_is_not_tracked() {
        // element on contig 0; region on contig 1
        let elements = vec![make_element(0, 100)];
        let mut pop = make_population(elements);

        identify_tracked_elements(&mut pop, &vec![(1, 0, 500)]);

        assert!(!pop.pop[0].seq[0].tracked);
    }

    #[test]
    fn test_element_partially_overlapping_region_is_tracked() {
        // element[0] at [0,150), element[1] at [150,250); region [200,500)
        // element[1] end (250) is inside region, so it overlaps
        let elements = vec![make_element(0, 150), make_element(0, 100)];
        let mut pop = make_population(elements);

        identify_tracked_elements(&mut pop, &vec![(0, 200, 500)]);

        assert!(!pop.pop[0].seq[0].tracked); // [0,150) — no overlap with [200,500)
        assert!(pop.pop[0].seq[1].tracked);  // [150,250) — overlaps [200,500)
    }

    #[test]
    fn test_element_spanning_entire_region_is_tracked() {
        // element at [0,1000); region [200,500) — element completely contains region
        let elements = vec![make_element(0, 1000)];
        let mut pop = make_population(elements);

        identify_tracked_elements(&mut pop, &vec![(0, 200, 500)]);

        assert!(pop.pop[0].seq[0].tracked);
    }
}