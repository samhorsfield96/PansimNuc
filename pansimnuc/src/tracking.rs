use crate::demography::MetaPopulation;
use crate::population::Population;
use std::fs::OpenOptions;
use std::path::Path;
use crate::population::NucElement;

// function to take information about NucElements and write to output file for tracking purposes
pub fn write_tracking_header(out_path: &str) {
    let mut wtr = csv::Writer::from_path(out_path).expect("Could not create output file for tracking.");
    
    // write header
    wtr.write_record(&["element_id", "feature_id", "generation", "population_id", "genome_id", "contig_id", "feature_type", "multiplier", "start", "end", "strand", "original_length", "length", "sequence", "log_selection_coefficients"]).expect("Could not write header to tracking output file.");
    
    wtr.flush().expect("Could not flush tracking output file.");
}

pub fn identify_tracked_element(element: &mut NucElement, element_start: usize, tracking_regions: &Vec<(String, usize, usize)>, contig_name_to_id: &Vec<String>)
{
    let element_end = element_start + element.seq.len(); 
    for (contig_id, start, end) in tracking_regions {
        // make sure there is an overlap between the element and the tracking region, and that they are on the same contig
        let element_contig_name = &contig_name_to_id[element.contig_id];
        
        if element_contig_name == contig_id && element_start <= *end && *start <= element_end {
            element.tracked = true;
            break; // stop checking other tracking regions once a match is found
        }
    }
}

pub fn identify_tracked_elements(population: &mut Population, tracking_regions: &Vec<(String, usize, usize)>, contig_name_to_id: &Vec<String>) {
    
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
            identify_tracked_element(element, current_start, tracking_regions, contig_name_to_id);
            current_start += element.seq.len(); 
        }
    }
}

pub fn write_tracking_output(out_path: &str, metapopulation: &MetaPopulation) {
    // Check if file exists; error if it doesn't
    if !Path::new(out_path).exists() {
        panic!("Tracking output file does not exist: {}", out_path);
    }
    
    // Open file in append mode
    let file = OpenOptions::new()
        .append(true)
        .open(out_path)
        .expect("Could not open tracking output file for appending.");
    
    let mut wtr = csv::Writer::from_writer(file);
    
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
                    let element_selection_coefficients = element.generate_selection_coefficients();

                    wtr.write_record(&[
                        element.element_id.to_string(),
                        element.feature_id.to_string(),
                        population.generation.to_string(),
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
    use crate::demography::MetaPopulation;
    use crate::mutation::{Distribution, MutationMap};
    use crate::population::{Genome, NucElement, Population};
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use std::collections::HashMap;
    use std::sync::Arc;

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
            mutation_map: Arc::new(MutationMap::new(0, 0, &seq, &dist, &mut rng)),
            seq: Arc::new(seq),
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
            total_exon_length: 0,
            total_intron_length: 0,
            total_intergenic_length: 0,
            total_te_cut_length: 0,
            total_te_copy_length: 0,
            total_tracking_length: 0,
            total_elements: 0,
            total_exon_elements: 0,
            total_intron_elements: 0,
            total_intergenic_elements: 0,
            total_te_cut_elements: 0,
            total_te_copy_elements: 0,
            total_tracking_elements: 0,
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
            augment_tracking: false,
            genome_size_penalty_per_bp: 0.0,
            optimal_genome_size: 0,
        }
    }

    #[test]
    fn test_element_inside_region_is_tracked() {
        // two elements at [0,100) and [100,200); region covers [0,500)
        let elements = vec![make_element(0, 100), make_element(0, 100)];
        let mut pop = make_population(elements);

        identify_tracked_elements(&mut pop, &vec![("0".to_string(), 0, 500)], &vec!["0".to_string()]);

        assert!(pop.pop[0].seq[0].tracked);
        assert!(pop.pop[0].seq[1].tracked);
    }

    #[test]
    fn test_element_outside_region_is_not_tracked() {
        // two elements at [0,100) and [100,200); region covers only [300,500)
        let elements = vec![make_element(0, 100), make_element(0, 100)];
        let mut pop = make_population(elements);

        identify_tracked_elements(&mut pop, &vec![("0".to_string(), 300, 500)], &vec!["0".to_string()]);

        assert!(!pop.pop[0].seq[0].tracked);
        assert!(!pop.pop[0].seq[1].tracked);
    }

    #[test]
    fn test_wrong_contig_is_not_tracked() {
        // element on contig 0; region on contig 1
        let elements = vec![make_element(0, 100)];
        let mut pop = make_population(elements);

        identify_tracked_elements(&mut pop, &vec![("1".to_string(), 0, 500)], &vec!["0".to_string()]);

        assert!(!pop.pop[0].seq[0].tracked);
    }

    #[test]
    fn test_element_partially_overlapping_region_is_tracked() {
        // element[0] at [0,150), element[1] at [150,250); region [200,500)
        // element[1] end (250) is inside region, so it overlaps
        let elements = vec![make_element(0, 150), make_element(0, 100)];
        let mut pop = make_population(elements);

        identify_tracked_elements(&mut pop, &vec![("0".to_string(), 200, 500)], &vec!["0".to_string()]);

        assert!(!pop.pop[0].seq[0].tracked); // [0,150) — no overlap with [200,500)
        assert!(pop.pop[0].seq[1].tracked);  // [150,250) — overlaps [200,500)
    }

    #[test]
    fn test_element_spanning_entire_region_is_tracked() {
        // element at [0,1000); region [200,500) — element completely contains region
        let elements = vec![make_element(0, 1000)];
        let mut pop = make_population(elements);

        identify_tracked_elements(&mut pop, &vec![("0".to_string(), 200, 500)], &vec!["0".to_string()]);

        assert!(pop.pop[0].seq[0].tracked);
    }

    #[test]
    fn test_write_tracking_output_appends_to_existing_file() {
        use std::fs;
        use std::io::Read;
        use crate::demography::MetaPopulation;

        // Create a temporary file path
        let temp_path = std::env::temp_dir().join(format!(
            "pansimnuc_tracking_test_{}.csv",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let temp_path_str = temp_path.to_string_lossy().into_owned();

        // Write header to the file
        write_tracking_header(&temp_path_str);

        // Create a population with a tracked element
        let mut elements = vec![make_element(0, 100)];
        elements[0].tracked = true;
        let mut pop = make_population(elements);
        pop.id = 0;
        pop.generation = 1;

        // Create a metapopulation
        let metapop = MetaPopulation {
            populations: vec![pop],
            population_split_config: crate::config::PopulationSplitConfig {
                population_splits: vec![1],
                generation_splits: vec![1],
                migration_rate: 0.0,
            },
            n_generations: 1,
            recombination_rate: 0.0,
            recombination_size_mean: 1.0,
            site_mutation_mus_vals: vec![],
            site_indel_mus_vals: vec![],
        };

        // Write tracking output (should append)
        write_tracking_output(&temp_path_str, &metapop);

        // Read the file and verify contents
        let mut content = String::new();
        fs::File::open(&temp_path_str)
            .expect("Could not open tracking file")
            .read_to_string(&mut content)
            .expect("Could not read tracking file");

        // Should have two lines: header and one data row
        let lines: Vec<&str> = content.lines().collect();
        assert!(lines.len() >= 2, "Expected at least 2 lines (header + data)");
        
        // Verify header line
        assert!(lines[0].contains("element_id"));
        assert!(lines[0].contains("generation"));
        
        // Verify data line contains expected values
        assert!(lines[1].contains("0")); // element_id
        assert!(lines[1].contains("1")); // generation (pop.generation = 1)

        // Clean up
        let _ = fs::remove_file(&temp_path_str);
    }

    #[test]
    #[should_panic(expected = "does not exist")]
    fn test_write_tracking_output_panics_if_file_not_exists() {
        use crate::demography::MetaPopulation;

        let nonexistent_path = "/tmp/this_file_should_not_exist_pansimnuc_12345.csv";
        let elements = vec![make_element(0, 100)];
        let pop = make_population(elements);
        let metapop = MetaPopulation {
            populations: vec![pop],
            population_split_config: crate::config::PopulationSplitConfig {
                population_splits: vec![1],
                generation_splits: vec![1],
                migration_rate: 0.0,
            },
            n_generations: 1,
            recombination_rate: 0.0,
            recombination_size_mean: 1.0,
            site_mutation_mus_vals: vec![],
            site_indel_mus_vals: vec![],
        };

        // Should panic because file doesn't exist
        write_tracking_output(nonexistent_path, &metapop);
    }
}