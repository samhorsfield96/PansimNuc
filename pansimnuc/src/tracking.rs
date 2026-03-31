use crate::demography::MetaPopulation;
use crate::population::Population;

// function to take information about NucElements and write to output file for tracking purposes
fn write_tracking_header(out_path: &str) {
    let mut wtr = csv::Writer::from_path(out_path).expect("Could not create output file for tracking.");
    
    // write header
    wtr.write_record(&["element_id", "feature_id", "population_id", "genome_id", "contig_id", "feature_type", "multiplier", "start", "end", "strand", "original_length", "length", "sequence", "log_selection_coefficients"]).expect("Could not write header to tracking output file.");
    
    wtr.flush().expect("Could not flush tracking output file.");
}

fn write_tracking_output(out_path: &str, metapopulation: &MetaPopulation) {
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