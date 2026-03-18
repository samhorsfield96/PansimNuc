use crate::gff::FeaturePos;
use crate::mutation::MutationMap;
use crate::mutation::Distribution as MutationDistribution;
use crate::structural::mutate_intra_genome;
use crate::structural::mutate_inter_genome;
use crate::structural::StructureMutationMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use rayon::prelude::*;
use logsumexp::LogSumExp;
use rand::distributions::{Distribution as RandDistribution, WeightedIndex};
use rand::rngs::StdRng;
use rand::{SeedableRng};

#[derive(Clone)]
pub struct NucElement {
    pub contig_id: usize,
    pub element_id: usize,
    pub feature_id: usize,
    pub feature_type: String,
    pub multiplier: f64,
    pub seq: Vec<u8>,
    pub mutation_map: MutationMap,
    pub strand: bool,
    pub structure_mutation_map: StructureMutationMap,
}

impl NucElement {
    fn element_selection_coefficient(&self, genome_identifier: &str) -> f64 {
        let mut element_log_sum = 0.0;
        for (site, allele) in self.seq.iter().enumerate() {
            let allele_shifted = 1 >> allele;
            if let Some(coeff) = self.mutation_map.get(*allele, site) {
                element_log_sum += coeff;
            } else {
                panic!(
                    "Failed to generate selection coefficient for genome {} allele {} (shifted {}) at site {}",
                    genome_identifier, allele, allele_shifted, site
                );
            }
        }
        element_log_sum
    }
}

pub struct Genome {
    pub identifier: String,
    pub genome_id: usize,
    pub contig_starts: Vec<usize>,
    pub parent: String,
    pub seq: Vec<NucElement>,
}

impl Genome {
    pub fn update_contig_starts(&mut self) {
        self.contig_starts.clear();
        let mut current_start = 0;
        for (idx, element) in self.seq.iter().enumerate() {
            if idx == 0 || element.contig_id != self.seq[idx - 1].contig_id {
                self.contig_starts.push(current_start);
                current_start = idx;
            }
        }
    }
}

pub struct Population{
    pub generation: usize,
    pub pop: Vec<Genome>,
    pub core_vec: Vec<Vec<u8>>,
    pub selection_dists: Vec<MutationDistribution>,
    pub mu_dists: Vec<MutationDistribution>,
    pub recombination_dists: Vec<MutationDistribution>,
    pub recombination_threshold: f64,
    pub homology_map: Vec<Vec<Vec<usize>>>, // Map from original element ID to positions of homologous regions in other genomes, outermost loop is the homology group, middle loop is genomes, inner loop is positions
    pub feature_map: HashMap<usize, Vec<usize>>, // Map from feature ID to genes that share same ID
}

impl Population {
    fn is_te_feature(feature_type: &str) -> bool {
        feature_type == "TE-CUT" || feature_type == "TE-COPY"
    }

    fn check_feature_order (&self, genome: &Genome, element_idx: usize, element: &NucElement) -> (bool, f64) {
        let mut feature_broken = false;
        let mut feature_multiplier = 1.0;
        
        // identify regions where order matters
        if element.feature_type == "exon" || element.feature_type == "intron" {
            let feature_map_entry = self.feature_map.get(&element.feature_id).expect("Entry missing from feature_map");
            let feature_map_entry_len = feature_map_entry.len();

            // get position of element in feature_map_entry
            let position = feature_map_entry.iter().position(|n| n == &element.element_id).expect("Element not found in feature_map_entry");
            
            // check upstream elements in feature_map_entry
            for feature_element_idx in 0..position {
                // get upstream feature in genome
                let actual_element = &genome.seq[element_idx - (position - feature_element_idx)];
                let actual_element_id = actual_element.element_id;

                // get expected element
                let expected_element_id = feature_map_entry[feature_element_idx];

                // TODO - check whether indexing makes sense with reversal

                // expected element and strand matches so continue
                if actual_element_id == expected_element_id && actual_element.strand == element.strand {
                    continue;
                } else {
                    // check if reversed order and strand matches, if so continue, as likely to be functional just reversed
                    // check if last feature matches in feature_map_entry
                    let reversed_feature_id = feature_map_entry[(feature_map_entry_len - 1) - feature_element_idx];
                    if actual_element_id == reversed_feature_id && actual_element.strand == element.strand {
                        continue;
                    } else {
                        // if not, set multiplier to 0, as likely to be non-functional
                        feature_broken = true;
                        break;
                    }
                }
            }

            // check downstream elements in feature_map_entry
            if !feature_broken {
                for feature_element_idx in position + 1..feature_map_entry_len {
                    // get downstream feature in genome
                    let actual_element = &genome.seq[element_idx + (feature_element_idx - position)];
                    let actual_element_id = actual_element.element_id;

                    // get expected element
                    let expected_element_id = feature_map_entry[feature_element_idx];

                    // expected element and strand matches so continue
                    if actual_element_id == expected_element_id && actual_element.strand == element.strand {
                        continue;
                    } else {
                        // check if reversed order and strand matches, if so continue, as likely to be functional just reversed
                        // check if last feature matches in feature_map_entry
                        let reversed_feature_id = feature_map_entry[(feature_map_entry_len - 1) - feature_element_idx];
                        if actual_element_id == reversed_feature_id && actual_element.strand == element.strand {
                            continue;
                        } else {
                            // if not, set multiplier to 0, as likely to be non-functional
                            feature_broken = true;
                            break;
                        }
                    }  
                }
            }

            // check if TE is present upstream or downstream, if so increase multiplier, as likely to increase expression of gene, and thus fitness contribution
            if !feature_broken {
                // check not at beginning of genome
                if element_idx >= position + 1 {
                    let upstream_element = &genome.seq[element_idx - (position + 1)];
                    if Self::is_te_feature(&upstream_element.feature_type) {
                        feature_multiplier = upstream_element.multiplier;
                    }
                }

                // check not at end of genome
                if element_idx + (feature_map_entry_len - position) < genome.seq.len() {
                    let downstream_element = &genome.seq[element_idx + (feature_map_entry_len - position)];
                    
                    if Self::is_te_feature(&downstream_element.feature_type) {
                        if feature_multiplier.abs() < downstream_element.multiplier.abs() {
                            feature_multiplier = downstream_element.multiplier;
                        }
                    }
                }
            }
        }
        
        (feature_broken, feature_multiplier)
    }

    fn genome_output_path(output_path: &str, genome_index: usize) -> io::Result<PathBuf> {
        let path = Path::new(output_path);
        let file_name = path.file_name().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Output path must include a file name: {output_path}"),
            )
        })?;

        let prefixed_name = format!("{genome_index}_{}", file_name.to_string_lossy());
        Ok(path.with_file_name(prefixed_name))
    }

    fn decode_base(base: u8) -> u8 {
        match base {
            1 => b'A',
            2 => b'C',
            4 => b'G',
            8 => b'T',
            16 => b'N',
            _ => panic!("Invalid base encoding: {}", base),
        }
    }

    fn decode_sequence(seq: &[u8]) -> String {
        seq.iter().map(|&base| Self::decode_base(base) as char).collect()
    }

    fn genome_selection_coefficient(&self, genome: &Genome) -> f64 {
        let mut log_sum = 0.0;

        for (element_idx, element) in genome.seq.iter().enumerate() {
            let mut element_log_sum = 0.0;

            let (feature_broken, feature_multiplier) = self.check_feature_order(genome, element_idx, element);
            
            // if feature not broken, add up sites
            if !feature_broken {
                element_log_sum += element.element_selection_coefficient(&genome.identifier);
            }

            log_sum += element_log_sum * feature_multiplier;
        }

        log_sum
    }

    pub fn new(
        root: Vec<Vec<FeaturePos>>,
        n_genomes: usize,
        selection_dists: Vec<MutationDistribution>,
        mu_dists: Vec<MutationDistribution>,
        recombination_dists: Vec<MutationDistribution>,
        recombination_threshold: f64,
        structural_dists: Vec<StructureMutationMap>,
        rng: &mut StdRng,
    ) -> Self {
        // initialise population
        let mut population: Vec<Genome> = Vec::new();
        let mut genome: Vec<NucElement> = Vec::new();

        // count for element ID, each NucElement gets own to signal it it's homology group
        let mut element_id: usize = 0;

        // initialise homology map, outermost loop is the homology group, middle loop is genomes, inner loop is positions
        let mut homology_map: Vec<Vec<Vec<usize>>> = Vec::new();

        // initialise feature map, maps feature ID to number of genes that should share same ID
        let mut feature_map: HashMap<usize, Vec<usize>> = HashMap::new();

        // generate starting genome
        for (contig_id, features) in root.iter().enumerate() {
            for feature in features {
                
                // TODO change this so that can specify different mutation rates per site
                // also add separate TE compartment
                let selection_dist_id:usize = match feature.feature_type.as_str() {
                    "exon" => 0,
                    "intron" => 1,
                    "intergenic" => 2,
                    "TE-CUT" => 3,
                    "TE-COPY" => 4,
                    _ => panic!("Unknown feature type: {}", feature.feature_type),
                };
                
                let mu_dist_id:usize = match feature.feature_type.as_str() {
                    "exon" => 0,
                    "intron" => 1,
                    "intergenic" => 2,
                    "TE-CUT" => 3,
                    "TE-COPY" => 4,
                    _ => panic!("Unknown feature type: {}", feature.feature_type),
                };

                let structural_map:StructureMutationMap = match feature.feature_type.as_str() {
                    "exon" => structural_dists[0].clone(),
                    "intron" => structural_dists[1].clone(),
                    "intergenic" => structural_dists[2].clone(),
                    "TE-CUT" => structural_dists[3].clone(),
                    "TE-COPY" => structural_dists[4].clone(),
                    _ => panic!("Unknown feature type: {}", feature.feature_type),
                };

                // Update feature_map, keeping track of how many features there are
                if feature.feature_id != 0 {
                    feature_map.entry(feature.feature_id).or_default().push(element_id);
                }

                genome.push(NucElement {
                    contig_id: contig_id,
                    element_id: element_id,
                    feature_id: feature.feature_id,
                    feature_type: feature.feature_type.clone(),
                    seq: feature.seq.clone(),
                    strand: feature.strand,
                    mutation_map: MutationMap::new(selection_dist_id, mu_dist_id, &feature.seq, &selection_dists[selection_dist_id], rng),
                    structure_mutation_map: structural_map.clone(),
                    multiplier: 1.0, // Initialize with a default value, can be updated later
                });
                element_id += 1;

                // generate homology map for this element, initially just self
                let mut element_homology_map: Vec<Vec<usize>> = Vec::new();
                for _ in 0..n_genomes {
                    element_homology_map.push(vec![element_id]);
                }
                homology_map.push(element_homology_map);
            }
        }

        // copy whole genome to start
        for i in 0..n_genomes {
            let mut genome_entry = Genome {
                identifier: format!("{}", i),
                genome_id: i,
                contig_starts: Vec::new(), // will be updated after mutations
                parent: "root".to_string(),
                seq: genome.clone(),
            };
            genome_entry.update_contig_starts();

            population.push(genome_entry);
        }
        
        let core_vec: Vec<Vec<u8>> =
            vec![vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4]];

        Self {
            generation: 0,
            pop: population,
            core_vec,
            selection_dists,
            mu_dists,
            recombination_dists,
            recombination_threshold,
            homology_map,
            feature_map
        }
    }

    // mutate individuals in the population according to their mutation maps and the provided distributions
    pub fn mutate (&mut self) {
        for genome in &mut self.pop {
            for element in &mut genome.seq {
                element.mutation_map.mutate(&self.core_vec, 
                    &mut element.seq, 
                    &self.selection_dists[element.mutation_map.selection_dist_id], 
                    &self.mu_dists[element.mutation_map.mu_dist_id],);
            }
        }
    }   

    pub fn structural_intra_genome(&mut self) {
        // probabilities for duplications
        let duplication_mu_dist = MutationDistribution::new_uniform(0.0, 1.0). expect("Failed to create uniform distribution for duplications"); // TODO all setting
        let duplication_pos_dist = MutationDistribution::new_poisson(1.0).expect("Failed to create poisson distribution for duplications");

        // make probabilities very high to favour only the most transposable of elements
        let translocation_mu_dist = MutationDistribution::new_uniform(0.9, 1.0).expect("Failed to create uniform distribution for translocations"); // TODO all setting of this to favour translocations
        let translocation_pos_dist = MutationDistribution::new_poisson(1000.0).expect("Failed to create poisson distribution for translocations"); // TODO all setting of this to favour translocations

        // duplications
        self.pop
            .par_iter_mut()
            .for_each(|genome| {
                
                // duplications
                mutate_intra_genome(genome, &duplication_mu_dist, &duplication_pos_dist);

                // translocations
                // clear homology map to enable fresh creation of groups
                mutate_intra_genome(genome, &translocation_mu_dist, &translocation_pos_dist);
            });
        
        // update homology map for all new elements
        for genome in &self.pop {
            self.homology_map.iter_mut().for_each(|element_homology_map| {
                    element_homology_map[genome.genome_id].clear();
            });
            
            for (element_idx, element) in genome.seq.iter().enumerate() {
                let element_id = element.element_id;
                let homology_group = &mut self.homology_map[element_id][genome.genome_id];
                homology_group.push(element_idx); // convert back to 0 indexed
            }
        }
    }

    pub fn structural_inter_genome(&mut self) {
        // recombination
        mutate_inter_genome(self);
    }

    // sample individuals using logsumexp normalisation to prevent underflow/overflow issues with very small/large weights
    pub fn sample_individuals (&mut self, rng: &mut StdRng) -> Vec<usize> {
        let mut selection_weights: Vec<f64> = vec![1.0; self.pop.len()];
        
        selection_weights = self.pop
            .par_iter()
            .map(|genome| {
                let log_sum = self.genome_selection_coefficient(genome);
                log_sum
            })
            .collect();

        // logsumexp normalization to prevent underflow/overflow issues with very small/large weights
        let logsumexp_value = selection_weights.iter().ln_sum_exp();

        selection_weights = selection_weights.into_iter()
            .map(|x| (x - logsumexp_value).exp()) // exp(log(w) - logsumexp)
            .collect();

        let sum_weights: f64 = selection_weights.iter().sum();
        selection_weights = selection_weights.iter().map(|&w| if w != std::f64::NEG_INFINITY {w / sum_weights} else {0.0}).collect();

        // Create a WeightedIndex distribution based on weights
        let sampling_dist = WeightedIndex::new(&selection_weights).expect("Failed to generate sampling index for population.");

        // Sample rows based on the distribution
        let sampled_indices: Vec<usize> = (0..self.pop.len()).map(|_| sampling_dist.sample(rng)).collect();

        #[cfg(debug_assertions)] {
            eprintln!("Selection weights: {:?}", selection_weights);
            eprintln!("Sampled indices: {:?}", sampled_indices);
        }

        sampled_indices
    }

    pub fn next_generation (&mut self, sampled_indices: Vec<usize>) {
        let mut new_pop: Vec<Genome> = Vec::new();
        let mut new_homology_map: Vec<Vec<Vec<usize>>> = Vec::new();
        let mut genome_id = 0;

        // update genomes and homology map
        for &selected_index in &sampled_indices {
            let selected_genome = &self.pop[selected_index];
            new_pop.push(Genome {
                identifier: format!("{}-{}", self.generation + 1, selected_genome.identifier),
                genome_id: genome_id,
                contig_starts: selected_genome.contig_starts.clone(),
                parent: selected_genome.identifier.clone(),
                seq: selected_genome.seq.clone(),
            });
            genome_id += 1;
        }

        // update homology map
        for element_homology_map in &self.homology_map {
            let mut new_element_homology_map: Vec<Vec<usize>> = Vec::new();
            for &selected_index in &sampled_indices {
                new_element_homology_map.push(element_homology_map[selected_index].clone());
            }
            new_homology_map.push(new_element_homology_map);
        }

        self.pop = new_pop;
        self.homology_map = new_homology_map;
        self.generation += 1;
    }

    pub fn write_fasta(&self, output_path: &str) -> io::Result<()> {
        for (genome_index, genome) in self.pop.iter().enumerate() {
            let genome_output_path = Self::genome_output_path(output_path, genome_index)?;
            let file = File::create(&genome_output_path)?;
            let mut writer = BufWriter::new(file);

            // Group element indices by seqname
            let mut contig_groups: HashMap<usize, Vec<usize>> = HashMap::new();
            for (idx, element) in genome.seq.iter().enumerate() {
                contig_groups.entry(element.contig_id)
                    .or_insert_with(Vec::new)
                    .push(idx);
            }

            // Write each contig group as a separate FASTA entry
            for (contig_id, indices) in contig_groups {
                writeln!(
                    writer,
                    ">{id}_contig{contig_id} parent={parent} generation={generation}",
                    id = genome.identifier,
                    contig_id = contig_id,
                    parent = genome.parent,
                    generation = self.generation
                )?;

                let mut wrapped_line_len = 0usize;
                for idx in indices {
                    for &base in &genome.seq[idx].seq {
                        writer.write_all(&[Self::decode_base(base)])?;
                        wrapped_line_len += 1;

                        if wrapped_line_len == 80 {
                            writer.write_all(b"\n")?;
                            wrapped_line_len = 0;
                        }
                    }
                }

                if wrapped_line_len > 0 {
                    writer.write_all(b"\n")?;
                }
            }

            writer.flush()?;
        }

        Ok(())
    }

    pub fn write_gff(&self, output_path: &str) -> io::Result<()> {
        for (genome_index, genome) in self.pop.iter().enumerate() {
            let genome_output_path = Self::genome_output_path(output_path, genome_index)?;
            let file = File::create(&genome_output_path)?;
            let mut writer = BufWriter::new(file);

            writeln!(writer, "##gff-version 3")?;

            let genome_selection_coefficient = self.genome_selection_coefficient(genome);
            let mut contig_offsets: HashMap<usize, usize> = HashMap::new();

            for element in &genome.seq {
                let offset = contig_offsets.entry(element.contig_id).or_insert(0);
                let start_0 = *offset;
                let end_0 = start_0 + element.seq.len();
                *offset = end_0;

                if start_0 >= end_0 {
                    continue;
                }

                let element_selection_coefficient = element.element_selection_coefficient(&genome.identifier);

                let seq_id = format!("contig_{}", element.contig_id + 1);
                let start_1based = start_0 + 1;
                let end_1based = end_0;
                let strand = if element.strand { "+" } else { "-" };
                let max_duplications = element
                    .structure_mutation_map
                    .max_duplications
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "none".to_string());

                let attributes = format!(
                    "genome_id={};element_id={};feature_type={};feature_id={};contig_id={};genome_identifier={};parent={};multiplier={:.6};sequence_length={};genome_selection_coefficient={:.6};element_selection_coefficient={:.6};sv_duplication_rate={:.6};sv_deletion_rate={:.6};sv_inversion_rate={:.6};sv_max_duplications={};sv_duplication_insertion_prob={:.6}",
                    genome.genome_id,
                    element.element_id,
                    element.feature_type,
                    element.feature_id,
                    element.contig_id,
                    genome.identifier,
                    genome.parent,
                    element.multiplier,
                    element.seq.len(),
                    genome_selection_coefficient,
                    element_selection_coefficient,
                    element.structure_mutation_map.duplication_rate,
                    element.structure_mutation_map.deletion_rate,
                    element.structure_mutation_map.inversion_rate,
                    max_duplications,
                    element.structure_mutation_map.duplication_insertion_prob,
                );

                writeln!(
                    writer,
                    "{}\tPansimNuc\t{}\t{}\t{}\t.\t{}\t.\t{}",
                    seq_id,
                    element.feature_type,
                    start_1based,
                    end_1based,
                    strand,
                    attributes
                )?;
            }

            writer.flush()?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gff::FeaturePos;
    use std::fs;
    use std::io::Read;

    fn default_structural_dists() -> Vec<StructureMutationMap> {
        vec![
            StructureMutationMap {
                duplication_rate: 0.0,
                deletion_rate: 0.0,
                inversion_rate: 0.0,
                max_duplications: None,
                duplication_insertion_prob: 0.5,
            },
            StructureMutationMap {
                duplication_rate: 0.0,
                deletion_rate: 0.0,
                inversion_rate: 0.0,
                max_duplications: None,
                duplication_insertion_prob: 0.5,
            },
            StructureMutationMap {
                duplication_rate: 0.0,
                deletion_rate: 0.0,
                inversion_rate: 0.0,
                max_duplications: None,
                duplication_insertion_prob: 0.5,
            },
            StructureMutationMap {
                duplication_rate: 0.0,
                deletion_rate: 0.0,
                inversion_rate: 0.0,
                max_duplications: None,
                duplication_insertion_prob: 0.5,
            },
            StructureMutationMap {
                duplication_rate: 0.0,
                deletion_rate: 0.0,
                inversion_rate: 0.0,
                max_duplications: None,
                duplication_insertion_prob: 0.5,
            },
        ]
    }

    #[test]
    fn test_population_new() {
        // Create test data
        let mut root = Vec::new();
        let features = vec![
            FeaturePos {
                contig_id: 0,
                feature_id: 0,
                feature_type: "exon".to_string(),
                start: 100,
                end: 200,
                strand: true,
                seq: vec![1, 2, 4, 8], // ACGT
            },
            FeaturePos {
                contig_id: 0,
                feature_id: 1,
                feature_type: "intron".to_string(),
                start: 300,
                end: 400,
                strand: false,
                seq: vec![8, 4, 2, 1], // TGCA
            },
        ];
        root.push(features);

        let n_genomes = 3;
        let exon_dist = MutationDistribution::new_double_exp(0.5, 2.0, 0.3).expect("Failed to create double exponential distribution for exon features");
        let intron_dist = MutationDistribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intron features");
        let intergenic_dist = MutationDistribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intergenic features");
        let site_mutation_dists = vec![exon_dist, intron_dist, intergenic_dist];

        let exon_mu = MutationDistribution::new_uniform(0.0, 1.0).expect("Failed to create double exponential distribution for exon features");
        let intron_mu = MutationDistribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intron features");
        let intergenic_mu = MutationDistribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intergenic features");
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];

        let recombination_prob_dist = MutationDistribution::new_poisson(1.0).expect("Failed to create uniform distribution for recombination");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0).expect("Failed to create uniform distribution for recombination");
        
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let pop = Population::new(
            root,
            n_genomes,
            site_mutation_dists,
            site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            &mut rng,
        );

        // Check population was created correctly
        assert_eq!(pop.generation, 0);
        assert_eq!(pop.pop.len(), n_genomes);

        // Check each genome
        for (i, genome) in pop.pop.iter().enumerate() {
            assert_eq!(genome.identifier, format!("{}", i));
            assert_eq!(genome.parent, "root");
            assert_eq!(genome.seq.len(), 2); // Two features

            // Check first feature
            assert_eq!(genome.seq[0].contig_id, 0);
            assert_eq!(genome.seq[0].feature_id, 0);
            assert_eq!(genome.seq[0].feature_type, "exon");

            // Check second feature
            assert_eq!(genome.seq[1].contig_id, 0);
            assert_eq!(genome.seq[1].feature_id, 1);
            assert_eq!(genome.seq[1].feature_type, "intron");
        }
    }

    #[test]
    fn test_write_fasta_decodes_nucleotides() {
        let mut root = Vec::new();
        let features = vec![FeaturePos {
            contig_id: 0,
            feature_id: 0,
            feature_type: "exon".to_string(),
            start: 0,
            end: 4,
            strand: true,
            seq: vec![1, 2, 4, 8],
        }];
        root.push(features);

        let exon_dist = MutationDistribution::new_double_exp(0.5, 2.0, 0.3)
            .expect("Failed to create double exponential distribution for exon features");
        let intron_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create uniform distribution for intron features");
        let intergenic_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create uniform distribution for intergenic features");
        let site_mutation_dists = vec![exon_dist, intron_dist, intergenic_dist];

        let exon_mu = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create double exponential distribution for exon features");
        let intron_mu = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create uniform distribution for intron features");
        let intergenic_mu = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create uniform distribution for intergenic features");
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];
        let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let pop = Population::new(
            root,
            1,
            site_mutation_dists,
            site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            &mut rng,
        );

        let temp_path = std::env::temp_dir().join(format!(
            "pansimnuc_pop_{}.fasta",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("system time before UNIX_EPOCH")
                .as_nanos()
        ));
        let output_path = temp_path.to_string_lossy().into_owned();
        let genome_output_path = Population::genome_output_path(&output_path, 0)
            .expect("failed to construct per-genome output path");

        pop.write_fasta(&output_path)
            .expect("failed to write test FASTA file");

        let mut content = String::new();
        fs::File::open(&genome_output_path)
            .expect("failed to open test FASTA file")
            .read_to_string(&mut content)
            .expect("failed to read test FASTA file");

        assert!(content.contains(">0_contig0 parent=root generation=0"));
        assert!(content.contains("ACGT"));

        let _ = fs::remove_file(genome_output_path);
    }

    #[test]
    fn test_write_gff_generates_expected_output() {
        let mut root = Vec::new();
        let features = vec![FeaturePos {
            contig_id: 0,
            feature_id: 0,
            feature_type: "intergenic".to_string(),
            start: 0,
            end: 4,
            strand: true,
            seq: vec![1, 2, 4, 8],
        }];
        root.push(features);

        let exon_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create distribution for exon features");
        let intron_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create distribution for intron features");
        let intergenic_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create distribution for intergenic features");
        let site_mutation_dists = vec![exon_dist, intron_dist, intergenic_dist];

        let exon_mu = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create mutation distribution for exon features");
        let intron_mu = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create mutation distribution for intron features");
        let intergenic_mu = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create mutation distribution for intergenic features");
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];

        let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create recombination probability distribution");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create recombination length distribution");
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let pop = Population::new(
            root,
            1,
            site_mutation_dists,
            site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            &mut rng,
        );

        let temp_path = std::env::temp_dir().join(format!(
            "pansimnuc_pop_{}.gff3",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("system time before UNIX_EPOCH")
                .as_nanos()
        ));
        let output_path = temp_path.to_string_lossy().into_owned();
        let genome_output_path = Population::genome_output_path(&output_path, 0)
            .expect("failed to construct per-genome output path");

        pop.write_gff(&output_path)
            .expect("failed to write test GFF file");

        let mut content = String::new();
        fs::File::open(&genome_output_path)
            .expect("failed to open test GFF file")
            .read_to_string(&mut content)
            .expect("failed to read test GFF file");

        assert!(content.contains("##gff-version 3"));
        assert!(content.contains("\tintergenic\t"));
        assert!(content.contains("feature_id=0"));
        assert!(content.contains("sequence=ACGT"));
        assert!(content.contains("overall_selection_coefficient="));
        assert!(content.contains("sv_duplication_rate="));
        assert!(content.contains("sv_deletion_rate="));
        assert!(content.contains("sv_inversion_rate="));
        assert!(content.contains("sv_duplication_insertion_prob="));

        let _ = fs::remove_file(genome_output_path);
    }

    #[test]
    fn test_mutate_changes_sequence_and_mutation_map() {
        let mut root = Vec::new();
        let features = vec![FeaturePos {
            contig_id: 0,
            feature_id: 0,
            feature_type: "exon".to_string(),
            start: 0,
            end: 4,
            strand: true,
            seq: vec![1, 2, 4, 8],
        }];
        root.push(features);

        let exon_selection_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create exon selection distribution");
        let intron_selection_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create intron selection distribution");
        let intergenic_selection_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create intergenic selection distribution");
        let site_mutation_dists = vec![
            exon_selection_dist,
            intron_selection_dist,
            intergenic_selection_dist,
        ];
        let force_mutation_dist = MutationDistribution::new_poisson(100.0)
            .expect("failed to create mutation distribution");
        let intron_mu_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create intron mutation distribution");
        let intergenic_mu_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create intergenic mutation distribution");
        let site_mutation_mus = vec![force_mutation_dist, intron_mu_dist, intergenic_mu_dist];
                let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let mut pop = Population::new(
            root,
            1,
            site_mutation_dists,
            site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            &mut rng,
        );

        let original_seq = pop.pop[0].seq[0].seq.clone();
        let original_map = pop.pop[0].seq[0].mutation_map.clone();

        pop.mutate();

        let mutated_element = &pop.pop[0].seq[0];
        assert_ne!(mutated_element.seq, original_seq);

        for (site, (&old_allele, &new_allele)) in original_seq
            .iter()
            .zip(mutated_element.seq.iter())
            .enumerate()
        {
            assert_ne!(new_allele, old_allele);
            assert!(mutated_element.mutation_map.get(new_allele, site).is_some());
        }
    }

    #[test]
    fn test_next_generation_creates_new_population() {
        let mut root = Vec::new();
        let features = vec![
            FeaturePos {
                contig_id: 0,
                feature_id: 0,
                feature_type: "exon".to_string(),
                start: 0,
                end: 4,
                strand: true,
                seq: vec![1, 2, 4, 8],
            },
            FeaturePos {
                contig_id: 0,
                feature_id: 1,
                feature_type: "intron".to_string(),
                start: 4,
                end: 8,
                strand: true,
                seq: vec![8, 4, 2, 1],
            },
        ];
        root.push(features);

        let exon_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create exon selection distribution");
        let intron_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create intron selection distribution");
        let intergenic_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create intergenic selection distribution");
        let site_mutation_dists = vec![exon_dist, intron_dist, intergenic_dist];

        let exon_mu = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create exon mutation distribution");
        let intron_mu = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create intron mutation distribution");
        let intergenic_mu = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create intergenic mutation distribution");
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];

        let mut pop = Population::new(
            root,
            3,
            site_mutation_dists,
            site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            &mut rng,
        );

        let original_identifiers: Vec<String> = pop
            .pop
            .iter()
            .map(|genome| genome.identifier.clone())
            .collect();
        let original_parents: Vec<String> = pop
            .pop
            .iter()
            .map(|genome| genome.parent.clone())
            .collect();
        let original_sequences: Vec<Vec<Vec<u8>>> = pop
            .pop
            .iter()
            .map(|genome| genome.seq.iter().map(|element| element.seq.clone()).collect())
            .collect();

        pop.next_generation(vec![2, 0, 1]);

        assert_eq!(pop.generation, 1);
        assert_eq!(pop.pop.len(), 3);

        for (new_index, genome) in pop.pop.iter().enumerate() {
            let selected_index = [2usize, 0, 1][new_index];
            assert_eq!(genome.identifier, format!("1-{}", original_identifiers[selected_index]));
            assert_eq!(genome.parent, original_identifiers[selected_index]);
            assert_ne!(genome.identifier, original_identifiers[new_index]);
            assert_ne!(genome.parent, original_parents[new_index]);

            let new_sequences: Vec<Vec<u8>> = genome
                .seq
                .iter()
                .map(|element| element.seq.clone())
                .collect();
            assert_eq!(new_sequences, original_sequences[selected_index]);
        }
    }

    fn make_check_feature_order_population() -> Population {
        let mut root = Vec::new();
        let features = vec![
            FeaturePos {
                contig_id: 0,
                feature_id: 0,
                feature_type: "intergenic".to_string(),
                start: 0,
                end: 4,
                strand: true,
                seq: vec![1, 1, 1, 1],
            },
            FeaturePos {
                contig_id: 0,
                feature_id: 1,
                feature_type: "exon".to_string(),
                start: 4,
                end: 8,
                strand: true,
                seq: vec![2, 2, 2, 2],
            },
            FeaturePos {
                contig_id: 0,
                feature_id: 1,
                feature_type: "intron".to_string(),
                start: 8,
                end: 12,
                strand: true,
                seq: vec![4, 4, 4, 4],
            },
            FeaturePos {
                contig_id: 0,
                feature_id: 1,
                feature_type: "exon".to_string(),
                start: 12,
                end: 16,
                strand: true,
                seq: vec![8, 8, 8, 8],
            },
            FeaturePos {
                contig_id: 0,
                feature_id: 0,
                feature_type: "intergenic".to_string(),
                start: 16,
                end: 20,
                strand: true,
                seq: vec![1, 1, 1, 1],
            },
        ];
        root.push(features);

        let exon_dist = MutationDistribution::new_uniform(0.0, 0.1)
            .expect("failed to create exon selection distribution");
        let intron_dist = MutationDistribution::new_uniform(0.0, 0.1)
            .expect("failed to create intron selection distribution");
        let intergenic_dist = MutationDistribution::new_uniform(0.0, 0.1)
            .expect("failed to create intergenic selection distribution");
        let site_mutation_dists = vec![exon_dist, intron_dist, intergenic_dist];

        let exon_mu = MutationDistribution::new_uniform(0.0, 0.1)
            .expect("failed to create exon mutation distribution");
        let intron_mu = MutationDistribution::new_uniform(0.0, 0.1)
            .expect("failed to create intron mutation distribution");
        let intergenic_mu = MutationDistribution::new_uniform(0.0, 0.1)
            .expect("failed to create intergenic mutation distribution");
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];

        let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("failed to create recombination count distribution");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("failed to create recombination length distribution");
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];

        let mut rng: StdRng = StdRng::seed_from_u64(123);
        Population::new(
            root,
            1,
            site_mutation_dists,
            site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            &mut rng,
        )
    }

    fn genome_from_seq(seq: Vec<NucElement>) -> Genome {
        Genome {
            identifier: "test".to_string(),
            genome_id: 0,
            contig_starts: vec![0],
            parent: "test-parent".to_string(),
            seq,
        }
    }

    fn check_feature_one_intron(pop: &Population, genome: &Genome) -> (bool, f64) {
        let idx = genome
            .seq
            .iter()
            .position(|element| element.feature_id == 1 && element.feature_type == "intron")
            .expect("expected intron in feature block");
        pop.check_feature_order(genome, idx, &genome.seq[idx])
    }

    #[test]
    fn test_check_feature_order_identifies_unbroken_gene() {
        let pop = make_check_feature_order_population();
        let genome = &pop.pop[0];
        let (broken, multiplier) = check_feature_one_intron(&pop, genome);
        assert!(!broken);
        assert!((multiplier - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_check_feature_order_identifies_gene_broken_by_insertion() {
        let pop = make_check_feature_order_population();
        let mut seq = pop.pop[0].seq.clone();

        let mut inserted_te = seq[0].clone();
        inserted_te.feature_type = "TE-CUT".to_string();
        inserted_te.feature_id = 0;
        inserted_te.multiplier = 2.0;
        inserted_te.element_id = 10_000;
        seq.insert(2, inserted_te);

        let genome = genome_from_seq(seq);
        let (broken, _) = check_feature_one_intron(&pop, &genome);
        assert!(broken);
    }

    #[test]
    fn test_check_feature_order_identifies_gene_broken_by_deletion() {
        let pop = make_check_feature_order_population();
        let mut seq = pop.pop[0].seq.clone();
        seq.retain(|element| element.element_id != 2);

        let genome = genome_from_seq(seq);
        let exon_idx = genome
            .seq
            .iter()
            .position(|element| element.feature_id == 1 && element.feature_type == "exon")
            .expect("expected exon in feature block");
        let (broken, _) = pop.check_feature_order(&genome, exon_idx, &genome.seq[exon_idx]);
        assert!(broken);
    }

    #[test]
    fn test_check_feature_order_identifies_gene_broken_by_partial_inversion() {
        let pop = make_check_feature_order_population();
        let mut seq = pop.pop[0].seq.clone();
        let intron_idx = seq
            .iter()
            .position(|element| element.element_id == 2)
            .expect("expected intron with element_id=2");
        seq[intron_idx].strand = !seq[intron_idx].strand;

        let genome = genome_from_seq(seq);
        let (broken, _) = check_feature_one_intron(&pop, &genome);
        assert!(broken);
    }

    #[test]
    fn test_check_feature_order_identifies_whole_gene_reversal() {
        let pop = make_check_feature_order_population();
        let base_seq = &pop.pop[0].seq;

        let left_flank = base_seq[0].clone();
        let right_flank = base_seq[4].clone();
        let mut reversed_gene = base_seq[1..=3].to_vec();
        reversed_gene.reverse();
        for element in &mut reversed_gene {
            element.strand = !element.strand;
        }

        let mut seq = vec![left_flank];
        seq.extend(reversed_gene);
        seq.push(right_flank);

        let genome = genome_from_seq(seq);
        let (broken, multiplier) = check_feature_one_intron(&pop, &genome);
        assert!(!broken);
        assert!((multiplier - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_check_feature_order_identifies_te_upstream_and_downstream() {
        let pop = make_check_feature_order_population();
        let mut seq = pop.pop[0].seq.clone();

        let mut upstream_te = seq[0].clone();
        upstream_te.feature_type = "TE-CUT".to_string();
        upstream_te.feature_id = 0;
        upstream_te.multiplier = 2.0;
        upstream_te.element_id = 20_000;
        seq.insert(1, upstream_te);

        let mut downstream_te = seq[0].clone();
        downstream_te.feature_type = "TE-COPY".to_string();
        downstream_te.feature_id = 0;
        downstream_te.multiplier = 3.5;
        downstream_te.element_id = 20_001;
        seq.insert(5, downstream_te);

        let genome = genome_from_seq(seq);
        let (broken, multiplier) = check_feature_one_intron(&pop, &genome);
        assert!(!broken);
        assert!(multiplier == 3.5);
    }
}
