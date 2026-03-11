use crate::gff::FeaturePos;
use crate::mutation::MutationMap;
use crate::mutation::Distribution as MutationDistribution;
use crate::structural::mutate_intra_genome;
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
    pub seqname: String,
    pub element_id: usize,
    pub feature_id: usize,
    pub feature_type: String,
    pub selection_coefficient: f64,
    pub seq: Vec<u8>,
    pub mutation_map: MutationMap,
    pub strand: bool,
    pub structure_mutation_map: StructureMutationMap,
}

pub struct Genome {
    pub identifier: String,
    pub genome_id: usize,
    pub parent: String,
    pub seq: Vec<NucElement>
}

pub struct Population{
    pub generation: usize,
    pub pop: Vec<Genome>,
    pub core_vec: Vec<Vec<u8>>,
    pub selection_dists: Vec<MutationDistribution>,
    pub mu_dists: Vec<MutationDistribution>,
    // Map from original element ID to positions of homologous regions in other genomes
    //outermost loop is the homology group, middle loop is genomes, inner loop is positions
    pub homology_map: Vec<Vec<Vec<usize>>>
}

impl Population {
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

    pub fn new(
        root: HashMap<String, Vec<FeaturePos>>,
        n_genomes: usize,
        selection_dists: Vec<MutationDistribution>,
        mu_dists: Vec<MutationDistribution>,
        rng: &mut StdRng,
    ) -> Self {
        // initialise population
        let mut population: Vec<Genome> = Vec::new();
        let mut genome: Vec<NucElement> = Vec::new();

        // count for element ID, each NucElement gets own to signal it it's homology group
        let mut element_id: usize = 0;

        // initialise homology map, outermost loop is the homology group, middle loop is genomes, inner loop is positions
        let mut homology_map: Vec<Vec<Vec<usize>>> = Vec::new();

        // generate starting genome
        for (seqname, features) in &root {
            for feature in features {
                
                // TODO change this so that can specify different mutation rates per site
                // also add separate TE compartment
                let selection_dist_id:usize = match feature.feature_type.as_str() {
                    "exon" => 0,
                    "intron" => 1,
                    "intergenic" => 2,
                    _ => panic!("Unknown feature type: {}", feature.feature_type),
                };
                
                let mu_dist_id:usize = match feature.feature_type.as_str() {
                    "exon" => 0,
                    "intron" => 1,
                    "intergenic" => 2,
                    _ => panic!("Unknown feature type: {}", feature.feature_type),
                };

                genome.push(NucElement {
                    seqname: seqname.clone(),
                    element_id: element_id,
                    feature_id: feature.feature_id,
                    feature_type: feature.feature_type.clone(),
                    seq: feature.seq.clone(),
                    strand: feature.strand,
                    mutation_map: MutationMap::new(selection_dist_id, mu_dist_id, &feature.seq, &selection_dists[selection_dist_id], rng),
                    structure_mutation_map: StructureMutationMap {
                        duplication_rate: 0.0,
                        deletion_rate: 0.0,
                        inversion_rate: 0.0,
                        recombination_rate: 0.0,
                        max_duplications: None,
                        duplication_insertion_prob: 0.5,
                    },
                    selection_coefficient: 0.0, // Initialize with a default value, can be updated later
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
            population.push(Genome {
                identifier: format!("{}", i),
                genome_id: i,
                parent: "root".to_string(),
                seq: genome.clone(),
            });
        }
        
        let core_vec: Vec<Vec<u8>> =
            vec![vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4]];

        Self {
            generation: 0,
            pop: population,
            core_vec,
            selection_dists,
            mu_dists,
            homology_map
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
        for idx in 0..self.pop.len() {
            let genome = &mut self.pop[idx];

            // duplications
            mutate_intra_genome(genome, &mut self.homology_map, &duplication_mu_dist, &duplication_pos_dist);

            // translocations
            // clear homology map to enable fresh creation of groups
            self.homology_map.iter_mut().for_each(|element_homology_map| {
                element_homology_map[idx].clear();
            });
            mutate_intra_genome(genome, &mut self.homology_map, &translocation_mu_dist, &translocation_pos_dist);
        }
    }

    // sample individuals using logsumexp normalisation to prevent underflow/overflow issues with very small/large weights
    pub fn sample_individuals (&mut self, rng: &mut StdRng) -> Vec<usize> {
        let mut selection_weights: Vec<f64> = vec![1.0; self.pop.len()];
        
        selection_weights = self.pop
            .par_iter_mut()
            .map(|genome| {
                let mut log_sum = 0.0;
                for element in &mut genome.seq {
                    let mut element_log_sum = 0.0;
                    for (site, allele) in element.seq.iter().enumerate() {
                        let allele_shifted = 1 >> allele;
                        if let Some(coeff) = element.mutation_map.get(*allele, site) {
                            log_sum += coeff;
                            element_log_sum += coeff;
                        } else {
                            panic!(
                                "Failed to generate selection coefficient for genome {} allele {} (shifted {}) at site {}",
                                genome.identifier, allele, allele_shifted, site
                            );
                        }
                    }
                    
                    element.selection_coefficient = element_log_sum;
                }
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
            let mut seqname_groups: HashMap<String, Vec<usize>> = HashMap::new();
            for (idx, element) in genome.seq.iter().enumerate() {
                seqname_groups.entry(element.seqname.clone())
                    .or_insert_with(Vec::new)
                    .push(idx);
            }

            // Write each seqname group as a separate FASTA entry
            for (seqname, indices) in seqname_groups {
                writeln!(
                    writer,
                    ">{id}_{seqname} parent={parent} generation={generation}",
                    id = genome.identifier,
                    seqname = seqname,
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gff::FeaturePos;
    use std::fs;
    use std::io::Read;

    #[test]
    fn test_population_new() {
        // Create test data
        let mut root = HashMap::new();
        let features = vec![
            FeaturePos {
                seqname: "chr1".to_string(),
                feature_id: 0,
                feature_type: "exon".to_string(),
                start: 100,
                end: 200,
                strand: true,
                seq: vec![1, 2, 4, 8], // ACGT
            },
            FeaturePos {
                seqname: "chr1".to_string(),
                feature_id: 1,
                feature_type: "intron".to_string(),
                start: 300,
                end: 400,
                strand: false,
                seq: vec![8, 4, 2, 1], // TGCA
            },
        ];
        root.insert("chr1".to_string(), features);

        let n_genomes = 3;
        let exon_dist = MutationDistribution::new_double_exp(0.5, 2.0, 0.3).expect("Failed to create double exponential distribution for exon features");
        let intron_dist = MutationDistribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intron features");
        let intergenic_dist = MutationDistribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intergenic features");
        let site_mutation_dists = vec![exon_dist, intron_dist, intergenic_dist];

        let exon_mu = MutationDistribution::new_uniform(0.0, 1.0).expect("Failed to create double exponential distribution for exon features");
        let intron_mu = MutationDistribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intron features");
        let intergenic_mu = MutationDistribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intergenic features");
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let pop = Population::new(root, n_genomes, site_mutation_dists, site_mutation_mus, &mut rng);

        // Check population was created correctly
        assert_eq!(pop.generation, 0);
        assert_eq!(pop.pop.len(), n_genomes);

        // Check each genome
        for (i, genome) in pop.pop.iter().enumerate() {
            assert_eq!(genome.identifier, format!("{}", i));
            assert_eq!(genome.parent, "root");
            assert_eq!(genome.seq.len(), 2); // Two features

            // Check first feature
            assert_eq!(genome.seq[0].seqname, "chr1");
            assert_eq!(genome.seq[0].feature_id, 0);
            assert_eq!(genome.seq[0].feature_type, "exon");

            // Check second feature
            assert_eq!(genome.seq[1].seqname, "chr1");
            assert_eq!(genome.seq[1].feature_id, 1);
            assert_eq!(genome.seq[1].feature_type, "intron");
        }
    }

    #[test]
    fn test_write_fasta_decodes_nucleotides() {
        let mut root = HashMap::new();
        let features = vec![FeaturePos {
            seqname: "chr1".to_string(),
            feature_id: 0,
            feature_type: "exon".to_string(),
            start: 0,
            end: 4,
            strand: true,
            seq: vec![1, 2, 4, 8],
        }];
        root.insert("chr1".to_string(), features);

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

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let pop = Population::new(root, 1, site_mutation_dists, site_mutation_mus, &mut rng);

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

        assert!(content.contains(">0_chr1 parent=root generation=0"));
        assert!(content.contains("ACGT"));

        let _ = fs::remove_file(genome_output_path);
    }

    #[test]
    fn test_mutate_changes_sequence_and_mutation_map() {
        let mut root = HashMap::new();
        let features = vec![FeaturePos {
            seqname: "chr1".to_string(),
            feature_id: 0,
            feature_type: "exon".to_string(),
            start: 0,
            end: 4,
            strand: true,
            seq: vec![1, 2, 4, 8],
        }];
        root.insert("chr1".to_string(), features);

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

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let mut pop = Population::new(root, 1, site_mutation_dists, site_mutation_mus, &mut rng);

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
        let mut root = HashMap::new();
        let features = vec![
            FeaturePos {
                seqname: "chr1".to_string(),
                feature_id: 0,
                feature_type: "exon".to_string(),
                start: 0,
                end: 4,
                strand: true,
                seq: vec![1, 2, 4, 8],
            },
            FeaturePos {
                seqname: "chr1".to_string(),
                feature_id: 1,
                feature_type: "intron".to_string(),
                start: 4,
                end: 8,
                strand: true,
                seq: vec![8, 4, 2, 1],
            },
        ];
        root.insert("chr1".to_string(), features);

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
        let mut pop = Population::new(root, 3, site_mutation_dists, site_mutation_mus, &mut rng);

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
}
