use crate::gff::FeaturePos;
use crate::mutation::Distribution as MutationDistribution;
use crate::mutation::MutationMap;
use crate::structural::mutate_inter_genome;
use crate::structural::mutate_intra_genome;
use crate::tracking::identify_tracked_element;
use logsumexp::LogSumExp;
use rand::SeedableRng;
use rand::distributions::{Distribution as RandDistribution, WeightedIndex};
use rand::rngs::StdRng;
use rand::rngs::ThreadRng;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

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
    pub original_length: usize,
    pub frameshift: bool,
    pub tracked: bool
}

impl NucElement {
    fn element_selection_coefficient(&self, genome_identifier: &str) -> f64 {
        let mut element_log_sum = 0.0;

        for (site, allele) in self.seq.iter().enumerate() {
            let allele_shifted = 1 >> allele;
            if let Some(coeff) = self.mutation_map.get(*allele, site) {
                let log_coeff = (1.0 + coeff).ln(); // add log of coefficient to log sum
                if log_coeff == std::f64::NEG_INFINITY {
                    // if coefficient is -1, set log sum to -inf and break loop, as any other mutations won't change this
                    element_log_sum = std::f64::NEG_INFINITY;
                    break;
                }
                element_log_sum += log_coeff;
            } else {
                panic!(
                    "Failed to generate selection coefficient for genome {} allele {} (shifted {}) at site {}",
                    genome_identifier, allele, allele_shifted, site
                );
            }
        }
        element_log_sum
    }

    pub fn generate_selection_coefficients(&self) -> Vec<f64> {
        self.seq
            .iter()
            .enumerate()
            .map(|(site, allele)| {
                if let Some(coeff) = self.mutation_map.get(*allele, site) {
                    let log_coeff = (1.0 + coeff).ln(); // add log of coefficient to log sum
                    log_coeff
                } else {
                    panic!(
                        "Failed to generate selection coefficient for allele {} at site {}",
                        allele, site
                    );
                }
            })
            .collect()
    }
}

#[derive(Clone)]
pub struct Genome {
    pub identifier: String,
    pub genome_id: usize,
    pub contig_starts: Vec<usize>,
    pub parent: String,
    pub seq: Vec<NucElement>,
    pub seq_length: usize,
}

impl Genome {
    pub fn update_contig_starts(&mut self) {
        self.contig_starts.clear();
        let mut total_length = 0;
        for (idx, element) in self.seq.iter().enumerate() {
            if idx == 0 || element.contig_id != self.seq[idx - 1].contig_id {
                self.contig_starts.push(idx);
            }
            total_length += element.seq.len();
        }
        self.seq_length = total_length;
    }
}

#[derive(Clone)]
pub struct Population {
    pub id: usize,
    pub generation: usize,
    pub pop: Vec<Genome>,
    pub core_vec: Vec<Vec<u8>>,
    pub selection_dists: Vec<MutationDistribution>,
    pub mu_dists: Vec<MutationDistribution>,
    pub indel_dists: Vec<MutationDistribution>,
    pub structural_mu_dists: Vec<Vec<MutationDistribution>>,
    pub recombination_dists: Vec<MutationDistribution>,
    pub recombination_threshold: f64,
    pub homology_map: Vec<Vec<Vec<usize>>>, // Map from original element ID to positions of homologous regions in other genomes, outermost loop is the homology group, middle loop is genomes, inner loop is positions
    pub feature_map: HashMap<usize, Vec<usize>>, // Map from feature ID to genes that share same ID
    pub max_multiplier_dist: usize,
    pub n_generations: usize,
    pub verbose: bool,
    pub augment_tracking: bool,
    pub genome_size_penalty_per_bp: f64,
    pub optimal_genome_size: usize,
}

impl Population {
    fn total_seq_length(&self) -> usize {
        self.pop.iter().map(|genome| genome.seq_length).sum()
    }

    fn is_te_feature(feature_type: &str) -> bool {
        feature_type == "TE-CUT" || feature_type == "TE-COPY"
    }

    fn check_feature_order(
        &self,
        genome: &Genome,
        element_idx: usize,
        element: &NucElement,
    ) -> (bool, f64) {
        let mut feature_broken = false;
        let mut feature_multiplier = 1.0;

        let max_multiplier_dist = self.max_multiplier_dist;

        // identify regions where order matters
        if element.feature_type == "exon" || element.feature_type == "intron" {
            let feature_map_entry = self
                .feature_map
                .get(&element.feature_id)
                .expect("Entry missing from feature_map");
            let feature_map_entry_len = feature_map_entry.len();

            // get position of element in feature_map_entry
            let position = feature_map_entry
                .iter()
                .position(|n| n == &element.element_id)
                .expect("Element not found in feature_map_entry");

            // check upstream elements in feature_map_entry
            for feature_element_idx in 0..position {
                // get upstream feature in genome
                let actual_element = &genome.seq[element_idx - (position - feature_element_idx)];
                let actual_element_id = actual_element.element_id;

                // get expected element
                let expected_element_id = feature_map_entry[feature_element_idx];

                // expected element and strand matches so continue
                if actual_element_id == expected_element_id
                    && actual_element.strand == element.strand
                {
                    // check if frameshift occurred, if so then feature is broken, as likely to be non-functional
                    if actual_element.feature_type == "exon" && actual_element.frameshift {
                        feature_broken = true;
                        break;
                    } else {
                        continue;
                    }
                } else {
                    // check if reversed order and strand matches, if so continue, as likely to be functional just reversed
                    // check if last feature matches in feature_map_entry
                    let reversed_feature_id =
                        feature_map_entry[(feature_map_entry_len - 1) - feature_element_idx];
                    if actual_element_id == reversed_feature_id
                        && actual_element.strand == element.strand
                    {
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
                for feature_element_idx in (position + 1)..feature_map_entry_len {
                    // check if at end of genome, if so then feature is broken, as missing downstream elements
                    if element_idx + (feature_element_idx - position) >= genome.seq.len() {
                        feature_broken = true;
                        break;
                    }

                    // get downstream feature in genome
                    let actual_element =
                        &genome.seq[element_idx + (feature_element_idx - position)];
                    let actual_element_id = actual_element.element_id;

                    // get expected element
                    let expected_element_id = feature_map_entry[feature_element_idx];

                    // expected element and strand matches so continue
                    if actual_element_id == expected_element_id
                        && actual_element.strand == element.strand
                    {
                        // check if frameshift occurred, if so then feature is broken, as likely to be non-functional
                        if actual_element.feature_type == "exon" && actual_element.frameshift {
                            feature_broken = true;
                            break;
                        } else {
                            continue;
                        }
                    } else {
                        // check if reversed order and strand matches, if so continue, as likely to be functional just reversed
                        // check if last feature matches in feature_map_entry
                        let reversed_feature_id =
                            feature_map_entry[(feature_map_entry_len - 1) - feature_element_idx];
                        if actual_element_id == reversed_feature_id
                            && actual_element.strand == element.strand
                        {
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
                    let upstream_size = upstream_element.seq.len();
                    if Self::is_te_feature(&upstream_element.feature_type) {
                        feature_multiplier = upstream_element.multiplier;
                    }
                    // check two down and see if TE
                    else if element_idx >= position + 2 && upstream_size <= max_multiplier_dist {
                        let upstream_element2 = &genome.seq[element_idx - (position + 2)];
                        if Self::is_te_feature(&upstream_element2.feature_type) {
                            if feature_multiplier.abs() < upstream_element2.multiplier.abs() {
                                feature_multiplier = upstream_element2.multiplier;
                            }
                        }
                    }
                }

                // check not at end of genome
                if element_idx + (feature_map_entry_len - position) < genome.seq.len() {
                    let downstream_element =
                        &genome.seq[element_idx + (feature_map_entry_len - position)];
                    let downstream_size = downstream_element.seq.len();
                    if Self::is_te_feature(&downstream_element.feature_type) {
                        if feature_multiplier.abs() < downstream_element.multiplier.abs() {
                            feature_multiplier = downstream_element.multiplier;
                        }
                    }
                    // check two down and see if TE
                    else if element_idx + ((feature_map_entry_len - position) + 1)
                        < genome.seq.len()
                        && downstream_size <= max_multiplier_dist
                    {
                        let downstream_element2 =
                            &genome.seq[element_idx + ((feature_map_entry_len - position) + 1)];
                        if Self::is_te_feature(&downstream_element2.feature_type) {
                            if feature_multiplier.abs() < downstream_element2.multiplier.abs() {
                                feature_multiplier = downstream_element2.multiplier;
                            }
                        }
                    }
                }
            }
        }

        (feature_broken, feature_multiplier)
    }

    fn genome_output_path(output_path: &str, prefix: &str) -> io::Result<PathBuf> {
        let path = Path::new(output_path);
        let file_name = path.file_name().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Output path must include a file name: {output_path}"),
            )
        })?;

        let prefixed_name = format!("{}{}", prefix, file_name.to_string_lossy());
        Ok(path.with_file_name(prefixed_name))
    }

    pub fn decode_base(base: u8) -> u8 {
        match base {
            1 => b'A',
            2 => b'C',
            4 => b'G',
            8 => b'T',
            16 => b'N',
            _ => panic!("Invalid base encoding: {}", base),
        }
    }

    fn genome_selection_coefficient(&self, genome: &Genome) -> f64 {
        let mut log_sum = 0.0;

        for (element_idx, element) in genome.seq.iter().enumerate() {
            let mut element_log_sum = 0.0;

            let (feature_broken, feature_multiplier) =
                self.check_feature_order(genome, element_idx, element);

            // if feature not broken, add up sites
            if !feature_broken {
                element_log_sum += element.element_selection_coefficient(&genome.identifier);
                
                // if any value is zero, product is zero so genome selection coefficient is zero
                if element_log_sum == std::f64::NEG_INFINITY {
                    log_sum = std::f64::NEG_INFINITY;
                    break; // if no contribution to selection coefficient, skip multiplier calculation to avoid unnecessary calculations and potential floating point issues
                }
            }

            log_sum += element_log_sum * feature_multiplier;
        }
        log_sum
    }

    fn log_sum_exp(&self) -> (Vec<f64>, f64) {
        let selection_weights = self
            .pop
            .par_iter()
            .map(|genome| {
                let mut log_sum = self.genome_selection_coefficient(genome);

                if log_sum == std::f64::NEG_INFINITY {
                    // if selection coefficient is zero, set log sum to -inf to prevent underflow issues when exponentiating later, and skip multiplier calculation to avoid unnecessary calculations and potential floating point issues
                    log_sum = 0.0;
                }
                log_sum
            })
            .collect::<Vec<f64>>();

        // logsumexp normalization to prevent underflow/overflow issues with very small/large weights
        let logsumexp_value = selection_weights.iter().ln_sum_exp();
        (selection_weights, logsumexp_value)
    }

    pub fn update_homology_map(&mut self) {
        // update homology map for all new elements
        for genome in &self.pop {
            self.homology_map
                .iter_mut()
                .for_each(|element_homology_map| {
                    element_homology_map[genome.genome_id].clear();
                });

            for (element_idx, element) in genome.seq.iter().enumerate() {
                let element_id = element.element_id;
                let homology_group = &mut self.homology_map[element_id][genome.genome_id];
                homology_group.push(element_idx); // convert back to 0 indexed
            }
        }
    }

    pub fn new(
        root: Vec<Vec<FeaturePos>>,
        n_genomes: usize,
        selection_dists: Vec<MutationDistribution>,
        mu_dist_vals: &Vec<f64>,
        indel_dist_vals: &Vec<f64>,
        recombination_dists: Vec<MutationDistribution>,
        recombination_threshold: f64,
        structural_mu_dists: Vec<Vec<MutationDistribution>>,
        max_multiplier_dist: usize,
        multiplier_dists: Vec<MutationDistribution>,
        n_generations: usize,
        rng: &mut StdRng,
        verbose: bool,
        contig_name_to_id: &Vec<String>,
        tracking_regions: &Vec<(String, usize, usize)>,
        augment_tracking: bool,
        genome_size_penalty_per_bp: f64,
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

        // determine if tracking enabled
        let is_tracking = !tracking_regions.is_empty();

        // track optimal genome size
        let mut optimal_genome_size = 0;

        // generate starting genome
        for (contig_id, features) in root.iter().enumerate() {
            let mut current_start = 0;
            for feature in features {
                let selection_dist_id: usize = match feature.feature_type.as_str() {
                    "exon" => 0,
                    "intron" => 1,
                    "intergenic" => 2,
                    "TE-CUT" => 3,
                    "TE-COPY" => 4,
                    _ => panic!("Unknown feature type: {}", feature.feature_type),
                };

                let mu_dist_id: usize = match feature.feature_type.as_str() {
                    "exon" => 0,
                    "intron" => 1,
                    "intergenic" => 2,
                    "TE-CUT" => 3,
                    "TE-COPY" => 4,
                    _ => panic!("Unknown feature type: {}", feature.feature_type),
                };

                let multiplier_dist: &MutationDistribution = match feature.feature_type.as_str() {
                    "exon" => &multiplier_dists[0],
                    "intron" => &multiplier_dists[1],
                    "intergenic" => &multiplier_dists[2],
                    "TE-CUT" => &multiplier_dists[3],
                    "TE-COPY" => &multiplier_dists[4],
                    _ => panic!("Unknown feature type: {}", feature.feature_type),
                };

                let multiplier = match feature.feature_type.as_str() {
                    "exon" => 1.0,
                    "intron" => 1.0,
                    "intergenic" => 1.0,
                    "TE-CUT" => multiplier_dist.sample(rng),
                    "TE-COPY" => multiplier_dist.sample(rng),
                    _ => panic!("Unknown feature type: {}", feature.feature_type),
                };

                // Update feature_map, keeping track of how many features there are
                if feature.feature_id != 0 {
                    feature_map
                        .entry(feature.feature_id)
                        .or_default()
                        .push(element_id);
                }

                let mut element = NucElement {
                    contig_id: contig_id,
                    element_id: element_id,
                    feature_id: feature.feature_id,
                    feature_type: feature.feature_type.clone(),
                    seq: feature.seq.clone(),
                    strand: feature.strand,
                    mutation_map: MutationMap::new(
                        selection_dist_id,
                        mu_dist_id,
                        &feature.seq,
                        &selection_dists[selection_dist_id],
                        rng,
                    ),
                    multiplier: multiplier,
                    original_length: feature.seq.len(),
                    frameshift: false,
                    tracked: false,
                };

                // determine if element is tracked and if so update mutation map accordingly
                if is_tracking {
                    identify_tracked_element(&mut element, current_start, &tracking_regions, &contig_name_to_id);
                    if element.tracked && augment_tracking {
                        // update selection maps etc
                        element.mutation_map = MutationMap::new(
                            5, // new selection dist ID for tracked elements
                            5, // new mu dist ID for tracked elements
                            &element.seq,
                            &selection_dists[5], // new selection distribution for tracked elements
                            rng,
                        );
                        element.multiplier = multiplier_dists[5].sample(rng);
                    }
                }

                genome.push(element);
                element_id += 1;
                current_start += feature.seq.len();

                optimal_genome_size += feature.seq.len();

                // generate homology map for this element, initially just self
                let mut element_homology_map: Vec<Vec<usize>> = Vec::new();
                for _ in 0..n_genomes {
                    element_homology_map.push(vec![element_id]);
                }
                homology_map.push(element_homology_map);
            }
        }

        // copy whole genome to start
        let mut total_length = 0;
        for i in 0..n_genomes {
            let mut genome_entry = Genome {
                identifier: format!("{}", i),
                genome_id: i,
                contig_starts: Vec::new(), // will be updated after mutations
                parent: "root".to_string(),
                seq: genome.clone(),
                seq_length: 0, // will be updated after mutations
            };
            genome_entry.update_contig_starts();
            total_length += genome_entry.seq_length;
            population.push(genome_entry);
        }

        // generate mu_dists based on total population size and total sequence length, so that mutation rates are per base per genome generation
        let mu_dists = mu_dist_vals
            .into_iter()
            .map(|mu| {
                MutationDistribution::new_poisson(mu * (total_length as f64) * (n_genomes as f64) * n_generations as f64)
                    .expect("Failed to create poisson distribution for mutation rates")
            })
            .collect();
            
        let indel_dists = indel_dist_vals
            .into_iter()
            .map(|mu| {
                MutationDistribution::new_poisson(mu * (total_length as f64) * (n_genomes as f64) * n_generations as f64)
                    .expect("Failed to create poisson distribution for indel rates")
            })
            .collect();

        let core_vec: Vec<Vec<u8>> =
            vec![vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4], vec![1, 2, 4, 8, 16]];

        Self {
            id: 0,
            generation: 0,
            pop: population,
            core_vec,
            selection_dists,
            mu_dists,
            indel_dists,
            structural_mu_dists: structural_mu_dists,
            recombination_dists,
            recombination_threshold,
            homology_map,
            feature_map,
            max_multiplier_dist,
            n_generations,
            verbose,
            augment_tracking,
            genome_size_penalty_per_bp,
            optimal_genome_size,
        }
    }

    // mutate individuals in the population according to their mutation maps and the provided distributions
    pub fn mutate(&mut self) -> (usize, usize) {
        let core_vec = &self.core_vec;
        let selection_dists = &self.selection_dists;
        let mu_dists = &self.mu_dists;
        let indel_dists = &self.indel_dists;

        let (total_snps, total_indels): (usize, usize) = self
            .pop
            .par_iter_mut()
            .map(|genome| {
                genome
                    .seq
                    .iter_mut()
                    .fold((0, 0), |(snps, indels), element| {
                        let (s, i) = element.mutation_map.mutate(
                            core_vec,
                            &mut element.seq,
                            element.original_length,
                            &mut element.frameshift,
                            &selection_dists[element.mutation_map.selection_dist_id],
                            &mu_dists[element.mutation_map.mu_dist_id],
                            &indel_dists[element.mutation_map.mu_dist_id],
                        );
                        (snps + s, indels + i)
                    })
            })
            .reduce(|| (0, 0), |(a0, a1), (b0, b1)| (a0 + b0, a1 + b1));

        if self.verbose {
            println!("Total SNPs: {}", total_snps);
            println!("Total indels: {}", total_indels);
        }
        (total_snps, total_indels)
    }

    pub fn update_mu_dists(&mut self, mu_dist_vals: &Vec<f64>) {
        let total_length = self.total_seq_length() as f64;
        let n_genomes = self.pop.len() as f64;
        let n_generations = self.n_generations as f64;

        let new_mu_dists: Vec<MutationDistribution> = mu_dist_vals
            .into_iter()
            .map(|mu| {
                MutationDistribution::new_poisson(mu * total_length * n_genomes * n_generations)
                    .expect("Failed to create poisson distribution for mutation rates")
            })
            .collect();
        self.mu_dists = new_mu_dists;
    }

    pub fn structural_intra_genome(&mut self) {
        // probabilities for structural variations
        let pos_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create poisson distribution for structural variations");

        // iterate over genomes
        let totals = self
            .pop
            .par_iter_mut()
            .map(|genome| mutate_intra_genome(genome, &self.structural_mu_dists, &pos_dist, self.augment_tracking))
            .reduce(
                || (0usize, 0usize, 0usize, 0usize, 0usize, 0usize, 0usize),
                |a, b| (
                    a.0 + b.0, // total_non_te_duplications
                    a.1 + b.1, // total_non_te_deletions
                    a.2 + b.2, // te_cut_duplications
                    a.3 + b.3, // te_copy_duplications
                    a.4 + b.4, // te_cut_deletions
                    a.5 + b.5, // te_copy_deletions
                    a.6 + b.6, // total_inversions
                ),
            );

        let (
            total_non_te_duplications,
            total_non_te_deletions,
            te_cut_duplications,
            te_copy_duplications,
            te_cut_deletions,
            te_copy_deletions,
            total_inversions,
        ) = totals;
        
        if self.verbose { 
            println!("Total non-TE duplications: {}", total_non_te_duplications);
            println!("Total non-TE deletions: {}", total_non_te_deletions);
            println!("Total TE-CUT duplications: {}", te_cut_duplications);
            println!("Total TE-COPY duplications: {}", te_copy_duplications);
            println!("Total TE-CUT deletions: {}", te_cut_deletions);
            println!("Total TE-COPY deletions: {}", te_copy_deletions);
            println!("Total inversions: {}", total_inversions);
        }

        // update homology map for all new elements
        self.update_homology_map();
    }

    pub fn structural_inter_genome(&mut self, recombination_rate: f64, total_sites: usize, recombination_size_mean: f64) {
        // generate recombination distributions
        let average_recombinations_per_generation = (recombination_rate * total_sites as f64) / recombination_size_mean;
        if self.verbose {
            println!("Average recombinations per generation: {}", average_recombinations_per_generation);
        }
        self.recombination_dists[0] = MutationDistribution::new_poisson(average_recombinations_per_generation)
            .expect("Failed to create poisson distribution for recombination rates");

        let (n_recombinations, total_donor_length, total_recipient_length) = mutate_inter_genome(self);
        if self.verbose {
            println!("Total recombinations: {}", n_recombinations);
            println!("Total donor length: {}", total_donor_length);
            println!("Total recipient length: {}", total_recipient_length);
        }
    }

    // sample individuals using logsumexp normalisation to prevent underflow/overflow issues with very small/large weights
    pub fn sample_individuals(&mut self, rng: &mut ThreadRng) -> Vec<usize> {
        let (mut selection_weights, logsumexp_value) = self.log_sum_exp();

        #[cfg(debug_assertions)]
        {
            eprintln!("Selection pre-norm weights: {:?}", selection_weights);
        }

        selection_weights = selection_weights
            .into_iter()
            .map(|x| (x - logsumexp_value).exp()) // exp(log(w) - logsumexp)
            .collect();

        #[cfg(debug_assertions)]
        {
            eprintln!("Selection post-norm weights: {:?}", selection_weights);
        }

        // update selection weights by genome size penalty
        selection_weights = selection_weights
            .into_iter()
            .enumerate()
            .map(|(i, w)| {
                let genome_size = self.pop[i].seq_length;
                let size_penalty = self.genome_size_penalty_per_bp * ((genome_size as isize - self.optimal_genome_size as isize).abs() as f64);
                (w * size_penalty).max(std::f64::NEG_INFINITY) // ensure weights don't become negative due to penalty
            })
            .collect();

        let sum_weights: f64 = selection_weights.iter().sum();
        selection_weights = selection_weights
            .iter()
            .map(|&w| {
                if w != std::f64::NEG_INFINITY {
                    w / sum_weights
                } else {
                    0.0
                }
            })
            .collect();


        // Create a WeightedIndex distribution based on weights
        let sampling_dist = WeightedIndex::new(&selection_weights)
            .expect("Failed to generate sampling index for population.");

        // Sample rows based on the distribution
        let sampled_indices: Vec<usize> = (0..self.pop.len())
            .map(|_| sampling_dist.sample(rng))
            .collect();

        #[cfg(debug_assertions)]
        {
            println!("Selection weights: {:?}", selection_weights);
            println!("Sampled indices: {:?}", sampled_indices);
        }

        sampled_indices
    }

    pub fn next_generation(&mut self, sampled_indices: Vec<usize>) {
        let new_pop: Vec<Genome> = sampled_indices
            .par_iter()
            .enumerate()
            .map(|(genome_id, &selected_index)| {
                let selected_genome = &self.pop[selected_index];
                Genome {
                    identifier: format!("{}-{}", self.generation + 1, selected_genome.identifier),
                    genome_id,
                    contig_starts: selected_genome.contig_starts.clone(),
                    parent: selected_genome.identifier.clone(),
                    seq: selected_genome.seq.clone(),
                    seq_length: selected_genome.seq_length,
                }
            })
            .collect();

        let new_homology_map: Vec<Vec<Vec<usize>>> = self
            .homology_map
            .par_iter()
            .map(|element_homology_map| {
                let mut new_element_homology_map = Vec::with_capacity(sampled_indices.len());
                for &selected_index in &sampled_indices {
                    new_element_homology_map.push(element_homology_map[selected_index].clone());
                }
                new_element_homology_map
            })
            .collect();

        self.pop = new_pop;
        self.homology_map = new_homology_map;
        self.generation += 1;
    }

    pub fn write_fasta(&self, output_path: &str, root_genome: bool) -> io::Result<()> {
        for (genome_index, genome) in self.pop.iter().enumerate() {
            let prefix = if root_genome { "root".to_string() } else { format!("pop_{}_gen_{}_genome_{}", self.id, self.generation, genome_index) };
            let genome_output_path = Self::genome_output_path(output_path, &prefix)?;
            let file = File::create(&genome_output_path)?;
            let mut writer = BufWriter::new(file);

            // Group element indices by seqname
            let mut contig_groups: HashMap<usize, Vec<usize>> = HashMap::new();
            for (idx, element) in genome.seq.iter().enumerate() {
                contig_groups
                    .entry(element.contig_id)
                    .or_insert_with(Vec::new)
                    .push(idx);
            }

            // Write each contig group as a separate FASTA entry, in ascending contig order
            let mut sorted_contig_groups: Vec<(usize, Vec<usize>)> = contig_groups.into_iter().collect();
            sorted_contig_groups.sort_by_key(|&(id, _)| id);
            for (contig_id, indices) in sorted_contig_groups {
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

            if root_genome {
                break; // only write root genome if specified
            }
        }

        Ok(())
    }

    pub fn write_gff(&self, output_path: &str, root_genome: bool) -> io::Result<()> {
        // calculate selection coefficients for all genomes once to avoid redundant calculations when writing attributes
        let (mut selection_weights, logsumexp_value) = self.log_sum_exp();

        selection_weights = selection_weights
            .into_iter()
            .map(|x| (x - logsumexp_value).exp()) // exp(log(w) - logsumexp)
            .collect();

        let sum_weights: f64 = selection_weights.iter().sum();
        selection_weights = selection_weights
            .iter()
            .map(|&w| {
                if w != std::f64::NEG_INFINITY {
                    w / sum_weights
                } else {
                    0.0
                }
            })
            .collect();

        for (genome_index, genome) in self.pop.iter().enumerate() {
            let prefix = if root_genome { "root".to_string() } else { format!("pop_{}_gen_{}_genome_{}", self.id, self.generation, genome_index) };
            let genome_output_path = Self::genome_output_path(output_path, &prefix)?;
            let file = File::create(&genome_output_path)?;
            let mut writer = BufWriter::new(file);

            writeln!(writer, "##gff-version 3")?;

            let genome_selection_coefficient = selection_weights[genome_index].ln();
            let mut contig_offsets: HashMap<usize, usize> = HashMap::new();

            // Sort element indices by (contig_id, seq position) so GFF is grouped
            // by contig and entries within each contig appear in physical order.
            let mut sorted_seq_indices: Vec<usize> = (0..genome.seq.len()).collect();
            sorted_seq_indices.sort_by_key(|&i| (genome.seq[i].contig_id, i));

            for seq_idx in sorted_seq_indices {
                let element = &genome.seq[seq_idx];
                let offset = contig_offsets.entry(element.contig_id).or_insert(0);
                let start_0 = *offset;
                let end_0 = start_0 + element.seq.len();
                *offset = end_0;

                if start_0 >= end_0 {
                    continue;
                }

                // calculate element selection coefficient
                let log_element_selection_coefficient =
                    element.element_selection_coefficient(&genome.identifier);
                let element_selection_coefficient =
                    if log_element_selection_coefficient == std::f64::NEG_INFINITY {
                        0.0
                    } else {
                        ((log_element_selection_coefficient - logsumexp_value).exp() / sum_weights).ln() // exp(log(w) - logsumexp)
                    };

                let seq_id = format!("contig_{}", element.contig_id + 1);
                let start_1based = start_0 + 1;
                let end_1based = end_0;
                let strand = if element.strand { "+" } else { "-" };

                let attributes = format!(
                    "genome_id={};element_id={};feature_type={};feature_id={};contig_id={};genome_identifier={};parent={};multiplier={:.6};sequence_length={};log_genome_selection_coefficient={:.6};log_element_selection_coefficient={:.6}",
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
                );

                writeln!(
                    writer,
                    "{}\tPansimNuc\t{}\t{}\t{}\t.\t{}\t.\t{}",
                    seq_id, element.feature_type, start_1based, end_1based, strand, attributes
                )?;
            }

            writer.flush()?;

            if root_genome {
                break; // only write root genome if specified
            }
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

    fn default_structural_dists() -> Vec<Vec<MutationDistribution>> {
        let mut structural_dists = Vec::new();
        for _ in 0..5 {
            structural_dists.push(vec![
                MutationDistribution::new_poisson(0.1).expect("Failed to create poisson distribution for duplication"),
                MutationDistribution::new_poisson(0.1).expect("Failed to create poisson distribution for deletions"),
                MutationDistribution::new_poisson(0.1).expect("Failed to create poisson distribution for inversions"),
            ]);
        }

        structural_dists
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
        let exon_dist = MutationDistribution::new_double_exp(0.5, 2.0, 0.3)
            .expect("Failed to create double exponential distribution for exon features");
        let intron_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create uniform distribution for intron features");
        let intergenic_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("Failed to create uniform distribution for intergenic features");
        let site_mutation_dists = vec![exon_dist, intron_dist, intergenic_dist];

        let exon_mu = 1.0;
        let intron_mu = 1.0;
        let intergenic_mu = 1.0;
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];

        let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");

        let multiplier_dist = MutationDistribution::new_uniform(0.5, 1.5)
            .expect("Failed to create uniform distribution for TE multipliers");
        let multiplier_dists = vec![multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone()];

        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let pop = Population::new(
            root,
            n_genomes,
            site_mutation_dists,
            &site_mutation_mus,
            &site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            10, // max_multiplier_dist
            multiplier_dists,
            10,
            &mut rng,
            true, // verbose
            &vec!["chr1".to_string()], // contig_name_to_id
            &vec![("chr1".to_string(), 150, 350)], // tracking_regions
            true, // augment_tracking
            0.01, // genome_size_penalty_per_bp
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

        let exon_mu = 1.0;
        let intron_mu = 1.0;
        let intergenic_mu = 1.0;
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];
        let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];

          let multiplier_dist = MutationDistribution::new_uniform(0.5, 1.5)
            .expect("Failed to create uniform distribution for TE multipliers");
        let multiplier_dists = vec![multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone()];


        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let pop = Population::new(
            root,
            1,
            site_mutation_dists,
            &site_mutation_mus,
            &site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            10, // max_multiplier_dist
            multiplier_dists,
            10,
            &mut rng,
            true, // verbose
            &vec!["chr1".to_string()], // contig_name_to_id
            &vec![("chr1".to_string(), 150, 350)], // tracking
            true, // augment_tracking
            0.01, // genome_size_penalty_per_bp
        );

        let temp_path = std::env::temp_dir().join(format!(
            "pansimnuc_pop_{}.fasta",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("system time before UNIX_EPOCH")
                .as_nanos()
        ));
        let output_path = temp_path.to_string_lossy().into_owned();
        let genome_output_path = Population::genome_output_path(&output_path, "pop_0_gen_0_genome_0")
            .expect("failed to construct per-genome output path");

        pop.write_fasta(&output_path, false)
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

        let exon_mu = 1.0;
        let intron_mu = 1.0;
        let intergenic_mu = 1.0;
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];

        let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create recombination probability distribution");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create recombination length distribution");
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];

        let multiplier_dist = MutationDistribution::new_uniform(0.5, 1.5)
            .expect("Failed to create uniform distribution for TE multipliers");
        let multiplier_dists = vec![multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone()];


        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let pop = Population::new(
            root,
            1,
            site_mutation_dists,
            &site_mutation_mus,
            &site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            10, // max_multiplier_dist
            multiplier_dists,
            10,
            &mut rng,
            true, // verbose
            &vec!["chr1".to_string()], // contig_name_to_id
            &vec![("chr1".to_string(), 150, 350)], // tracking
            true, // augment_tracking
            0.01, // genome_size_penalty_per_bp

        );

        let temp_path = std::env::temp_dir().join(format!(
            "pansimnuc_pop_{}.gff3",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("system time before UNIX_EPOCH")
                .as_nanos()
        ));
        let output_path = temp_path.to_string_lossy().into_owned();
        let genome_output_path = Population::genome_output_path(&output_path, "pop_0_gen_0_genome_0")
            .expect("failed to construct per-genome output path");

        pop.write_gff(&output_path, false)
            .expect("failed to write test GFF file");

        let mut content = String::new();
        fs::File::open(&genome_output_path)
            .expect("failed to open test GFF file")
            .read_to_string(&mut content)
            .expect("failed to read test GFF file");

        assert!(content.contains("##gff-version 3"));
        assert!(content.contains("\tintergenic\t"));
        assert!(content.contains("feature_id=0"));
        assert!(content.contains("sequence_length=4"));
        assert!(content.contains("genome_selection_coefficient="));
        assert!(content.contains("element_selection_coefficient="));

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
        let force_mutation_dist = 1000.0;
        let intron_mu_dist = 1.0;
        let intergenic_mu_dist = 1.0;
        let site_mutation_mus = vec![force_mutation_dist, intron_mu_dist, intergenic_mu_dist];
        let indel_mus = vec![1e-12, 1e-12, 1e-12]; // effectively no indels for this test
        let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];
        let multiplier_dist = MutationDistribution::new_uniform(0.5, 1.5)
            .expect("Failed to create uniform distribution for TE multipliers");
        let multiplier_dists = vec![multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone()];


        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let mut pop = Population::new(
            root,
            1,
            site_mutation_dists,
            &site_mutation_mus,
            &indel_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            10, // max_multiplier_dist
            multiplier_dists,
            10,
            &mut rng,
            true, // verbose
            &vec!["chr1".to_string()], // contig_name_to_id
            &vec![("chr1".to_string(), 150, 350)], // tracking
            true, // augment_tracking
            0.01, // genome_size_penalty_per_bp
        );

        let original_seq = pop.pop[0].seq[0].seq.clone();

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

        let exon_mu = 1.0;
        let intron_mu = 1.0;
        let intergenic_mu = 1.0;
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("Failed to create uniform distribution for recombination");
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];
        let multiplier_dist = MutationDistribution::new_uniform(0.5, 1.5)
            .expect("Failed to create uniform distribution for TE multipliers");
        let multiplier_dists = vec![multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone()];


        let mut pop = Population::new(
            root,
            3,
            site_mutation_dists,
            &site_mutation_mus,
            &site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            10, // max_multiplier_dist
            multiplier_dists,
            10,
            &mut rng,
            true, // verbose
            &vec!["chr1".to_string()], // contig_name_to_id
            &vec![("chr1".to_string(), 150, 350)], // tracking
            true, // augment_tracking
            0.01, // genome_size_penalty_per_bp
        );

        let original_identifiers: Vec<String> = pop
            .pop
            .iter()
            .map(|genome| genome.identifier.clone())
            .collect();
        let original_parents: Vec<String> =
            pop.pop.iter().map(|genome| genome.parent.clone()).collect();
        let original_sequences: Vec<Vec<Vec<u8>>> = pop
            .pop
            .iter()
            .map(|genome| {
                genome
                    .seq
                    .iter()
                    .map(|element| element.seq.clone())
                    .collect()
            })
            .collect();

        pop.next_generation(vec![2, 0, 1]);

        assert_eq!(pop.generation, 1);
        assert_eq!(pop.pop.len(), 3);

        for (new_index, genome) in pop.pop.iter().enumerate() {
            let selected_index = [2usize, 0, 1][new_index];
            assert_eq!(
                genome.identifier,
                format!("1-{}", original_identifiers[selected_index])
            );
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

        let exon_mu = 1.0;
        let intron_mu = 1.0;
        let intergenic_mu = 1.0;
        let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];

        let recombination_prob_dist = MutationDistribution::new_poisson(1.0)
            .expect("failed to create recombination count distribution");
        let recombination_len_dist = MutationDistribution::new_poisson(1.0)
            .expect("failed to create recombination length distribution");
        let recombination_dists = vec![recombination_prob_dist, recombination_len_dist];
        let multiplier_dist = MutationDistribution::new_uniform(0.5, 1.5)
            .expect("Failed to create uniform distribution for TE multipliers");
        let multiplier_dists = vec![multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone(), multiplier_dist.clone()];


        let mut rng: StdRng = StdRng::seed_from_u64(123);
        Population::new(
            root,
            1,
            site_mutation_dists,
            &site_mutation_mus,
            &site_mutation_mus,
            recombination_dists,
            1.0,
            default_structural_dists(),
            10, // max_multiplier_dist
            multiplier_dists,
            10,
            &mut rng,
            true, // verbose
            &vec!["chr1".to_string()], // contig_name_to_id
            &vec![("chr1".to_string(), 150, 350)], // tracking
            true, // augment_tracking
            0.01, // genome_size_penalty_per_bp
        )
    }

    fn genome_from_seq(seq: Vec<NucElement>) -> Genome {
        Genome {
            identifier: "test".to_string(),
            genome_id: 0,
            contig_starts: vec![0],
            parent: "test-parent".to_string(),
            seq,
            seq_length: 0,
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

    #[test]
    fn test_check_feature_order_identifies_te_2_upstream_and_downstream() {
        let mut pop = make_check_feature_order_population();
        let mut seq = pop.pop[0].seq.clone();

        pop.max_multiplier_dist = 10; // ensure the TEs we add are within the max multiplier distance

        let mut upstream_te = seq[0].clone();
        upstream_te.feature_type = "TE-CUT".to_string();
        upstream_te.feature_id = 0;
        upstream_te.multiplier = 2.0;
        upstream_te.element_id = 20_000;
        seq.insert(0, upstream_te);

        let mut downstream_te = seq[0].clone();
        downstream_te.feature_type = "TE-COPY".to_string();
        downstream_te.feature_id = 0;
        downstream_te.multiplier = 3.5;
        downstream_te.element_id = 20_001;
        seq.push(downstream_te);

        println!(
            "Sequence: {:?}",
            seq.iter()
                .map(|e| format!("{}-{}", e.feature_type, e.feature_id))
                .collect::<Vec<String>>()
        );

        let genome = genome_from_seq(seq);
        let (broken, multiplier) = check_feature_one_intron(&pop, &genome);
        assert!(!broken);
        assert!(multiplier == 3.5);
    }

    #[test]
    fn test_check_feature_order_non_identifies_te_2_upstream_and_downstream() {
        let mut pop = make_check_feature_order_population();
        let mut seq = pop.pop[0].seq.clone();

        pop.max_multiplier_dist = 2; // ensure the TEs we add are not within the max multiplier distance

        let mut upstream_te = seq[0].clone();
        upstream_te.feature_type = "TE-CUT".to_string();
        upstream_te.feature_id = 0;
        upstream_te.multiplier = 2.0;
        upstream_te.element_id = 20_000;
        seq.insert(0, upstream_te);

        let mut downstream_te = seq[0].clone();
        downstream_te.feature_type = "TE-COPY".to_string();
        downstream_te.feature_id = 0;
        downstream_te.multiplier = 3.5;
        downstream_te.element_id = 20_001;
        seq.push(downstream_te);

        println!(
            "Sequence: {:?}",
            seq.iter()
                .map(|e| format!("{}-{}", e.feature_type, e.feature_id))
                .collect::<Vec<String>>()
        );

        let genome = genome_from_seq(seq);
        let (broken, multiplier) = check_feature_one_intron(&pop, &genome);
        assert!(!broken);
        assert!(multiplier == 1.0);
    }

    fn make_test_element_with_coefficients(
        feature_id: usize,
        seq: Vec<u8>,
        coefficients: &[(usize, u8, f64)],
    ) -> NucElement {
        let mut rng: StdRng = StdRng::seed_from_u64(99);
        let seed_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create uniform distribution for seeded mutation map");
        let mut mutation_map = MutationMap::new(0, 0, &seq, &seed_dist, &mut rng);

        for (site, allele, coeff) in coefficients {
            mutation_map.set_for_test(*allele, *site, *coeff);
        }

        NucElement {
            contig_id: 0,
            element_id: 0,
            feature_id: feature_id,
            feature_type: "exon".to_string(),
            multiplier: 1.0,
            seq: seq.clone(),
            mutation_map,
            strand: true,
            original_length: seq.len(),
            frameshift: false,
            tracked: false,
        }
    }

    #[test]
    fn test_element_selection_coefficient_returns_negative_infinity_when_any_site_is_minus_one() {
        let coefficients = vec![(0, 1, 0.2), (1, 2, -1.0), (2, 4, 0.5)];
        let element = make_test_element_with_coefficients(
            0,
            vec![1, 2, 4],
            &coefficients,
        );

        let log_sum = element.element_selection_coefficient("test-genome");
        assert_eq!(log_sum, std::f64::NEG_INFINITY);
    }

    #[test]
    fn test_element_selection_coefficient_is_deterministic_for_predefined_coefficients() {
        let coefficients = vec![(0, 1, 0.2), (1, 2, 0.5), (2, 4, 0.1)];
        let element = make_test_element_with_coefficients(
            1,
            vec![1, 2, 4],
            &coefficients,
        );

        let val1 = 1.0 + coefficients[0].2; // 1.0 + 0.2
        let val2 = 1.0 + coefficients[1].2; // 1.0 + 0.5
        let val3 = 1.0 + coefficients[2].2; // 1.0 + 0.1 

        let expected = val1.ln() + val2.ln() + val3.ln();

        let log_sum = element.element_selection_coefficient("test-genome");

        assert!(log_sum == expected);
    }

    #[test]
    fn test_genome_selection_coefficient_neg_infinity_on_lethal_mutation() {
        let pop = make_check_feature_order_population();

        let neutral = {
            let mut e = make_test_element_with_coefficients(1,vec![1u8], &[(0, 1u8, 0.0)]);
            e.feature_type = "intergenic".to_string();
            e
        };
        let lethal = {
            let mut e = make_test_element_with_coefficients(2,vec![2u8], &[(0, 2u8, -1.0)]);
            e.feature_type = "intergenic".to_string();
            e
        };

        let genome = genome_from_seq(vec![neutral, lethal]);
        assert_eq!(
            pop.genome_selection_coefficient(&genome),
            std::f64::NEG_INFINITY
        );
    }

    #[test]
    fn test_genome_selection_coefficient_correct_sum() {
        let mut pop = make_check_feature_order_population();
        // Override feature_map so feature_id=1 maps to a single exon (element_id=1)
        pop.feature_map.insert(1, vec![1]);

        let e1 = {
            let mut e = make_test_element_with_coefficients(
                0,
                vec![1u8, 2u8],
                &[(0, 1u8, 0.5), (1, 2u8, 0.3)],
            );
            e.feature_type = "intergenic".to_string();
            e.element_id = 0;
            e
        };
        let e2 = {
            let mut e = make_test_element_with_coefficients(1,vec![4u8], &[(0, 4u8, 0.2)]);
            e.feature_type = "exon".to_string();
            e.element_id = 1;
            e
        };

        let e1_val1: f64 = 1.0 + 0.5;
        let e1_val2: f64 = 1.0 + 0.3;
        let e1_val = e1_val1.ln() + e1_val2.ln(); // ln((1.0 + 0.5) * (1.0 + 0.3))
        let e2_val1: f64 = 1.0 + 0.2;
        let e2_val = e2_val1.ln(); // ln(1.2)

        let genome = genome_from_seq(vec![e1, e2]);
        let expected = e1_val + e2_val;
        let actual = pop.genome_selection_coefficient(&genome);
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_genome_selection_coefficient_applies_te_multiplier() {
        let mut pop = make_check_feature_order_population();

        let e1 = {
            let mut e = make_test_element_with_coefficients(
                0,
                vec![1u8, 2u8],
                &[(0, 1u8, 0.5), (1, 2u8, 0.3)],
            );
            e.feature_type = "intergenic".to_string();
            e.element_id = 0;
            e
        };
        let e2 = {
            let mut e = make_test_element_with_coefficients(1, vec![4u8], &[(0, 4u8, 0.2)]);
            e.feature_type = "exon".to_string();
            e.element_id = 1;
            e
        };

        // Override feature_map: feature_id=1 is a single-exon gene with element_id=1
        pop.feature_map.insert(1, vec![1]);

        let e1_val1: f64 = 1.0 + 0.5;
        let e1_val2: f64 = 1.0 + 0.3;
        let e1_val = e1_val1.ln() + e1_val2.ln(); // ln((1.0 + 0.5) * (1.0 + 0.3))
        let e2_val1: f64 = 1.0 + 0.2;
        let e2_val = e2_val1.ln(); // ln(1.2)

        // initial check of equivalence
        let genome = genome_from_seq(vec![e1.clone(), e2.clone()]);
        let expected_pre = e1_val + e2_val;
        let actual_pre = pop.genome_selection_coefficient(&genome);
        assert!(actual_pre == expected_pre);

        // Insert a neutral TE-CUT (coeff=0, so te_coeff=0) directly between e1 and e2.
        // Being 1 position upstream of the exon, its multiplier=2.0 scales e2's contribution.
        let te = {
            let mut e = make_test_element_with_coefficients(0, vec![1u8], &[(0, 1u8, 1.0)]);
            e.feature_type = "TE-CUT".to_string();
            e.element_id = 99_999;
            e.multiplier = 2.0;
            e
        };
        let seq_with_te = vec![e1.clone(), te, e2.clone()];
        // feature_map is unchanged: feature_id=1 still maps to [element_id=1]
        let genome_with_te = genome_from_seq(seq_with_te);

        let te_val1: f64 = 1.0 + 1.0; // TE-CUT with coeff=1.0
        let te_val = te_val1.ln(); // ln(2.0)

        let coeff_with_te = pop.genome_selection_coefficient(&genome_with_te);

        // te_coeff = ln(1.0) = 0.0; e2 is scaled by 2.0
        let expected_with_te = e1_val + 0.0 + e2_val * 2.0 + te_val;
        assert_eq!(
            coeff_with_te as f32,
            expected_with_te as f32,
            "Expected e2 contribution to be scaled by TE multiplier"
        );
    }

    // -----------------------------------------------------------------------
    // Helper shared by the tracking tests below
    // -----------------------------------------------------------------------
    fn make_two_exon_population_with_tracking(
        structural_mu_dists: Vec<Vec<MutationDistribution>>,
        tracking_region: (String, usize, usize),
        seq_len: usize,
    ) -> Population {
        let root = vec![vec![
            FeaturePos {
                contig_id: 0,
                feature_id: 1,
                feature_type: "exon".to_string(),
                start: 0,
                end: seq_len,
                strand: true,
                seq: vec![1u8; seq_len],
            },
            FeaturePos {
                contig_id: 0,
                feature_id: 2,
                feature_type: "exon".to_string(),
                start: seq_len,
                end: seq_len * 2,
                strand: true,
                seq: vec![1u8; seq_len],
            },
        ]];

        // 6 selection dists: indices 0-4 for feature types, 5 for tracked elements
        let selection_dists: Vec<MutationDistribution> = (0..6)
            .map(|_| MutationDistribution::new_uniform(0.0, 1.0).unwrap())
            .collect();
        // 6 mu/indel rates (index 5 accessed for tracked elements)
        let mu_dist_vals = vec![1.0f64; 6];
        let indel_dist_vals = vec![1e-6f64; 6];
        let recombination_dists = vec![
            MutationDistribution::new_poisson(0.01).unwrap(),
            MutationDistribution::new_poisson(1.0).unwrap(),
        ];
        // 6 multiplier dists (index 5 for tracked elements)
        let multiplier_dists: Vec<MutationDistribution> = (0..6)
            .map(|_| MutationDistribution::new_uniform(0.5, 1.5).unwrap())
            .collect();

        let mut rng = StdRng::seed_from_u64(0);
        Population::new(
            root,
            1,
            selection_dists,
            &mu_dist_vals,
            &indel_dist_vals,
            recombination_dists,
            1.0,
            structural_mu_dists,
            10,
            multiplier_dists,
            1,
            &mut rng,
            false,
            &vec!["chr1".to_string()],
            &vec![tracking_region],
            true,
            0.01
        )
    }

    fn zero_structural_dists() -> Vec<Vec<MutationDistribution>> {
        let zero = MutationDistribution::new_uniform(0.0, 0.5).unwrap();
        (0..6)
            .map(|_| vec![zero.clone(), zero.clone(), zero.clone()])
            .collect()
    }

    #[test]
    fn test_population_new_tracked_exon_uses_dist_id_5_untracked_keeps_type_dist_id() {
        // Two exon features at [0, 100) and [100, 200).
        // Tracking region [0, 99] covers only the first exon.
        // After Population::new the tracked exon should have selection_dist_id=5
        // and mu_dist_id=5; the untracked exon should keep id=0 (exon default).
        let pop = make_two_exon_population_with_tracking(
            zero_structural_dists(),
            ("chr1".to_string(), 0, 99),
            100,
        );

        let exon0 = &pop.pop[0].seq[0]; // first exon — inside tracking region
        let exon1 = &pop.pop[0].seq[1]; // second exon — outside tracking region

        assert!(exon0.tracked, "first exon inside region should be tracked");
        assert_eq!(exon0.mutation_map.selection_dist_id, 5,
            "tracked exon should sample selection coefficients from tracked dist (id=5)");
        assert_eq!(exon0.mutation_map.mu_dist_id, 5,
            "tracked exon should use tracked mutation rate dist (id=5)");

        assert!(!exon1.tracked, "second exon outside region should not be tracked");
        assert_eq!(exon1.mutation_map.selection_dist_id, 0,
            "untracked exon should keep exon selection dist (id=0)");
        assert_eq!(exon1.mutation_map.mu_dist_id, 0,
            "untracked exon should keep exon mutation rate dist (id=0)");
    }

    #[test]
    fn test_population_new_tracked_exon_uses_tracked_structural_dist() {
        // Two exon features at [0, 10) and [10, 20).
        // Tracking region [0, 9] covers only the first exon so it is tracked.
        // structural_mu_dists[0] (exon path) guarantees 10 duplications;
        // structural_mu_dists[5] (tracked path) guarantees 0 duplications.
        // After structural_intra_genome the genome should contain exactly 1 tracked
        // copy and 11 untracked copies (original + 10 dups).
        let zero = MutationDistribution::new_uniform(0.0, 0.5).unwrap(); // → 0 when cast to usize
        let ten_dups = MutationDistribution::new_uniform(10.0, 10.5).unwrap(); // → 10 when cast to usize
        let structural_mu_dists = vec![
            vec![ten_dups, zero.clone(), zero.clone()], // 0: exon — 10 dups, 0 dels, 0 inv
            vec![zero.clone(), zero.clone(), zero.clone()], // 1: intron
            vec![zero.clone(), zero.clone(), zero.clone()], // 2: intergenic
            vec![zero.clone(), zero.clone(), zero.clone()], // 3: TE-CUT
            vec![zero.clone(), zero.clone(), zero.clone()], // 4: TE-COPY
            vec![zero.clone(), zero.clone(), zero.clone()], // 5: tracked — 0 dups
        ];

        let mut pop = make_two_exon_population_with_tracking(
            structural_mu_dists,
            ("chr1".to_string(), 0, 9),
            10,
        );

        pop.structural_intra_genome();

        let tracked_count = pop.pop[0].seq.iter().filter(|e| e.tracked).count();
        let untracked_count = pop.pop[0].seq.iter().filter(|e| !e.tracked).count();

        assert_eq!(tracked_count, 1,
            "tracked exon should not be duplicated (structural_mu_dists[5] has 0 dups)");
        assert_eq!(untracked_count, 11,
            "untracked exon should be duplicated 10 times (structural_mu_dists[0] has 10 dups)");
    }

    #[test]
    fn test_genome_size_penalty_impacts_selection_weights() {
        let mut root = Vec::new();
        root.push(vec![FeaturePos {
            contig_id: 0,
            feature_id: 0,
            feature_type: "intergenic".to_string(),
            start: 0,
            end: 4,
            strand: true,
            seq: vec![1, 2, 4, 8], // 4 bp
        }]);

        let selection_dist = MutationDistribution::new_uniform(0.0, 1.0)
            .expect("failed to create selection dist");
        let site_mutation_dists = vec![
            selection_dist.clone(), // exon
            selection_dist.clone(), // intron
            selection_dist,         // intergenic
        ];
        let mu_vals = vec![1e-12, 1e-12, 1e-12]; // effectively no mutations

        let recomb_prob = MutationDistribution::new_poisson(1.0)
            .expect("failed to create recombination prob dist");
        let recomb_len = MutationDistribution::new_poisson(1.0)
            .expect("failed to create recombination len dist");
        let recombination_dists = vec![recomb_prob, recomb_len];

        let multiplier_dist = MutationDistribution::new_uniform(0.5, 1.5)
            .expect("failed to create multiplier dist");
        let multiplier_dists = vec![
            multiplier_dist.clone(),
            multiplier_dist.clone(),
            multiplier_dist.clone(),
            multiplier_dist.clone(),
            multiplier_dist,
        ];

        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let mut pop = Population::new(
            root,
            2,
            site_mutation_dists,
            &mu_vals,
            &mu_vals,
            recombination_dists,
            1.0,
            default_structural_dists(),
            10,
            multiplier_dists,
            10,
            &mut rng,
            false,
            &vec!["chr1".to_string()],
            &vec![],  // no tracking
            false,
            1.0, // genome_size_penalty_per_bp
        );

        assert_eq!(pop.optimal_genome_size, 4);

        // Give genome 0 a large deviation and genome 1 a small deviation.
        // The penalty multiplies each weight by (penalty_per_bp * |deviation|), so
        // genome 0 (deviation=1000) will have ~10x the weight of genome 1 (deviation=100).
        pop.pop[0].seq_length = pop.optimal_genome_size + 1000;
        pop.pop[1].seq_length = pop.optimal_genome_size + 100;

        let mut thread_rng = rand::thread_rng();
        let n_rounds = 500;
        let mut counts = [0usize; 2];
        for _ in 0..n_rounds {
            for idx in pop.sample_individuals(&mut thread_rng) {
                counts[idx] += 1;
            }
        }

        // Genome 0 has a 10x larger deviation so it should be sampled far more often.
        assert!(
            counts[0] > counts[1],
            "genome with larger size deviation should have higher selection weight; counts = {:?}",
            counts
        );
    }
}
