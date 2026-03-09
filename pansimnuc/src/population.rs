use crate::gff::FeaturePos;
use crate::mutation::MutationMap;
use crate::mutation::Distribution as MutationDistribution;
use std::collections::HashMap;
use rayon::prelude::*;
use logsumexp::LogSumExp;
use rand::distributions::{Distribution as RandDistribution, WeightedIndex};
use rand::rngs::StdRng;

#[derive(Clone)]
pub struct NucElement {
    pub seqname: String,
    pub feature_id: usize,
    pub feature_type: String,
    pub seq: Vec<u8>,
    pub mutation_map: MutationMap
}

pub struct Genome {
    pub identifier: String,
    pub parent: String,
    pub seq: Vec<NucElement>
}

pub struct Population{
    pub generation: usize,
    pub pop: Vec<Genome>,
    pub core_vec: Vec<Vec<u8>>,
    pub selection_dists: Vec<MutationDistribution>,
    pub mu_dists: Vec<MutationDistribution>,
}

impl Population {
    pub fn new(
        root: HashMap<String, Vec<FeaturePos>>,
        n_genomes: usize,
        selection_dists: Vec<MutationDistribution>,
        mu_dists: Vec<MutationDistribution>
    ) -> Self {
        let mut population: Vec<Genome> = Vec::new();

        for i in 0..n_genomes {
            let mut genome: Vec<NucElement> = Vec::new();
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
                        feature_id: feature.feature_id,
                        feature_type: feature.feature_type.clone(),
                        seq: feature.seq.clone(),
                        mutation_map: MutationMap::new(selection_dist_id, mu_dist_id),
                    });
                }
            }
            population.push(Genome {
                identifier: format!("{}", i),
                parent: "root".to_string(),
                seq: genome,
            });
        }
        
        let core_vec: Vec<Vec<u8>> =
            vec![vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4]];

        Self {
            generation: 0,
            pop: population,
            core_vec,
            selection_dists,
            mu_dists
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

    // sample individuals using logsumexp normalisation to prevent underflow/overflow issues with very small/large weights
    pub fn sample_individuals (&self, rng: &mut StdRng) -> Vec<usize> {
        let mut selection_weights: Vec<f64> = vec![1.0; self.pop.len()];
        
        selection_weights = self.pop
            .par_iter()
            .map(|genome| {
                let mut log_sum = 0.0;
                for element in &genome.seq {
                    for (site, allele) in element.seq.iter().enumerate() {
                        if let Some(coeff) = element.mutation_map.get(*allele, site) {
                            log_sum += coeff;
                        } else {
                            panic!(
                                "Failed to generate selection coefficient for genome {} allele {} at site {}",
                                genome.identifier, allele, site
                            );
                        }
                    }
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

        sampled_indices
    }

    pub fn next_generation (&mut self, sampled_indices: Vec<usize>) {
        let mut new_pop: Vec<Genome> = Vec::new();

        for &selected_index in &sampled_indices {
            let selected_genome = &self.pop[selected_index];
            new_pop.push(Genome {
                identifier: format!("{}-{}", self.generation + 1, selected_genome.identifier),
                parent: selected_genome.identifier.clone(),
                seq: selected_genome.seq.clone(),
            });
        }
        self.pop = new_pop;
        self.generation += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gff::FeaturePos;

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


        let pop = Population::new(root, n_genomes, site_mutation_dists, site_mutation_mus);

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
}
