use crate::gff::FeaturePos;
use crate::mutation::MutationMap;
use crate::mutation::Distribution;
use std::collections::HashMap;
use statrs::distribution::Poisson;

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
    pub selection_dists: Vec<Distribution>,
    pub mu_dists: Vec<Distribution>,
}

impl Population {
    pub fn new(
        root: HashMap<String, Vec<FeaturePos>>,
        n_genomes: usize,
        selection_dists: Vec<Distribution>,
        mu_dists: Vec<Distribution>
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

    pub fn mutate (&mut self) {
        for genome in &mut self.pop {
            for element in &mut genome.seq {
                element.mutation_map.mutate(&self.core_vec, &mut element.seq, &self.selection_dists[element.mutation_map.selection_dist_id], &self.mu_dists[element.mutation_map.mu_dist_id],);


                // apply mutations to element.seq based on element.mutation_map and self.distributions
            }
        }
        // placeholder for mutation function, which will apply mutations to each genome in the population based on the mutation map of each NucElement and the specified distributions
    }

    pub fn next_generation (&mut self) {

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
        let exon_dist = Distribution::new_double_exp(0.5, 2.0, 0.3).expect("Failed to create double exponential distribution for exon features");
        let intron_dist = Distribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intron features");
        let intergenic_dist = Distribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intergenic features");
        let site_mutation_dists = vec![exon_dist, intron_dist, intergenic_dist];

        let exon_mu = Distribution::new_uniform(0.0, 1.0).expect("Failed to create double exponential distribution for exon features");
        let intron_mu = Distribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intron features");
        let intergenic_mu = Distribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intergenic features");
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
