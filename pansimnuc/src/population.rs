use crate::gff::FeaturePos;
use crate::mutation::MutationMap;
use crate::mutation::Distribution;
use std::collections::HashMap;

pub struct NucElement<'a> {
    pub seqname: String,
    pub feature_id: usize,
    pub feature_type: String,
    pub seq: Vec<u8>,
    pub mutation_map: MutationMap<'a>
}

pub struct Genome<'a> {
    pub identifier: String,
    pub parent: String,
    pub seq: Vec<NucElement<'a>>
}

pub struct Population<'a>{
    pub generation: usize,
    pub pop: Vec<Genome<'a>>,
    pub core_vec: Vec<Vec<u8>>
}

impl<'a> Population<'a> {
    pub fn new(
        root: HashMap<String, Vec<FeaturePos>>,
        n_genomes: usize,
        // TODO add way of specifying number of mutation distribution types, which can then be added by reference to mutation_map
    ) -> Self {
        let mut population: Vec<Genome> = Vec::new();

        // placeholder for mutation map, which will be added to each NucElement in the genome
        let exon_dist = Distribution::new_double_exp(0.0, 1.0, 0.5).expect("Failed to create double exponential distribution for exon features");
        let intron_dist = Distribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intron features");
        let intergenic_dist = Distribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intergenic features");

        for i in 0..n_genomes {
            let mut genome: Vec<NucElement> = Vec::new();
            for (seqname, features) in &root {
                for feature in features {
                    let distribution = match feature.feature_type.as_str() {
                        "exon" => exon_dist,
                        "intron" => intron_dist,
                        "intergenic" => intergenic_dist,
                        _ => panic!("Unknown feature type: {}", feature.feature_type),
                    };

                    genome.push(NucElement {
                        seqname: feature.seqname.clone(),
                        feature_id: feature.feature_id,
                        feature_type: feature.feature_type.clone(),
                        mutation_map: MutationMap::new(distribution, &mut feature.seq),
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
            core_vec
        }
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
                feature_type: "gene".to_string(),
                start: 100,
                end: 200,
                strand: true,
                seq: vec![0, 1, 2, 3], // ACGT
            },
            FeaturePos {
                seqname: "chr1".to_string(),
                feature_id: 1,
                feature_type: "gene".to_string(),
                start: 300,
                end: 400,
                strand: false,
                seq: vec![3, 2, 1, 0], // TGCA
            },
        ];
        root.insert("chr1".to_string(), features);

        let n_genomes = 3;
        let pop = Population::new(root, n_genomes);

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
            assert_eq!(genome.seq[0].feature_type, "gene");

            // Check second feature
            assert_eq!(genome.seq[1].seqname, "chr1");
            assert_eq!(genome.seq[1].feature_id, 1);
            assert_eq!(genome.seq[1].feature_type, "gene");
        }
    }
}
