use crate::gff::FeaturePos;
use std::collections::HashMap;
pub struct NucElement {
    pub seqname: String,
    pub feature_id: usize,
    pub feature_type: String,
    seq: Vec<u8>,
}

// impl NucElement {
//     fn new (
//         seqname: String,
//         feature_id: usize,
//         feature_type: String,
//         seq: Vec<u8>
//     ) -> Self {
//         Self {
//             seqname,
//             feature_id,
//             feature_type,
//             seq
//         }
//     }
// }

pub struct Genome {
    pub identifier: String,
    pub parent: String,
    pub seq: Vec<NucElement>
}

// impl Genome {
//     fn new (
//         identifier: String,
//         parent: String,
//         seq: Vec<NucElement>
//     ) -> Self {
//         Self {
//             identifier,
//             parent,
//             seq
//         }
//     }
// }

pub struct Population{
    pub generation: usize,
    pub pop: Vec<Genome>
}

impl Population {
    pub fn new(
        root: HashMap<String, Vec<FeaturePos>>,
        n_genomes: usize,
    ) -> Self {
        let mut population: Vec<Genome> = Vec::new();

        for i in 0..n_genomes {
            let mut genome: Vec<NucElement> = Vec::new();
            for (seqname, features) in &root {
                for feature in features {
                    genome.push(NucElement {
                        seqname: feature.seqname.clone(),
                        feature_id: feature.feature_id,
                        feature_type: feature.feature_type.clone(),
                        seq: feature.seq.clone()
                    });
                }
            }
            population.push(Genome {
                identifier: format!("{}", i),
                parent: "root".to_string(),
                seq: genome,
            });
        }

        Self {
            generation: 0,
            pop: population
        }
    }
}
