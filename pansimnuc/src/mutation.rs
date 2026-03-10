// TODO - each site has a selection coefficient drawn at start time
// if a new mutation arises, store the selection coefficient for that allele
// so that repeated mutations go back to the same selection coefficient
// means don't have to store all selection coefficients until mutation to that one occurs
// then selection coefficient for the whole block is based on the product of all
// variants in the block

// also need to come up with way of insertion into/out of block e.g. by TE
// which negates all selection coefficient of gene
// could be beneficial or deleterious with insertion, depending on effect of the 
// given gene in first place.

// recombination is done by block (could do whole blocks) - needs to be homologous, could have shared
// ancestry that is held in parent/id? Or calculate edit distance? Probably will take too long

use rustc_hash::FxHashMap;
use rand::rngs::StdRng;
use rand::{Rng};
use rand::seq::IteratorRandom;
use rand_distr::{Normal, Uniform, Exp, Distribution as RandDist};
use std::fmt;
use statrs::distribution::Poisson;

#[derive(Debug)]
pub enum DistributionError {
    InvalidNormalParameters,
    InvalidUniformParameters,
    InvalidExponentialParameters,
    InvalidDoubleExpParameters,
    InvalidPoissonParameters,
}

impl fmt::Display for DistributionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            DistributionError::InvalidNormalParameters => {
                write!(f, "Invalid Normal distribution: std_dev must be positive")
            }
            DistributionError::InvalidUniformParameters => {
                write!(f, "Invalid Uniform distribution: low must be less than high")
            }
            DistributionError::InvalidExponentialParameters => {
                write!(f, "Invalid Exponential distribution: lambda must be positive")
            }
            DistributionError::InvalidDoubleExpParameters => {
                write!(f, "Invalid DoubleExp distribution: lambdas must be positive and cutoff must be between 0 and 1")
            }
            DistributionError::InvalidPoissonParameters => {
                write!(f, "Invalid Poisson distribution: lambda must be positive")
            }
        }
    }
}

impl std::error::Error for DistributionError {}

pub struct DoubleExponential {
    exp1: Exp<f64>,
    exp2: Exp<f64>,
    cutoff: f64,
    weight: f64,
}

impl DoubleExponential {
    pub fn new(lambda1: f64, lambda2: f64, cutoff: f64, weight: f64) -> Result<Self, DistributionError> {
        if lambda1 <= 0.0 || lambda2 <= 0.0 || cutoff < 0.0 || cutoff > 1.0 {
            return Err(DistributionError::InvalidDoubleExpParameters);
        }
        
        let exp1 = Exp::new(lambda1).map_err(|_| DistributionError::InvalidDoubleExpParameters)?;
        let exp2 = Exp::new(lambda2).map_err(|_| DistributionError::InvalidDoubleExpParameters)?;
        
        Ok(Self { exp1, exp2, cutoff, weight })
    }

    pub fn weight<R: Rng>(&mut self, uniform: &Uniform<f64>, rng: &mut R) {
        self.weight = uniform.sample(rng);
    }

    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        if self.weight <= self.cutoff {
            // Higher selected gene
            self.exp2.sample(rng)
        } else {
            // lower selected gene
            self.exp1.sample(rng)
        }
    }
}

pub enum Distribution {
    Normal(Normal<f64>),
    Uniform(Uniform<f64>),
    Exp(Exp<f64>),
    DoubleExp(DoubleExponential),
    Poisson(Poisson)
}

impl Distribution {
    pub fn new_normal(mean: f64, std_dev: f64) -> Result<Self, DistributionError> {
        Normal::new(mean, std_dev)
            .map(Distribution::Normal)
            .map_err(|_| DistributionError::InvalidNormalParameters)
    }

    pub fn new_uniform(low: f64, high: f64) -> Result<Self, DistributionError> {
        if low >= high {
            return Err(DistributionError::InvalidUniformParameters);
        }
        Ok(Distribution::Uniform(Uniform::new(low, high)))
    }

    pub fn new_exp(lambda: f64) -> Result<Self, DistributionError> {
        Exp::new(lambda)
            .map(Distribution::Exp)
            .map_err(|_| DistributionError::InvalidExponentialParameters)
    }

    pub fn new_double_exp(lambda1: f64, lambda2: f64, cutoff: f64) -> Result<Self, DistributionError> {
        DoubleExponential::new(lambda1, lambda2, cutoff, 0.0)
            .map(Distribution::DoubleExp)
    }

    pub fn new_poisson(lambda: f64) -> Result<Self, DistributionError> {
        Poisson::new(lambda)
            .map(Distribution::Poisson)
            .map_err(|_| DistributionError::InvalidPoissonParameters)
    }

    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        match self {
            Distribution::Normal(d) => d.sample(rng),
            Distribution::Uniform(d) => d.sample(rng),
            Distribution::Exp(d) => d.sample(rng),
            Distribution::DoubleExp(d) => d.sample(rng),
            Distribution::Poisson(d) => d.sample(rng),
        }
    }
}

#[derive(Clone)]
pub struct MutationMap {
    pub selection_dist_id: usize,
    pub mu_dist_id: usize,
    data: [FxHashMap<usize, f64>; 5],
}

impl MutationMap {
    fn allele_to_index(level: u8) -> Option<usize> {
        // Convert one-hot allele code to zero-based index using bit shifting:
        // 1 -> 0, 2 -> 1, 4 -> 2, 8 -> 3.
        // N (16) is represented at index 4.
        if level == 16 {
            return Some(4);
        }

        if level == 0 || (level & (level - 1)) != 0 {
            panic!("Allele code must be one-hot (1, 2, 4, 8); got {}", level);
        }

        let mut idx = 0usize;
        let mut value = level;
        while value > 1 {
            value >>= 1;
            idx += 1;
        }

        if idx >= 4 {
            panic!("Allele code out of range for A/C/G/T map: {}", level);
        }

        Some(idx)
    }

    pub fn new(selection_dist_id: usize, mu_dist_id: usize, seq: &Vec<u8>, selection_dist: &Distribution, rng: &mut StdRng) -> Self {
        let mut data = std::array::from_fn(|_| FxHashMap::default());
        
        for (site, allele) in seq.iter().enumerate() {
            let allele_index = Self::allele_to_index(*allele)
                .expect("Allele code conversion failed while building mutation map");
            data[allele_index].insert(site, selection_dist.sample(rng));
        }

        Self {selection_dist_id, mu_dist_id, data}
    }

    fn insert(&mut self, level: u8, key: usize, value: f64) {
        let allele_index = Self::allele_to_index(level)
            .expect("Cannot insert selection coefficient for invalid allele code");
        self.data[allele_index].insert(key, value);
    }

    pub fn get(&self, level: u8, key: usize) -> Option<&f64> {
        let allele_index = Self::allele_to_index(level)
            .expect("Cannot lookup selection coefficient for invalid allele code");
        self.data[allele_index].get(&key)
    }    

    pub fn mutate (& mut self, core_vec: &Vec<Vec<u8>>, seq: &mut Vec<u8>, selection_dist: &Distribution, mu_dist: &Distribution) {
        // thread-specific random number generator
        let mut thread_rng = rand::thread_rng();

        // sample from Poisson distribution for number of sites to mutate in this isolate
        let n_sites = mu_dist.sample(&mut thread_rng) as usize;
        let seq_len = seq.len();
        let sampled_sites = (0..seq_len).choose_multiple(&mut thread_rng, n_sites);

        // iterate for number of mutations required to reach mutation rate
        for mutant_site in sampled_sites {
            // sample new site to mutate
            let value = seq[mutant_site];

            // N is stored in the map but should not be mutated; skip
            if value == 16 {
                continue;
            }
            
            let allele_index = Self::allele_to_index(value)
                .expect("Sampled non-mutable N allele for mutation");
            let values = &core_vec[allele_index];

            // sample new allele
            let new_allele = values.iter().choose_multiple(&mut thread_rng, 1)[0];
            let mut selection_coefficient: f64 = 0.0;

            // generate new selection coefficient for this mutation if necessary, otherwise retrieve existing one
            if let Some(coeff) = self.get(*new_allele, mutant_site) {
                // value exists
                selection_coefficient = *coeff;
            } else {
                selection_coefficient = selection_dist.sample(&mut thread_rng);
            }

            // set value in place
            seq[mutant_site] = *new_allele;
            self.insert(*new_allele, mutant_site, selection_coefficient);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn test_normal_distribution_creation() {
        let dist = Distribution::new_normal(0.0, 1.0);
        assert!(dist.is_ok());
    }

    #[test]
    fn test_uniform_distribution_creation() {
        let dist = Distribution::new_uniform(-1.0, 1.0);
        assert!(dist.is_ok());
    }

    #[test]
    fn test_uniform_distribution_invalid_range() {
        let dist = Distribution::new_uniform(1.0, 1.0);
        assert!(dist.is_err());
        
        let dist2 = Distribution::new_uniform(2.0, 1.0);
        assert!(dist2.is_err());
    }

    #[test]
    fn test_exp_distribution_creation() {
        let dist = Distribution::new_exp(1.5);
        assert!(dist.is_ok());
    }

    #[test]
    fn test_double_exp_distribution_creation() {
        let dist = Distribution::new_double_exp(0.5, 2.0, 0.3);
        assert!(dist.is_ok());
    }

    #[test]
    fn test_poisson_creation() {
        let dist = Distribution::new_poisson(1.1);
        assert!(dist.is_ok());
    }

    #[test]
    fn test_poisson_invalid_params() {
        let dist = Distribution::new_poisson(-1.1);
        assert!(dist.is_err());
    }

    #[test]
    fn test_double_exp_distribution_invalid_params() {
        // Invalid lambda1
        let dist1 = Distribution::new_double_exp(-0.5, 2.0, 0.3);
        assert!(dist1.is_err());
        
        // Invalid lambda2
        let dist2 = Distribution::new_double_exp(0.5, 0.0, 0.3);
        assert!(dist2.is_err());
        
        // Invalid cutoff
        let dist3 = Distribution::new_double_exp(0.5, 2.0, 1.5);
        assert!(dist3.is_err());
        
        let dist4 = Distribution::new_double_exp(0.5, 2.0, -0.1);
        assert!(dist4.is_err());
    }

    #[test]
    fn test_distribution_sampling() {
        let mut rng = StdRng::seed_from_u64(42);
        
        let normal = Distribution::new_normal(0.0, 1.0).unwrap();
        let sample = normal.sample(&mut rng);
        assert!(sample.is_finite());
        
        let uniform = Distribution::new_uniform(-1.0, 1.0).unwrap();
        let sample = uniform.sample(&mut rng);
        assert!(sample >= -1.0 && sample <= 1.0);
        
        let exp = Distribution::new_exp(1.0).unwrap();
        let sample = exp.sample(&mut rng);
        assert!(sample >= 0.0);

        let poisson = Distribution::new_poisson(1.0).unwrap();
        let sample = poisson.sample(&mut rng);
        assert!(sample >= 0.0)
    }

    #[test]
    fn test_mutation_map_creation() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let test_dist = Distribution::new_double_exp(0.5, 2.0, 0.3).expect("Failed to create double exponential distribution for exon features");
        let test_seq = vec![1, 1, 4, 8, 2, 1, 2, 4];

        let map = MutationMap::new(1, 1, &test_seq,  &test_dist, &mut rng);
        assert_eq!(map.selection_dist_id, 1);
    }

    #[test]
    fn test_mutation_map_insert_and_get() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let test_dist = Distribution::new_double_exp(0.5, 2.0, 0.3).expect("Failed to create double exponential distribution for exon features");
        let test_seq = vec![1, 1, 4, 8, 2, 1, 2, 4];
        let mut map = MutationMap::new(0, 0, &test_seq, &test_dist, &mut rng);
        
        map.insert(1, 100, 0.5);
        let value = map.get(1, 100);
        assert_eq!(value, Some(&0.5));
        
        let missing = map.get(2, 100);
        assert_eq!(missing, None);
    }

    #[test]
    fn test_mutation_map_multiple_levels() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let test_dist = Distribution::new_double_exp(0.5, 2.0, 0.3).expect("Failed to create double exponential distribution for exon features");
        let test_seq = vec![1, 1, 4, 8, 2, 1, 2, 4];
        let mut map = MutationMap::new(0, 0, &test_seq, &test_dist, &mut rng);
        
        map.insert(1, 10, 0.1);
        map.insert(2, 10, 0.2);
        map.insert(4, 10, 0.3);
        map.insert(8, 10, 0.4);
        
        assert_eq!(map.get(1, 10), Some(&0.1));
        assert_eq!(map.get(4, 10), Some(&0.3));
        assert_eq!(map.get(2, 10), Some(&0.2));
        assert_eq!(map.get(8, 10), Some(&0.4));
    }

    #[test]
    fn test_n_allele_is_non_mutable_and_unmapped() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let selection_dist = Distribution::new_uniform(0.0, 1.0)
            .expect("failed to create selection distribution");
        let mu_dist = Distribution::new_poisson(10.0)
            .expect("failed to create mutation-rate distribution");
        let core_vec: Vec<Vec<u8>> = vec![vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4]];

        let mut seq = vec![16, 1, 16, 2, 4, 8, 16];
        let original_n_sites: Vec<u8> = seq.iter().copied().filter(|&x| x == 16).collect();

        let mut map = MutationMap::new(0, 0, &seq, &selection_dist, &mut rng);
        map.mutate(&core_vec, &mut seq, &selection_dist, &mu_dist);

        assert_eq!(seq[0], 16);
        assert_eq!(seq[2], 16);
        assert_eq!(seq[6], 16);
        assert!(map.get(16, 0).is_some());
        assert!(map.get(16, 2).is_some());
        assert!(map.get(16, 6).is_some());

        let post_n_sites: Vec<u8> = seq.iter().copied().filter(|&x| x == 16).collect();
        assert_eq!(original_n_sites, post_n_sites);
    }

    #[test]
    fn test_double_exponential_weight_update() {
        let mut double_exp = DoubleExponential::new(1.0, 2.0, 0.5, 0.0).unwrap();
        let uniform = Uniform::new(0.0, 1.0);
        let mut rng = StdRng::seed_from_u64(123);
        
        let old_weight = double_exp.weight;
        double_exp.weight(&uniform, &mut rng);
        let new_weight = double_exp.weight;
        
        // Weight should be in valid range
        assert!(new_weight >= 0.0 && new_weight <= 1.0);
        // With seeded RNG, weight should have changed (unless by coincidence)
        assert_ne!(old_weight, new_weight);
    }
}

