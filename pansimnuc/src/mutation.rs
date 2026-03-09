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
        // TODO generate uniform distribution before this that is passed between samples, avoiding regenerating each time
        //let uniform = Uniform::new(0.0, 1.0).map_err(|_| DistributionError::InvalidDoubleExpParameters)?;
        
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

    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        match self {
            Distribution::Normal(d) => d.sample(rng),
            Distribution::Uniform(d) => d.sample(rng),
            Distribution::Exp(d) => d.sample(rng),
            Distribution::DoubleExp(d) => d.sample(rng),
        }
    }
}

pub struct MutationMap {
    pub distribution_id: usize,
    data: [FxHashMap<usize, f64>; 4],
}

impl MutationMap {
    pub fn new(distribution_id: usize) -> Self {
        Self {distribution_id, data: std::array::from_fn(|_| FxHashMap::default()) }
    }

    fn insert(&mut self, level: u8, key: usize, value: f64) {
        self.data[level as usize].insert(key, value);
    }

    fn get(&self, level: u8, key: usize) -> Option<&f64> {
        self.data[level as usize].get(&key)
    }    

    fn mutate (& mut self, poisson: &mut Poisson, core_vec: &Vec<Vec<u8>>, seq: &mut Vec<u8>, distribution: &Distribution) {
        // thread-specific random number generator
        let mut thread_rng = rand::thread_rng();

        let seq_len = seq.len();

        // sample from Poisson distribution for number of sites to mutate in this isolate
        let n_sites = poisson.sample(&mut thread_rng) as usize;
        let sampled_sites = (0..seq_len).choose_multiple(&mut thread_rng, n_sites);

        // iterate for number of mutations required to reach mutation rate
        for mutant_site in sampled_sites {
            // sample new site to mutate
            let value = seq[mutant_site];
            
            let values = &core_vec[value as usize];

            // sample new allele
            let new_allele = values.iter().choose_multiple(&mut thread_rng, 1)[0];
            let mut selection_coefficient: f64 = 0.0;

            // generate new selection coefficient for this mutation if necessary, otherwise retrieve existing one
            if let Some(coeff) = self.get(*new_allele, mutant_site) {
                // value exists
                selection_coefficient = *coeff;
            } else {
                selection_coefficient = distribution.sample(&mut thread_rng);
            }

            // set value in place
            seq[mutant_site] = *new_allele;
            self.insert(*new_allele, mutant_site, selection_coefficient);

        }
    }
}

