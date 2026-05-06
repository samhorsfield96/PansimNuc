// TODO each element mutation map needs to be shared globally, 
//at the moment it is only found in a specific individual and so if a mutation goes
// extinct, we lose tracking of that allele, effect should be maintained globally

use rand::Rng;
use rand::rngs::{StdRng, ThreadRng};
use rand::seq::IteratorRandom;
use rand_distr::{Distribution as RandDist, Exp, Normal, Uniform, Gamma};
use statrs::distribution::{Poisson, NegativeBinomial};
use std::collections::HashMap;
use std::fmt;
use std::os::unix::thread;
use crate::population::NucElement;

#[derive(Debug)]
pub enum DistributionError {
    InvalidNormalParameters,
    InvalidUniformParameters,
    InvalidExponentialParameters,
    InvalidDoubleExpParameters,
    InvalidPoissonParameters,
    InvalidGammaParameters,
    InvalidNegativeBinomialParameters,
}

impl fmt::Display for DistributionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            DistributionError::InvalidNormalParameters => {
                write!(f, "Invalid Normal distribution: std_dev must be positive")
            }
            DistributionError::InvalidUniformParameters => {
                write!(
                    f,
                    "Invalid Uniform distribution: low must be less than high"
                )
            }
            DistributionError::InvalidExponentialParameters => {
                write!(
                    f,
                    "Invalid Exponential distribution: lambda must be positive"
                )
            }
            DistributionError::InvalidDoubleExpParameters => {
                write!(
                    f,
                    "Invalid DoubleExp distribution: lambdas must be positive and cutoff must be between 0 and 1"
                )
            }
            DistributionError::InvalidPoissonParameters => {
                write!(f, "Invalid Poisson distribution: lambda must be positive")
            }
            DistributionError::InvalidNegativeBinomialParameters => {
                write!(f, "Invalid Negative Binomial distribution: parameters must be positive")
            }
            DistributionError::InvalidGammaParameters => {
                write!(f, "Invalid Gamma distribution: shape and scale must be positive")
            }
        }
    }
}

impl std::error::Error for DistributionError {}

#[derive(Debug)]
pub enum DistributionConfigError {
    MissingDistributionType {
        section: String,
        distribution_key: String,
        legacy_key: String,
    },
    MissingParameter {
        section: String,
        distribution: String,
        parameter: String,
    },
    InvalidParameterValue {
        section: String,
        key: String,
        value: String,
    },
    UnsupportedDistribution {
        section: String,
        distribution: String,
    },
    InvalidDistributionParameters {
        section: String,
        distribution: String,
        source: DistributionError,
    },
}

impl fmt::Display for DistributionConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DistributionConfigError::MissingDistributionType {
                section,
                distribution_key,
                legacy_key,
            } => write!(
                f,
                "Missing selection distribution for section '{}'. Set '{}' or use legacy key '{}'.",
                section, distribution_key, legacy_key
            ),
            DistributionConfigError::MissingParameter {
                section,
                distribution,
                parameter,
            } => write!(
                f,
                "Missing required parameter '{}' for selection distribution '{}' in section '{}'.",
                parameter, distribution, section
            ),
            DistributionConfigError::InvalidParameterValue {
                section,
                key,
                value,
            } => write!(
                f,
                "Invalid value '{}' for config key '{}' in section '{}'. Expected a float.",
                value, key, section
            ),
            DistributionConfigError::UnsupportedDistribution {
                section,
                distribution,
            } => write!(
                f,
                "Unsupported selection distribution '{}' in section '{}'. Supported values: normal, uniform, exp, double_exp, poisson, gamma.",
                distribution, section
            ),
            DistributionConfigError::InvalidDistributionParameters {
                section,
                distribution,
                source,
            } => write!(
                f,
                "Invalid parameters for selection distribution '{}' in section '{}': {}",
                distribution, section, source
            ),
        }
    }
}

impl std::error::Error for DistributionConfigError {}

#[derive(Debug, Clone)]
pub struct DoubleExponential {
    exp1: Exp<f64>,
    exp2: Exp<f64>,
    cutoff: f64,
    weight_rng: Uniform<f64>,
}

impl DoubleExponential {
    pub fn new(
        lambda1: f64,
        lambda2: f64,
        cutoff: f64,
    ) -> Result<Self, DistributionError> {
        if lambda1 <= 0.0 || lambda2 <= 0.0 || cutoff < 0.0 || cutoff > 1.0 {
            return Err(DistributionError::InvalidDoubleExpParameters);
        }

        let weight_rng = Uniform::new(0.0, 1.0);

        let exp1 = Exp::new(lambda1).map_err(|_| DistributionError::InvalidDoubleExpParameters)?;
        let exp2 = Exp::new(lambda2).map_err(|_| DistributionError::InvalidDoubleExpParameters)?;

        Ok(Self {
            exp1,
            exp2,
            cutoff,
            weight_rng,
        })
    }

    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        let weight = self.weight_rng.sample(rng);

        if weight <= self.cutoff {
            // Higher selected gene
            return self.exp2.sample(rng);
        } else {
            // lower selected gene
            // keep sampling until <= 1, mimicking selection coefficient sampling
            let mut selection_coefficient = 2.0;

            while selection_coefficient > 1.0 {
                selection_coefficient = self.exp1.sample(rng);
            }
            selection_coefficient *= -1.0; // make negative for deleterious mutations
            return selection_coefficient;
        }
    }
}

#[derive(Clone)]
pub enum Distribution {
    Normal(Normal<f64>),
    Uniform(Uniform<f64>),
    Exp(Exp<f64>),
    DoubleExp(DoubleExponential),
    Poisson(Poisson),
    Gamma(Gamma<f64>),
    NegativeBinomial(NegativeBinomial),
}

impl Distribution {
    fn parse_required_f64(
        configuration: &HashMap<String, String>,
        section: &str,
        distribution: &str,
        key: &str,
    ) -> Result<f64, DistributionConfigError> {
        let value =
            configuration
                .get(key)
                .ok_or_else(|| DistributionConfigError::MissingParameter {
                    section: section.to_string(),
                    distribution: distribution.to_string(),
                    parameter: key.rsplit('.').next().unwrap_or(key).to_string(),
                })?;

        value
            .parse::<f64>()
            .map_err(|_| DistributionConfigError::InvalidParameterValue {
                section: section.to_string(),
                key: key.to_string(),
                value: value.clone(),
            })
    }

    pub fn from_selection_config(
        configuration: &HashMap<String, String>,
        section: &str,
    ) -> Result<Self, DistributionConfigError> {
        let distribution_key = format!("{}.selection_distribution", section);
        let legacy_key = format!("{}.selection_coefficient", section);

        let raw_distribution = if let Some(distribution) = configuration.get(&distribution_key) {
            distribution.clone()
        } else if configuration.contains_key(&legacy_key) {
            "exp".to_string()
        } else {
            return Err(DistributionConfigError::MissingDistributionType {
                section: section.to_string(),
                distribution_key,
                legacy_key,
            });
        };

        let normalized_distribution = raw_distribution
            .trim()
            .to_ascii_lowercase()
            .replace('-', "_");

        match normalized_distribution.as_str() {
            "normal" => {
                let mean = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_mean", section),
                )?;
                let std_dev = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_std_dev", section),
                )?;

                Self::new_normal(mean, std_dev).map_err(|source| {
                    DistributionConfigError::InvalidDistributionParameters {
                        section: section.to_string(),
                        distribution: raw_distribution,
                        source,
                    }
                })
            }
            "uniform" => {
                let low = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_low", section),
                )?;
                let high = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_high", section),
                )?;

                Self::new_uniform(low, high).map_err(|source| {
                    DistributionConfigError::InvalidDistributionParameters {
                        section: section.to_string(),
                        distribution: raw_distribution,
                        source,
                    }
                })
            }
            "exp" | "exponential" => {
                let lambda_key = if configuration.contains_key(&distribution_key) {
                    format!("{}.selection_lambda", section)
                } else {
                    legacy_key.clone()
                };
                let lambda = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &lambda_key,
                )?;

                Self::new_exp(lambda).map_err(|source| {
                    DistributionConfigError::InvalidDistributionParameters {
                        section: section.to_string(),
                        distribution: raw_distribution,
                        source,
                    }
                })
            }
            "double_exp" | "doubleexponential" | "double_exponential" => {
                let lambda1 = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_lambda1", section),
                )?;
                let lambda2 = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_lambda2", section),
                )?;
                let cutoff = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_cutoff", section),
                )?;

                Self::new_double_exp(lambda1, lambda2, cutoff).map_err(|source| {
                    DistributionConfigError::InvalidDistributionParameters {
                        section: section.to_string(),
                        distribution: raw_distribution,
                        source,
                    }
                })
            }
            "poisson" => {
                let lambda = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_lambda", section),
                )?;

                Self::new_poisson(lambda).map_err(|source| {
                    DistributionConfigError::InvalidDistributionParameters {
                        section: section.to_string(),
                        distribution: raw_distribution,
                        source,
                    }
                })
            }
            "negbinom" | "negativebinomial" | "negative_binomial" => {
                let r = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_r", section),
                )?;
                let p = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_p", section),
                )?;

                Self::new_negative_binomial(r, p).map_err(|source| {
                    DistributionConfigError::InvalidDistributionParameters {
                        section: section.to_string(),
                        distribution: raw_distribution,
                        source,
                    }
                })
            }
            "gamma" => {
                let shape = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_shape", section),
                )?;
                let scale = Self::parse_required_f64(
                    configuration,
                    section,
                    &raw_distribution,
                    &format!("{}.selection_scale", section),
                )?;

                Self::new_gamma(shape, scale).map_err(|source| {
                    DistributionConfigError::InvalidDistributionParameters {
                        section: section.to_string(),
                        distribution: raw_distribution,
                        source,
                    }
                })
            }
            _ => Err(DistributionConfigError::UnsupportedDistribution {
                section: section.to_string(),
                distribution: raw_distribution,
            }),
        }
    }

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

    pub fn new_double_exp(
        lambda1: f64,
        lambda2: f64,
        cutoff: f64,
    ) -> Result<Self, DistributionError> {
        DoubleExponential::new(lambda1, lambda2, cutoff).map(Distribution::DoubleExp)
    }

    pub fn new_poisson(lambda: f64) -> Result<Self, DistributionError> {
        Poisson::new(lambda)
            .map(Distribution::Poisson)
            .map_err(|_| DistributionError::InvalidPoissonParameters)
    }

    pub fn new_negative_binomial(r: f64, p: f64) -> Result<Self, DistributionError> {
        if r <= 0.0 || p <= 0.0 || p >= 1.0 {
            return Err(DistributionError::InvalidNegativeBinomialParameters);
        }
        NegativeBinomial::new(r, p)
            .map(Distribution::NegativeBinomial)
            .map_err(|_| DistributionError::InvalidNegativeBinomialParameters)
    }

    pub fn new_gamma(shape: f64, scale: f64) -> Result<Self, DistributionError> {
        Gamma::new(shape, scale)
            .map(Distribution::Gamma)
            .map_err(|_| DistributionError::InvalidGammaParameters)
    }

    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        match self {
            Distribution::Normal(d) => d.sample(rng),
            Distribution::Uniform(d) => d.sample(rng),
            Distribution::Exp(d) => d.sample(rng),
            Distribution::DoubleExp(d) => d.sample(rng),
            Distribution::Poisson(d) => d.sample(rng),
            Distribution::Gamma(d) => d.sample(rng),
            Distribution::NegativeBinomial(d) => d.sample(rng) as f64,
        }
    }
}

#[derive(Clone)]
pub struct MutationMap {
    pub selection_dist_id: usize,
    pub mu_dist_id: usize,
    data: [Vec<Option<f64>>; 5],
}

impl MutationMap {
    fn allele_to_index(level: u8) -> Option<usize> {
        match level {
            1 => Some(0),
            2 => Some(1),
            4 => Some(2),
            8 => Some(3),
            16 => Some(4),
            _ => panic!("Allele code must be one-hot (1, 2, 4, 8, 16); got {}", level),
        }
    }

    pub fn new(
        selection_dist_id: usize,
        mu_dist_id: usize,
        seq: &Vec<u8>,
        selection_dist: &Distribution,
        rng: &mut StdRng,
    ) -> Self {
        let mut data = std::array::from_fn(|_| vec![None; seq.len()]);

        for (site, allele) in seq.iter().enumerate() {
            let allele_index = Self::allele_to_index(*allele)
                .expect("Allele code conversion failed while building mutation map");
            data[allele_index][site] = Some(selection_dist.sample(rng));
        }

        Self {
            selection_dist_id,
            mu_dist_id,
            data,
        }
    }

    fn update_data (&mut self, site: usize, is_insertion: bool, sequence_length: usize) {
        // if is_insertion, duplicate the current site entry and shift all later sites up by one
        if is_insertion {
            for allele_map in self.data.iter_mut() {
                let existing_value = allele_map.get(site).copied().flatten();
                allele_map.insert(site, existing_value);
            }
        } else {
            // deletion: remove this site and shift later sites down by one
            for allele_map in self.data.iter_mut() {
                if site < sequence_length {
                    allele_map.remove(site);
                }
            }
        }
    }

    fn insert(&mut self, level: u8, key: usize, value: f64) {
        let allele_index = Self::allele_to_index(level)
            .expect("Cannot insert selection coefficient for invalid allele code");
        if self.data[allele_index].len() <= key {
            self.data[allele_index].resize(key + 1, None);
        }
        self.data[allele_index][key] = Some(value);
    }

    pub fn get(&self, level: u8, key: usize) -> Option<&f64> {
        let allele_index = Self::allele_to_index(level)
            .expect("Cannot lookup selection coefficient for invalid allele code");
        self.data[allele_index].get(key).and_then(Option::as_ref)
    }

    #[cfg(test)]
    pub(crate) fn set_for_test(&mut self, level: u8, key: usize, value: f64) {
        self.insert(level, key, value);
    }

    fn mutate_snps(
        &mut self,
        core_vec: &Vec<Vec<u8>>,
        seq: &mut Vec<u8>,
        selection_dist: &Distribution,
        n_sites: usize,
        thread_rng: &mut ThreadRng
    ) -> usize {
        let seq_len = seq.len();
        if seq_len == 0 {
            return 0;
        }

        let n_draws = n_sites.min(seq_len);
        let sampled_sites = (0..seq_len).choose_multiple(thread_rng, n_draws);

        // iterate for number of mutations required to reach mutation rate
        for mutant_site in sampled_sites {
            // sample new site to mutate
            let value = seq[mutant_site];

            if value == 16 {
                continue;
            }

            let allele_index =
                Self::allele_to_index(value).expect("Sampled non-mutable N allele for mutation");
            let values = &core_vec[allele_index];

            // sample new allele
            let new_allele = values[thread_rng.gen_range(0..values.len())];

            // generate new selection coefficient for this mutation if necessary, otherwise retrieve existing one
            if let Some(_) = self.get(new_allele, mutant_site) {
                // value exists
                continue;
            } else {
                let selection_coefficient = selection_dist.sample(thread_rng);
                self.insert(new_allele, mutant_site, selection_coefficient);
            }

            // set value in place
            seq[mutant_site] = new_allele;
        }
        n_draws
    }

    fn mutate_indels (
        &mut self,
        core_vec: &Vec<Vec<u8>>,
        seq: &mut Vec<u8>,
        selection_dist: &Distribution,
        n_sites: usize,
        thread_rng: &mut ThreadRng
    ) -> usize {
        // iterate for number of mutations required to reach mutation rate
        // needs to be dynamic to allow for changes in sequence length from indels
        for _ in 0..n_sites {
            // sample new site to mutate
            let seq_len = seq.len();

            let mut mutant_site = 0;
            let mut value = 1;
            
            // protect against empty sequence edge case where no deletions can occur
            if seq_len != 0 {
                mutant_site = thread_rng.gen_range(0..seq_len);
                value = seq[mutant_site];
            }

            // N is stored in the map but should not be mutated; skip
            if value == 16 {
                continue;
            }

            // determine whether will be insertion or deletion, set to insertion if empty already
            let mut is_insertion = true;
            if seq_len != 0 {
                is_insertion = thread_rng.gen_bool(0.5);
            }

            // update data for selection coefficient
            self.update_data(mutant_site, is_insertion, seq_len);

            if is_insertion {
                // insertion: sample new allele to insert
                let new_allele = core_vec[4][thread_rng.gen_range(0..core_vec[4].len())];
                seq.insert(mutant_site, new_allele);

                // add selection coefficient
                let selection_coefficient = selection_dist.sample(thread_rng);
                self.insert(new_allele, mutant_site, selection_coefficient);
            } else {
                // deletion: remove allele at mutant_site
                seq.remove(mutant_site);
            }
        }
        n_sites
    }

    pub fn mutate(
        &mut self,
        core_vec: &Vec<Vec<u8>>,
        seq: &mut Vec<u8>,
        original_length: usize,
        frameshift: &mut bool,
        selection_dist: &Distribution,
        n_snps: usize,
        n_indels: usize,
        thread_rng: &mut ThreadRng
    ) -> (usize, usize) {
        // mutate SNPs
        let n_snps = self.mutate_snps(core_vec, seq, selection_dist, n_snps, thread_rng);
        let n_indels = self.mutate_indels(core_vec, seq, selection_dist, n_indels, thread_rng);

        // update frameshift status
        if n_indels > 0 {
            let max_length = usize::max(original_length, seq.len());
            let min_length = usize::min(original_length, seq.len());
            if (max_length - min_length) % 3 != 0 {
                *frameshift = true;
            } else {
                *frameshift = false;
            }
        }
        

        (n_snps, n_indels)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use std::collections::HashMap;

    fn selection_config(section: &str, entries: &[(&str, &str)]) -> HashMap<String, String> {
        let mut config = HashMap::new();
        for (key, value) in entries {
            config.insert(format!("{}.{}", section, key), (*value).to_string());
        }
        config
    }

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
    fn test_negative_binomial_creation() {
        let dist = Distribution::new_negative_binomial(5.0, 0.5);
        assert!(dist.is_ok());
    }

    #[test]
    fn test_negative_binomial_invalid_params() {
        let dist = Distribution::new_negative_binomial(-5.0, 0.5);
        assert!(dist.is_err());

        let dist2 = Distribution::new_negative_binomial(5.0, -0.5);
        assert!(dist2.is_err());
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
    fn test_double_exp_distribution_generates_negative_weights() {
        let mut rng = StdRng::seed_from_u64(42);
        let dist = Distribution::new_double_exp(0.5, 0.5, 0.0).unwrap();

        // Generate a large number of samples to ensure we get some from the exp1 distribution
        let samples: Vec<f64> = (0..20).map(|_| dist.sample(&mut rng)).collect();

        println!("Sampled values: {:?}", &samples[0..20]); // Print first 20 samples for inspection

        // Check that we have some negative values, which would indicate sampling from exp1
        assert!(samples.iter().all(|&x| x < 0.0));
    }

    #[test]
    fn test_gamma_distribution_creation() {
        let dist = Distribution::new_gamma(2.0, 3.0);
        assert!(dist.is_ok());
    }

    #[test]
    fn test_gamma_distribution_invalid_params() {
        let dist = Distribution::new_gamma(-2.0, 3.0);
        assert!(dist.is_err());

        let dist2 = Distribution::new_gamma(2.0, -3.0);
        assert!(dist2.is_err());
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
        assert!(sample >= 0.0);

        let gamma = Distribution::new_gamma(2.0, 3.0).unwrap();
        let sample = gamma.sample(&mut rng);
        assert!(sample >= 0.0);

        let neg_binomial = Distribution::new_negative_binomial(5.0, 0.5).unwrap();
        let sample = neg_binomial.sample(&mut rng);
        assert!(sample >= 0.0);
    }

    #[test]
    fn test_selection_distribution_from_config_normal() {
        let config = selection_config(
            "exons",
            &[
                ("selection_distribution", "normal"),
                ("selection_mean", "0.1"),
                ("selection_std_dev", "0.2"),
            ],
        );

        let result = Distribution::from_selection_config(&config, "exons");
        assert!(matches!(result, Ok(Distribution::Normal(_))));
    }

    #[test]
    fn test_selection_distribution_from_config_uniform() {
        let config = selection_config(
            "introns",
            &[
                ("selection_distribution", "uniform"),
                ("selection_low", "0.0"),
                ("selection_high", "1.0"),
            ],
        );

        let result = Distribution::from_selection_config(&config, "introns");
        assert!(matches!(result, Ok(Distribution::Uniform(_))));
    }

    #[test]
    fn test_selection_distribution_from_config_exp_legacy_key() {
        let config = selection_config("intergenic", &[("selection_coefficient", "0.02")]);

        let result = Distribution::from_selection_config(&config, "intergenic");
        assert!(matches!(result, Ok(Distribution::Exp(_))));
    }

    #[test]
    fn test_selection_distribution_from_config_double_exp() {
        let config = selection_config(
            "TE-CUT",
            &[
                ("selection_distribution", "double_exp"),
                ("selection_lambda1", "0.5"),
                ("selection_lambda2", "2.0"),
                ("selection_cutoff", "0.3"),
            ],
        );

        let result = Distribution::from_selection_config(&config, "TE-CUT");
        assert!(matches!(result, Ok(Distribution::DoubleExp(_))));
    }

    #[test]
    fn test_selection_distribution_from_config_poisson() {
        let config = selection_config(
            "TE-COPY",
            &[
                ("selection_distribution", "poisson"),
                ("selection_lambda", "1.0"),
            ],
        );

        let result = Distribution::from_selection_config(&config, "TE-COPY");
        assert!(matches!(result, Ok(Distribution::Poisson(_))));
    }

    #[test]
    fn test_selection_distribution_from_config_negative_binomial() {
        let config = selection_config(
            "TE-COPY",
            &[
                ("selection_distribution", "negative_binomial"),
                ("selection_r", "5.0"),
                ("selection_p", "0.5"),
            ],
        );

        let result = Distribution::from_selection_config(&config, "TE-COPY");
        assert!(matches!(result, Ok(Distribution::NegativeBinomial(_))));
    }

    #[test]
    fn test_selection_distribution_from_config_missing_parameter() {
        let config = selection_config(
            "exons",
            &[
                ("selection_distribution", "normal"),
                ("selection_mean", "0.1"),
            ],
        );

        let result = Distribution::from_selection_config(&config, "exons");
        assert!(matches!(
            result,
            Err(DistributionConfigError::MissingParameter { .. })
        ));
    }

    #[test]
    fn test_mutation_map_creation() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let test_dist = Distribution::new_double_exp(0.5, 2.0, 0.3)
            .expect("Failed to create double exponential distribution for exon features");
        let test_seq = vec![1, 1, 4, 8, 2, 1, 2, 4];

        let map = MutationMap::new(1, 1, &test_seq, &test_dist, &mut rng);
        assert_eq!(map.selection_dist_id, 1);
    }

    #[test]
    fn test_mutation_map_insert_and_get() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let test_dist = Distribution::new_double_exp(0.5, 2.0, 0.3)
            .expect("Failed to create double exponential distribution for exon features");
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
        let test_dist = Distribution::new_double_exp(0.5, 2.0, 0.3)
            .expect("Failed to create double exponential distribution for exon features");
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
        let mut thread_rng = rand::thread_rng();
        let selection_dist =
            Distribution::new_uniform(0.0, 1.0).expect("failed to create selection distribution");
        let mu_dist =
            Distribution::new_poisson(10.0).expect("failed to create mutation-rate distribution");
        let indel_dist =
            Distribution::new_poisson(1e-12).expect("failed to create mutation-rate distribution");
        let core_vec: Vec<Vec<u8>> =
            vec![vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4]];

        let mut seq = vec![16, 1, 16, 2, 4, 8, 16];
        let original_length = seq.len();
        let original_n_sites: Vec<u8> = seq.iter().copied().filter(|&x| x == 16).collect();

        let n_snps = mu_dist.sample(&mut rng) as usize;
        let n_indels = indel_dist.sample(&mut rng) as usize;

        let mut map = MutationMap::new(0, 0, &seq, &selection_dist, &mut rng);
        map.mutate(&core_vec, &mut seq, original_length, &mut false, &selection_dist, n_snps, n_indels, &mut thread_rng);

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
    // test that map is updated many times, and that selection cofficients generated remain same after many mutations
    fn test_mutation_map_consistency_after_mutations() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let mut thread_rng = rand::thread_rng();
        let selection_dist =
            Distribution::new_uniform(1e-10, 1e10).expect("failed to create selection distribution");
        let mu_dist =
            Distribution::new_poisson(1e3).expect("failed to create mutation-rate distribution");
        let indel_dist =
            Distribution::new_poisson(1e-30).expect("failed to create mutation-rate distribution");
        let core_vec: Vec<Vec<u8>> =
            vec![vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4]];

        let mut seq = vec![1, 1, 4, 8, 2, 1, 2, 4];
        let original_length = seq.len();

        let mut map = MutationMap::new(0, 0, &seq, &selection_dist, &mut rng);
        let original_map_state = map.data.clone();

        // mutate many times with very low mutation rates to ensure map is updated but sequence does not change
        for _ in 0..100 {
            let n_snps = mu_dist.sample(&mut rng) as usize;
            let n_indels = indel_dist.sample(&mut rng) as usize;
            map.mutate(&core_vec, &mut seq, original_length, &mut false, &selection_dist, n_snps, n_indels, &mut thread_rng);
            assert_eq!(seq.len(), original_length);
            assert_eq!(map.data.len(), original_map_state.len());
            for (allele_map, original_allele_map) in map.data.iter().zip(original_map_state.iter()) {
                // for each pre-existing entry, check that the same key still has the same value
                for (key, value) in original_allele_map.iter().enumerate() {
                    if let Some(value) = value {
                        assert_eq!(allele_map.get(key), Some(&Some(*value)));
                    }
                }
            }
        }
    }

    #[test]
    fn test_indels_change_sequence_length() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let mut thread_rng = rand::thread_rng();
        let dist = Distribution::new_uniform(0.0, 1.0).unwrap();
        // Large sequence so many indels fire; zero SNP rate so only indels mutate
        let seq: Vec<u8> = vec![1u8; 100];
        let mut map = MutationMap::new(0, 0, &seq, &dist, &mut rng);
        let mut seq_mut = seq.clone();

        let mu_dist = Distribution::new_poisson(1e-12).unwrap();
        // Force many insertions by biasing gen_bool via a deterministic seed that
        // reliably produces insertions; use a very high rate to guarantee length change
        let indel_dist = Distribution::new_poisson(10.0).unwrap();
        let core_vec: Vec<Vec<u8>> = vec![
            vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4], vec![1, 2, 4, 8, 16],
        ];

        let n_snps = mu_dist.sample(&mut rng) as usize;
        let n_indels = indel_dist.sample(&mut rng) as usize;
        map.mutate(&core_vec, &mut seq_mut, seq.len(), &mut false, &dist, n_snps, n_indels, &mut thread_rng);

        // With 1000 indels on a seed that produces ~50% insertions, final length must differ
        assert_ne!(seq_mut.len(), seq.len());
    }

    #[test]
    fn test_insertion_shifts_selection_coefficients_up() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let dist = Distribution::new_uniform(0.0, 1.0).unwrap();
        // All-A sequence so data[0] has a dense entry at every site
        let seq = vec![1u8, 1, 1, 1];
        let mut map = MutationMap::new(0, 0, &seq, &dist, &mut rng);
        map.set_for_test(1, 0, 0.10);
        map.set_for_test(1, 1, 0.20);
        map.set_for_test(1, 2, 0.30);
        map.set_for_test(1, 3, 0.40);

        map.set_for_test(2, 0, 0.10);
        map.set_for_test(2, 1, 0.20);
        map.set_for_test(2, 2, 0.30);
        map.set_for_test(2, 3, 0.40);

        map.set_for_test(4, 0, 0.10);
        map.set_for_test(4, 1, 0.20);
        map.set_for_test(4, 2, 0.30);
        map.set_for_test(4, 3, 0.40);

        map.set_for_test(8, 0, 0.10);
        map.set_for_test(8, 1, 0.20);
        map.set_for_test(8, 2, 0.30);
        map.set_for_test(8, 3, 0.40);

        // Inserting at site 1 should shift coefficients at sites 1+ up by one
        map.update_data(1, true, 4);

        // check insertion correct
        assert_ne!(map.get(1, 1), Some(&0.10)); // site 1 changed
        assert_ne!(map.get(2, 1), Some(&0.10)); // site 1 changed
        assert_ne!(map.get(4, 1), Some(&0.10)); // site 1 changed
        assert_ne!(map.get(8, 1), Some(&0.10)); // site 1 changed

        assert_ne!(map.get(1, 1), None); // site 1 changed
        assert_ne!(map.get(2, 1), None); // site 1 changed
        assert_ne!(map.get(4, 1), None); // site 1 changed
        assert_ne!(map.get(8, 1), None); // site 1 changed


        assert_eq!(map.get(1, 0), Some(&0.10)); // site 0 unchanged
        assert_eq!(map.get(1, 2), Some(&0.20)); // site 1 shifted to 2
        assert_eq!(map.get(1, 3), Some(&0.30)); // site 2 shifted to 3
        assert_eq!(map.get(1, 4), Some(&0.40)); // site 3 shifted to 4
        assert_eq!(map.get(2, 0), Some(&0.10)); // site 0 unchanged
        assert_eq!(map.get(2, 2), Some(&0.20)); // site 1 shifted to 2
        assert_eq!(map.get(2, 3), Some(&0.30)); // site 2 shifted to 3
        assert_eq!(map.get(2, 4), Some(&0.40)); // site 3 shifted to 4
        assert_eq!(map.get(4, 0), Some(&0.10)); // site 0 unchanged
        assert_eq!(map.get(4, 2), Some(&0.20)); // site 1 shifted to 2
        assert_eq!(map.get(4, 3), Some(&0.30)); // site 2 shifted to 3
        assert_eq!(map.get(4, 4), Some(&0.40)); // site 3 shifted to 4
        assert_eq!(map.get(8, 0), Some(&0.10)); // site 0 unchanged
        assert_eq!(map.get(8, 2), Some(&0.20)); // site 1 shifted to 2
        assert_eq!(map.get(8, 3), Some(&0.30)); // site 2 shifted to 3
        assert_eq!(map.get(8, 4), Some(&0.40)); // site 3 shifted to 4

    }

    #[test]
    fn test_deletion_removes_selection_coefficient_at_deleted_site() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let dist = Distribution::new_uniform(0.0, 1.0).unwrap();
        // All-A sequence so data[0] has a dense entry at every site
        let seq = vec![1u8, 1, 1, 1];
        let mut map = MutationMap::new(0, 0, &seq, &dist, &mut rng);
        map.set_for_test(1, 0, 0.10);
        map.set_for_test(1, 1, 0.20);
        map.set_for_test(1, 2, 0.30);
        map.set_for_test(1, 3, 0.40);

        map.set_for_test(2, 0, 0.10);
        map.set_for_test(2, 1, 0.20);
        map.set_for_test(2, 2, 0.30);
        map.set_for_test(2, 3, 0.40);

        map.set_for_test(4, 0, 0.10);
        map.set_for_test(4, 1, 0.20);
        map.set_for_test(4, 2, 0.30);
        map.set_for_test(4, 3, 0.40);

        map.set_for_test(8, 0, 0.10);
        map.set_for_test(8, 1, 0.20);
        map.set_for_test(8, 2, 0.30);
        map.set_for_test(8, 3, 0.40);

        // Deleting site 2 should remove its coefficient from the map
        map.update_data(2, false, 4);

        assert_eq!(map.get(1, 0), Some(&0.10)); // site 0 unchanged
        assert_eq!(map.get(1, 1), Some(&0.20)); // site 1 unchanged
        assert_eq!(map.get(1, 2), Some(&0.40)); // site 3 coefficient shifted down to site 2
        assert_eq!(map.get(1, 3), None); // site 3 deleted, should be removed from map

        assert_eq!(map.get(2, 0), Some(&0.10)); // site 0 unchanged
        assert_eq!(map.get(2, 1), Some(&0.20)); // site 1 unchanged
        assert_eq!(map.get(2, 2), Some(&0.40)); // site 3 coefficient shifted down to site 2
        assert_eq!(map.get(2, 3), None); // site 3 deleted, should be removed from map

        assert_eq!(map.get(4, 0), Some(&0.10)); // site 0 unchanged
        assert_eq!(map.get(4, 1), Some(&0.20)); // site 1 unchanged
        assert_eq!(map.get(4, 2), Some(&0.40)); // site 3 coefficient shifted down to site 2
        assert_eq!(map.get(4, 3), None); // site 3 deleted, should be removed from map

        assert_eq!(map.get(8, 0), Some(&0.10)); // site 0 unchanged
        assert_eq!(map.get(8, 1), Some(&0.20)); // site 1 unchanged
        assert_eq!(map.get(8, 2), Some(&0.40)); // site 3 coefficient shifted down to site 2
        assert_eq!(map.get(8, 3), None); // site 3 deleted, should be removed from map
    }

    #[test]
    fn test_frameshift_flag_matches_length_change_mod3() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let mut thread_rng = rand::thread_rng();
        let dist = Distribution::new_uniform(0.0, 1.0).unwrap();
        let seq: Vec<u8> = vec![1u8; 100];
        let mut map = MutationMap::new(0, 0, &seq, &dist, &mut rng);
        let mut seq_mut = seq.clone();
        let original_length = seq.len();

        let mu_dist = Distribution::new_poisson(1e-12).unwrap();
        let indel_dist = Distribution::new_poisson(10.0).unwrap();
        let core_vec: Vec<Vec<u8>> = vec![
            vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4], vec![1, 2, 4, 8, 16],
        ];

        let mut frameshift = true;
        let n_snps = mu_dist.sample(&mut thread_rng) as usize;
        let n_indels = indel_dist.sample(&mut thread_rng) as usize;
        map.mutate(&core_vec, &mut seq_mut, original_length, &mut frameshift, &dist, n_snps, n_indels, &mut thread_rng);

        let expected = seq_mut.len().abs_diff(original_length) % 3 != 0;
        println!("Original length: {}, New length: {}, Frameshift: {}, Expected frameshift: {}", original_length, seq_mut.len(), frameshift, expected);
        assert_eq!(frameshift, expected);
    }

    #[test]
    fn test_frameshift_not_updated_when_no_indels() {
        let mut rng: StdRng = StdRng::seed_from_u64(42);
        let mut thread_rng = rand::thread_rng();
        let dist = Distribution::new_uniform(0.0, 1.0).unwrap();
        let seq: Vec<u8> = vec![1u8; 12];
        let mut map = MutationMap::new(0, 0, &seq, &dist, &mut rng);
        let mut seq_mut = seq.clone();

        let mu_dist = Distribution::new_poisson(1e-12).unwrap();
        let indel_dist = Distribution::new_poisson(1e-12).unwrap();
        let core_vec: Vec<Vec<u8>> = vec![
            vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4], vec![1, 2, 4, 8, 16],
        ];

        // Start with frameshift already set; expect it to remain unchanged when no indels fire
        let mut frameshift = true;
        let n_snps = mu_dist.sample(&mut thread_rng) as usize;
        let n_indels = indel_dist.sample(&mut thread_rng) as usize;
        map.mutate(&core_vec, &mut seq_mut, seq.len(), &mut frameshift, &dist, n_snps, n_indels, &mut thread_rng);

        assert!(frameshift);
    }
}
