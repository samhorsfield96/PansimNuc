use rand::Rng;
use rand::rngs::StdRng;
use rand::seq::IteratorRandom;
use rand_distr::{Distribution as RandDist, Exp, Normal, Uniform, Gamma};
use rustc_hash::FxHashMap;
use statrs::distribution::Poisson;
use std::collections::HashMap;
use std::fmt;

#[derive(Debug)]
pub enum DistributionError {
    InvalidNormalParameters,
    InvalidUniformParameters,
    InvalidExponentialParameters,
    InvalidDoubleExpParameters,
    InvalidPoissonParameters,
    InvalidGammaParameters,
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
        let mut selection_coefficient = 0.0;

        if weight <= self.cutoff {
            // Higher selected gene
            selection_coefficient = self.exp2.sample(rng);
        } else {
            // lower selected gene
            // keep sampling until below 1, mimicking selection coefficient sampling
            while selection_coefficient >= 1.0 {
                selection_coefficient = self.exp1.sample(rng);
            }
        }

        selection_coefficient
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

    pub fn new(
        selection_dist_id: usize,
        mu_dist_id: usize,
        seq: &Vec<u8>,
        selection_dist: &Distribution,
        rng: &mut StdRng,
    ) -> Self {
        let mut data = std::array::from_fn(|_| FxHashMap::default());

        for (site, allele) in seq.iter().enumerate() {
            let allele_index = Self::allele_to_index(*allele)
                .expect("Allele code conversion failed while building mutation map");
            data[allele_index].insert(site, selection_dist.sample(rng));
        }

        Self {
            selection_dist_id,
            mu_dist_id,
            data,
        }
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

    pub fn mutate(
        &mut self,
        core_vec: &Vec<Vec<u8>>,
        seq: &mut Vec<u8>,
        selection_dist: &Distribution,
        mu_dist: &Distribution,
    ) -> usize {
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

            let allele_index =
                Self::allele_to_index(value).expect("Sampled non-mutable N allele for mutation");
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
        n_sites
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
        let selection_dist =
            Distribution::new_uniform(0.0, 1.0).expect("failed to create selection distribution");
        let mu_dist =
            Distribution::new_poisson(10.0).expect("failed to create mutation-rate distribution");
        let core_vec: Vec<Vec<u8>> =
            vec![vec![2, 4, 8], vec![1, 4, 8], vec![1, 2, 8], vec![1, 2, 4]];

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
}
