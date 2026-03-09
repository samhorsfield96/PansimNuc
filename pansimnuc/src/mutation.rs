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
use rand::Rng;
use rand_distr::Normal;
use rand::distr::Uniform;
use statrs::distribution::Exp;
use rand::rngs::StdRng;

struct DoubleExponential {
    lambda1: f64,
    lambda2: f64,
    prop_lambda1: f64,
    cutoff: f64,
    sample: f64
}

impl DoubleExponential {
    pub fn new(lambda1: f64, lambda2: f64, prop_lambda1: f64, cutoff: f64, rng: &mut StdRng) -> Self {
        let uniform: Uniform<f64> = Uniform::new(0.0, 1.0).unwrap();
        let exponential1: Exp = Exp::new(lambda1).unwrap();
        let exponential2: Exp = Exp::new(lambda2).unwrap();

        let weight: f64 = uniform.sample(&mut rng) as f64;
        let mut sample: f64 = 0.0;

        // higher selected gene
        if weight <= cutoff {
            sample = exponential2.sample(&mut rng);
            //selection_coeffient = 100.0;
        } else {
            sample = exponential1.sample(&mut rng);

            while sample > 1.0 {
                sample = exponential1.sample(&mut rng);
            }

            sample = -1.0 * sample;
        }

        Self { lambda1, lambda2, prop_lambda1, cutoff, sample }
    }
    

}

#[derive(Clone)]
pub enum Distribution {
    Normal { mean: f64, std_dev: f64 },
    Uniform { low: f64, high: f64 },
    Exp { lambda: f64 },
    DoubleExp { lambda1: f64, lambda2: f64, prop_lambda1: f64, cutoff: f64, sample: f64 }
}

impl Distribution {
    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        match self {
            Distribution::Normal { mean, std_dev } => {
                Normal::new(*mean, *std_dev).unwrap().sample(rng)
            }
            Distribution::Uniform { low, high } => {
                Uniform::new(*low, *high).unwrap().sample(rng)
            }
            Distribution::Exp { lambda } => {
                Exp::new(*lambda).unwrap().sample(rng)
            }
        }
    }
}

pub struct MutationMap<V> {
    
    data: [FxHashMap<u32, V>; 4],
}

impl<V> MutationMap<V> {
    pub fn new() -> Self {
        Self { data: std::array::from_fn(|_| FxHashMap::default()) }
    }

    fn insert(&mut self, level: usize, key: u32, value: V) {
        self.data[level].insert(key, value);
    }

    fn get(&self, level: usize, key: u32) -> Option<&V> {
        self.data[level].get(&key)
    }
}