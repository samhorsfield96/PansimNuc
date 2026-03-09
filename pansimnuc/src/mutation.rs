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
use rand_distr::{Distribution as _, Normal, Uniform, Exponential, Laplace};

#[derive(Clone)]
pub enum Distribution {
    Normal { mean: f64, std_dev: f64 },
    Uniform { low: f64, high: f64 },
    Exponential { lambda: f64 },
    Laplace { loc: f64, scale: f64 },
}

impl Distribution {
    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        match self {
            Distribution::Normal { mean, std_dev } => {
                Normal::new(*mean, *std_dev).unwrap().sample(rng)
            }
            Distribution::Uniform { low, high } => {
                Uniform::new(*low, *high).sample(rng)
            }
            Distribution::Exponential { lambda } => {
                Exponential::new(*lambda).unwrap().sample(rng)
            }
            Distribution::Laplace { loc, scale } => {
                Laplace::new(*loc, *scale).unwrap().sample(rng)
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