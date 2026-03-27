use crate::population::Population;
use crate::config::PopulationSplitConfig;
use rand::Rng;
use rand::seq::IteratorRandom;

pub struct MetaPopulation {
    pub populations: Vec<Population>,
    pub population_split_config: PopulationSplitConfig,
}

impl MetaPopulation {
    pub fn new(population: Population, population_split_config: PopulationSplitConfig) -> Self {
        MetaPopulation {
            populations: vec![population],
            population_split_config,
        }
    }

    fn max_population_id(&self) -> usize {
        self.populations.iter().map(|p| p.id).max().unwrap_or(0)
    }

    pub fn split_population(&mut self) {
        // pick random population to duplicate
        let mut rng = rand::thread_rng();
        let population_index = rng.gen_range(0..self.populations.len());
        let mut new_population = self.populations[population_index].clone();

        // duplicate population and assign new ID
        new_population.id = self.max_population_id() + 1;

        self.populations.push(new_population);
    }

    pub fn merge_populations(&mut self) {
        if self.populations.len() < 2 {
            return; // Need at least two populations to merge
        }

        // pick two random populations to merge
        let mut rng = rand::thread_rng();
        let selected_indices = (0..self.populations.len()).choose_multiple(&mut rng, 2);

        let pop1_id = self.populations[selected_indices[0]].id;
        let pop2_id = self.populations[selected_indices[1]].id;
        let pop1 = &self.populations[selected_indices[0]];
        let pop2 = &self.populations[selected_indices[1]];

        // merge two populations, sampling by replacement to maintain population size
        // create vector of tuples of population ids and genome ids
        let mut merged_pop = Vec::new();
        for genome in pop1.pop.iter() {
            merged_pop.push((pop1.id, genome.genome_id));
        }
        for genome in pop2.pop.iter() {
            merged_pop.push((pop2.id, genome.genome_id));
        }

        // sample with replacement to create new population
        let mut new_pop_genomes = Vec::new();
        for new_genome_id in 0..pop1.pop.len() {
            let idx = rng.gen_range(0..merged_pop.len());
            let (pop_id, genome_id) = merged_pop[idx];
            let mut genome = if pop_id == pop1.id {
                pop1.pop[genome_id].clone()
            } else {
                pop2.pop[genome_id].clone()
            };
            genome.genome_id = new_genome_id; // assign new genome ID
            new_pop_genomes.push(genome);
        }

        let mut merged_population = pop1.clone();
        merged_population.pop = new_pop_genomes;
        merged_population.id = self.max_population_id() + 1;

        // remove old populations and add merged population, remove from back to front to avoid index issues
        if selected_indices[0] < selected_indices[1] {
            self.populations.remove(selected_indices[1]);
            self.populations.remove(selected_indices[0]);
        } else {
            self.populations.remove(selected_indices[0]);
            self.populations.remove(selected_indices[1]);
        }

        self.populations.push(merged_population);
    }

    pub fn migrate(&mut self) {
        if self.populations.len() < 2 {
            return; // Need at least two populations to migrate
        }

        let migration_rate = self.population_split_config.migration_rate;

        let mut rng = rand::thread_rng();
        let population_indices: Vec<usize> = (0..self.populations.len()).collect();

        for i in 0..self.populations.len() {
            let genomes_to_migrate: Vec<(usize, usize, _)> = {
                let source_pop = &self.populations[i];
                let target_indices: Vec<usize> = population_indices.iter().cloned().filter(|&idx| idx != i).collect();

                source_pop.pop.iter()
                    .filter_map(|genome| {
                        if rng.gen_bool(migration_rate) {
                            let target_pop_idx = target_indices[rng.gen_range(0..target_indices.len())];
                            let target_pop_length = self.populations[target_pop_idx].pop.len();
                            let target_genome_idx = self.populations[target_pop_idx].pop[rng.gen_range(0..target_pop_length)].genome_id;
                            Some((target_pop_idx, target_genome_idx, genome.clone()))
                        } else {
                            None
                        }
                    })
                    .collect()
            };

            for (target_pop_idx, target_genome_idx, migrated_genome) in genomes_to_migrate {
                self.populations[target_pop_idx].pop[target_genome_idx] = migrated_genome;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::population::Genome;

    fn test_split_config(migration_rate: f64) -> PopulationSplitConfig {
        PopulationSplitConfig {
            population_splits: vec![1],
            generation_splits: vec![1],
            migration_rate,
        }
    }

    fn test_genomes(n: usize, prefix: &str) -> Vec<Genome> {
        (0..n)
            .map(|idx| Genome {
                identifier: format!("{}_{}", prefix, idx),
                genome_id: idx,
                contig_starts: Vec::new(),
                parent: "root".to_string(),
                seq: Vec::new(),
                seq_length: 0,
            })
            .collect()
    }

    fn test_population(id: usize, n_genomes: usize, prefix: &str) -> Population {
        Population {
            id,
            generation: 0,
            pop: test_genomes(n_genomes, prefix),
            core_vec: Vec::new(),
            selection_dists: Vec::new(),
            mu_dists: Vec::new(),
            structural_mu_dists: Vec::new(),
            recombination_dists: Vec::new(),
            recombination_threshold: 0.0,
            homology_map: Vec::new(),
            feature_map: std::collections::HashMap::new(),
            max_multiplier_dist: 0,
            n_generations: 1,
            verbose: false,
        }
    }

    #[test]
    fn test_new_initialises_single_population() {
        let population = test_population(7, 2, "p0");
        let split_config = test_split_config(0.1);

        let meta = MetaPopulation::new(population.clone(), split_config.clone());

        assert_eq!(meta.populations.len(), 1);
        assert_eq!(meta.populations[0].id, 7);
        assert_eq!(meta.population_split_config, split_config);
    }

    #[test]
    fn test_max_population_id_returns_largest_id() {
        let population = test_population(2, 1, "p0");
        let mut meta = MetaPopulation::new(population, test_split_config(0.1));
        meta.populations.push(test_population(9, 1, "p1"));
        meta.populations.push(test_population(4, 1, "p2"));

        assert_eq!(meta.max_population_id(), 9);
    }

    #[test]
    fn test_split_population_adds_population_with_new_id() {
        let mut meta = MetaPopulation::new(test_population(3, 2, "p0"), test_split_config(0.1));

        meta.split_population();

        assert_eq!(meta.populations.len(), 2);
        let ids: Vec<usize> = meta.populations.iter().map(|p| p.id).collect();
        assert!(ids.contains(&3));
        assert!(ids.contains(&4));
    }

    #[test]
    fn test_merge_populations_merges_two_into_one() {
        let mut meta = MetaPopulation::new(test_population(1, 2, "p0"), test_split_config(0.1));
        meta.populations.push(test_population(5, 2, "p1"));

        meta.merge_populations();

        assert_eq!(meta.populations.len(), 1);
        assert_eq!(meta.populations[0].id, 6);
        assert_eq!(meta.populations[0].pop.len(), 2);
    }

    #[test]
    fn test_merge_populations_noop_when_single_population() {
        let mut meta = MetaPopulation::new(test_population(10, 2, "p0"), test_split_config(0.1));

        meta.merge_populations();

        assert_eq!(meta.populations.len(), 1);
        assert_eq!(meta.populations[0].id, 10);
    }

    #[test]
    fn test_migrate_noop_when_rate_zero() {
        let mut meta = MetaPopulation::new(test_population(1, 2, "a"), test_split_config(0.0));
        meta.populations.push(test_population(2, 2, "b"));

        let before: Vec<Vec<String>> = meta
            .populations
            .iter()
            .map(|p| p.pop.iter().map(|g| g.identifier.clone()).collect())
            .collect();

        meta.migrate();

        let after: Vec<Vec<String>> = meta
            .populations
            .iter()
            .map(|p| p.pop.iter().map(|g| g.identifier.clone()).collect())
            .collect();

        assert_eq!(before, after);
    }

    #[test]
    fn test_migrate_noop_when_single_population() {
        let mut meta = MetaPopulation::new(test_population(1, 2, "solo"), test_split_config(1.0));
        let before: Vec<String> = meta.populations[0]
            .pop
            .iter()
            .map(|g| g.identifier.clone())
            .collect();

        meta.migrate();

        let after: Vec<String> = meta.populations[0]
            .pop
            .iter()
            .map(|g| g.identifier.clone())
            .collect();
        assert_eq!(before, after);
    }
}