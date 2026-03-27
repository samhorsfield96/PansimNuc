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