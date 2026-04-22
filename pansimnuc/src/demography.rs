use crate::population::Population;
use crate::config::{self, PopulationSplitConfig};
use crate::tracking::write_tracking_output;
use rand::Rng;
use rand::seq::IteratorRandom;
use std::collections::HashSet;
use rayon::prelude::*;
use std::collections::HashMap;

pub struct MetaPopulation {
    pub populations: Vec<Population>,
    pub population_split_config: PopulationSplitConfig,
    pub n_generations: usize,
    pub recombination_rate: f64,
    pub recombination_size_mean: f64,
    pub site_mutation_mus_vals: Vec<f64>,
}

impl MetaPopulation {
    pub fn new(
        population: Population, 
        population_split_config: PopulationSplitConfig, 
        n_generations: usize,
        recombination_rate: f64,
        recombination_size_mean: f64,
        site_mutation_mus_vals: Vec<f64>,
    ) -> Self {
        MetaPopulation {
            populations: vec![population],
            population_split_config,
            n_generations,
            recombination_rate,
            recombination_size_mean,
            site_mutation_mus_vals,
        }
    }

    fn max_population_id(&self) -> usize {
        self.populations.iter().map(|p| p.id).max().unwrap_or(0)
    }

    fn split_population(&mut self) {
        // pick random population to duplicate
        let mut rng = rand::thread_rng();
        let population_index = rng.gen_range(0..self.populations.len());
        let mut new_population = self.populations[population_index].clone();

        // duplicate population and assign new ID
        new_population.id = self.max_population_id() + 1;

        self.populations.push(new_population);
    }

    fn merge_populations(&mut self) {
        if self.populations.len() < 2 {
            return; // Need at least two populations to merge
        }

        // pick two random populations to merge
        let mut rng = rand::thread_rng();
        let selected_indices = (0..self.populations.len()).choose_multiple(&mut rng, 2);

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

    fn migrate(&mut self) -> usize {
        if self.populations.len() < 2 {
            return 0; // Need at least two populations to migrate
        }

        let mut n_migrations = 0;

        let migration_rate = self.population_split_config.migration_rate;

        let mut rng = rand::thread_rng();
        let population_indices: Vec<usize> = (0..self.populations.len()).collect();

        let mut updated_populations: HashSet<usize> = HashSet::new();

        for i in 0..self.populations.len() {
            let genomes_to_migrate: Vec<(usize, usize, _)> = {
                let source_pop = &self.populations[i];
                let target_indices: Vec<usize> = population_indices.iter().cloned().filter(|&idx| idx != i).collect();

                source_pop.pop.iter()
                    .filter_map(|genome| {
                        if rng.gen_bool(migration_rate) {
                            let target_pop_idx = target_indices[rng.gen_range(0..target_indices.len())];

                            // keep track of updated populations
                            updated_populations.insert(target_pop_idx);

                            n_migrations += 1;

                            let target_pop_length = self.populations[target_pop_idx].pop.len();
                            let target_genome_idx = self.populations[target_pop_idx].pop[rng.gen_range(0..target_pop_length)].genome_id;
                            Some((target_pop_idx, target_genome_idx, genome.clone()))
                        } else {
                            None
                        }
                    })
                    .collect()
            };

            for (target_pop_idx, target_genome_idx, mut migrated_genome) in genomes_to_migrate {
                migrated_genome.genome_id = target_genome_idx; // assign genome ID of target genome to migrated genome to maintain population structure 
                self.populations[target_pop_idx].pop[target_genome_idx] = migrated_genome;
            }
        }

        // update homology maps for all populations after migration, to ensure they are consistent with the new population structure
        for population_idx in updated_populations {
            self.populations[population_idx].update_homology_map();
        }

        n_migrations
    }

    pub fn run_simulation(&mut self, is_tracking: bool, configuration: &HashMap<String, String>) {
        let verbose = self.populations[0].verbose; // assume all populations have same verbose setting

        let mut print_all_generations = false;
        if let Some(print_all_generations_str) = configuration.get("misc.print_all_generations") {
            print_all_generations = print_all_generations_str
                .parse::<bool>()
                .expect("print_all_generations must be a boolean (true/false).");
        }

        for generation in 1..=self.n_generations {
            self.populations.par_iter_mut().for_each(|population| {
                let mut rng = rand::thread_rng();
                
                // mutate at nucleotide level
                let (total_snps, total_indels) = population.mutate();

                // perform intragenome structural mutations
                population.structural_intra_genome();

                // perform intergenome structural mutations, based on number of SNPs
                population.structural_inter_genome(self.recombination_rate, total_snps, self.recombination_size_mean);

                // sample next generation
                let sampled_indices = population.sample_individuals(&mut rng);
                population.next_generation(sampled_indices);

                if generation < self.n_generations {
                    population.update_mu_dists(&self.site_mutation_mus_vals);
                }
            });

            // perform migration between populations
            let n_migrations = self.migrate();

            if verbose {
                println!("{} migration events", n_migrations);
            }

            // perform population splits and merges at specified generations
            let current_gen_size = self.populations.len();

            // determine if population split required at this generation by getting index of current generation in population_split_config.generation_splits, if it exists
            if let Some(split_index) = self.population_split_config.generation_splits.iter().position(|&g| g == generation) {
                let new_gen_size = self.population_split_config.population_splits[split_index];

                if new_gen_size > current_gen_size {
                    // perform population splits until we have the required number of populations

                    if verbose {
                        println!("Splitting populations from {} to {} at generation {}", current_gen_size, new_gen_size, generation);
                    }

                    while self.populations.len() < new_gen_size {
                        self.split_population();
                    }
                } else if new_gen_size < current_gen_size {

                    if verbose {
                        println!("Merging populations from {} to {} at generation {}", current_gen_size, new_gen_size, generation);
                    }

                    // perform population merges until we have the required number of populations
                    while self.populations.len() > new_gen_size {
                        self.merge_populations();
                    }
                }
            }
            
            // print tracking information for this generation if tracking enabled
            if is_tracking {
                if let Some(outdir) = configuration
                    .get("output.outdir") {
                    let output_path = format!("{}/tracking.csv", outdir);
                    write_tracking_output(&output_path, &self);
                };
            }

            if print_all_generations {
                self.write_output(configuration);
            }

            println!("Finished generation {}", generation);
        }
    }

    pub fn write_output(&self, configuration: &HashMap<String, String>) {
        for population in &self.populations {
            if let Some(outdir) = configuration.get("output.outdir") {
                let gff_path = format!("{}/.gff", outdir);
                if let Err(err) = population.write_gff(&gff_path, false) {
                    eprintln!("Failed to write final population GFF files: {err}");
                    std::process::exit(1);
                }

                let fasta_path = format!("{}/.fasta", outdir);
                if let Err(err) = population.write_fasta(&fasta_path, false) {
                    eprintln!("Failed to write final population FASTA files: {err}");
                    std::process::exit(1);
                }
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
                total_exon_length: 0,
                total_intron_length: 0,
                total_intergenic_length: 0,
                total_te_cut_length: 0,
                total_te_copy_length: 0,
                total_tracking_length: 0,
                total_elements: 0,
                total_exon_elements: 0,
                total_intron_elements: 0,
                total_intergenic_elements: 0,
                total_te_cut_elements: 0,
                total_te_copy_elements: 0,
                total_tracking_elements: 0,
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
            indel_dists: Vec::new(),
            structural_mu_dists: Vec::new(),
            recombination_dists: Vec::new(),
            recombination_threshold: 0.0,
            homology_map: Vec::new(),
            feature_map: std::collections::HashMap::new(),
            max_multiplier_dist: 0,
            n_generations: 1,
            verbose: false,
            augment_tracking: false,
            genome_size_penalty_per_bp: 0.0,
            optimal_genome_size: 0,
        }
    }

    #[test]
    fn test_new_initialises_single_population() {
        let population = test_population(7, 2, "p0");
        let split_config = test_split_config(0.1);

        let meta = MetaPopulation::new(population.clone(), split_config.clone(), 10, 0.01, 100.0, vec![0.01]);

        assert_eq!(meta.populations.len(), 1);
        assert_eq!(meta.populations[0].id, 7);
        assert_eq!(meta.population_split_config, split_config);
    }

    #[test]
    fn test_max_population_id_returns_largest_id() {
        let population = test_population(2, 1, "p0");
        let mut meta = MetaPopulation::new(population, test_split_config(0.1), 10, 0.01, 100.0, vec![0.01]);
        meta.populations.push(test_population(9, 1, "p1"));
        meta.populations.push(test_population(4, 1, "p2"));

        assert_eq!(meta.max_population_id(), 9);
    }

    #[test]
    fn test_split_population_adds_population_with_new_id() {
        let mut meta = MetaPopulation::new(test_population(3, 2, "p0"), test_split_config(0.1), 10, 0.01, 100.0, vec![0.01]);

        meta.split_population();

        assert_eq!(meta.populations.len(), 2);
        let ids: Vec<usize> = meta.populations.iter().map(|p| p.id).collect();
        assert!(ids.contains(&3));
        assert!(ids.contains(&4));
    }

    #[test]
    fn test_merge_populations_merges_two_into_one() {
        let mut meta = MetaPopulation::new(test_population(1, 2, "p0"), test_split_config(0.1), 10, 0.01, 100.0, vec![0.01]);
        meta.populations.push(test_population(5, 2, "p1"));

        meta.merge_populations();

        assert_eq!(meta.populations.len(), 1);
        assert_eq!(meta.populations[0].id, 6);
        assert_eq!(meta.populations[0].pop.len(), 2);
    }

    #[test]
    fn test_merge_populations_noop_when_single_population() {
        let mut meta = MetaPopulation::new(test_population(10, 2, "p0"), test_split_config(0.1), 10, 0.01, 100.0, vec![0.01]);

        meta.merge_populations();

        assert_eq!(meta.populations.len(), 1);
        assert_eq!(meta.populations[0].id, 10);
    }

    #[test]
    fn test_migrate_noop_when_rate_zero() {
        let mut meta = MetaPopulation::new(test_population(1, 2, "a"), test_split_config(0.0), 10, 0.01, 100.0, vec![0.01]);
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
        let mut meta = MetaPopulation::new(test_population(1, 2, "solo"), test_split_config(1.0), 10, 0.01, 100.0, vec![0.01]);
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

    #[test]
    fn test_migrate_moves_members_between_populations() {
        let mut meta = MetaPopulation::new(test_population(1, 2, "a"), test_split_config(1.0), 10, 0.01, 100.0, vec![0.01]);
        meta.populations.push(test_population(2, 2, "b"));

        let pop00_before = meta.populations[0].pop[0].identifier.clone();
        let pop01_before = meta.populations[0].pop[1].identifier.clone();
        let pop10_before = meta.populations[1].pop[0].identifier.clone();
        let pop11_before = meta.populations[1].pop[1].identifier.clone();

        meta.migrate();

        let pop00_after = meta.populations[0].pop[0].identifier.clone();
        let pop01_after = meta.populations[0].pop[1].identifier.clone();
        let pop10_after = meta.populations[1].pop[0].identifier.clone();
        let pop11_after = meta.populations[1].pop[1].identifier.clone();

        let moved_from_pop0_to_pop1 = pop10_after == pop00_before || pop11_after == pop01_before || pop10_after == pop01_before || pop11_after == pop00_before;
        let moved_from_pop1_to_pop0 = pop00_after == pop10_before || pop01_after == pop11_before || pop00_after == pop11_before || pop01_after == pop10_before;

        // assert that A genome moves from pop0 to pop1 and B genome moves from pop1 to pop0

        println!("Pop00 before: {}, Pop01 before: {}", pop00_before, pop01_before);
        println!("Pop10 before: {}, Pop11 before: {}", pop10_before, pop11_before);
        println!("Pop00 after: {}, Pop01 after: {}", pop00_after, pop01_after);
        println!("Pop10 after: {}, Pop11 after: {}", pop10_after, pop11_after);

        assert!(
            moved_from_pop0_to_pop1 || moved_from_pop1_to_pop0,
            "Expected at least one genome to move between populations during migration"
        );
    }
}