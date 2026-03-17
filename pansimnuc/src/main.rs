mod gff;
mod config;
mod population;
mod structural;
mod mutation;

use clap::Parser;
use gff::read_gff_lines;
#[cfg(debug_assertions)]
use gff::write_root_genome_gff;
use population::Population;
use config::Config;
use std::collections::HashMap;
use crate::mutation::Distribution;
use rand::rngs::StdRng;
use rand::{SeedableRng};

#[derive(Parser, Debug)]
#[command(name = "pansimnuc")]
#[command(about = "Forward simulation of base genomes at nucleotide level", long_about = None)]
struct Args {
	#[arg(long, help = "Optional path to config file")]
	config: Option<String>,
}

fn main() {
	let args = Args::parse();

	let mut configuration: HashMap<String, String> = HashMap::new();

	// Load config if provided
	if let Some(config_path) = &args.config {
		match Config::from_file(config_path) {
			Ok(config) => {
				println!("Loaded config from: {}", config_path);
				// Flatten config into a single HashMap for easy access
				configuration = config.flatten();
				println!("Configuration values:");
				for (key, value) in &configuration {
					println!("  {} = {}", key, value);
				}
			}
			Err(err) => {
				eprintln!("Failed to read config file: {err}");
				std::process::exit(1);
			}
		}
	}

	// enable multithreading
	if let Some(n_threads_str) = configuration.get("misc.threads") {
		let mut n_threads: usize = n_threads_str.parse::<usize>().expect("threads must be an integer.");
		if n_threads < 1 {
			n_threads = 1;
		}
		rayon::ThreadPoolBuilder::new()
			.num_threads(n_threads)
			.build_global()
			.expect("Failed to initialise rayon pool");
	}

	// initialise RNG with seed for reproducibility
 	let seed_str = configuration.get("misc.seed").expect("Require seed for reproducibility. Please specify in config file.");
	let seed: u64 = seed_str.parse::<u64>().expect("seed must be an integer.");
	let mut rng: StdRng = StdRng::seed_from_u64(seed);

	if let (Some(gff_path), Some(fasta_path), Some(earlgrey_gff_path)) = (
		configuration.get("input.gff_file"),
		configuration.get("input.fasta_file"),
		configuration.get("input.earlgrey_gff_file"),
	) {
		match read_gff_lines(&gff_path, &fasta_path, Some(&earlgrey_gff_path)) {
			Ok(features) => {
				println!("Loaded {} contigs with features", features.len());

				#[cfg(debug_assertions)]
				{
					let root_gff_path = configuration
						.get("output.root_gff_file")
						.cloned()
						.unwrap_or_else(|| "root_genome.debug.gff3".to_string());

					if let Err(err) = write_root_genome_gff(&features, &root_gff_path) {
						eprintln!("Failed to write debug root genome GFF: {err}");
					} else {
						println!("Wrote debug root genome GFF: {}", root_gff_path);
					}
				}

				if let (Some(n_individuals_str), Some(n_generation_str)) = (configuration.get("population.n_individuals"), configuration.get("population.n_generations")) {
					
					// generate distributions to draw mutations from
					// selection distrubtions
        			let exon_dist = Distribution::new_double_exp(1.0, 10.0, 0.5).expect("Failed to create selection distribution for exon features");
        			let intron_dist = Distribution::new_normal(0.0, 1.0).expect("Failed to create selection distribution for intron features");
        			let intergenic_dist = Distribution::new_exp(1.0).expect("Failed to create selection distribution for intergenic features");

					// mutation rate distributions
					let exon_mu = Distribution::new_poisson(0.000001).expect("Failed to create mu dist for exon features");
        			let intron_mu = Distribution::new_poisson(0.00001).expect("Failed to create mu dist for intron features");
        			let intergenic_mu = Distribution::new_poisson(0.0001).expect("Failed to create mu dist for intergenic features");

					// recombination distributions
					let recombination_prob_dist = Distribution::new_poisson(5.0).expect("Failed to create recombination probability distribution");
					let recombination_size_dist = Distribution::new_poisson(1000.0).expect("Failed to create recombination distance probability distribution");
					let recombination_threshold = 0.90;

					let site_mutation_dists = vec![exon_dist, intron_dist, intergenic_dist];
					let site_mutation_mus = vec![exon_mu, intron_mu, intergenic_mu];
					let recombination_dists = vec![recombination_prob_dist, recombination_size_dist];

					// generate initial population
					let n_individuals: usize = n_individuals_str.parse::<usize>().expect("n_individuals must be an integer.");
					println!("Initialising population...");
					let mut population = Population::new(features, n_individuals, site_mutation_dists, site_mutation_mus, recombination_dists, recombination_threshold, &mut rng);
					println!("Finished initialising population...");

					// mutate population
					let n_generation: usize = n_generation_str.parse::<usize>().expect("n_generation must be an integer.");
										
					for generation in 0..n_generation {
						// mutate at nucleotide level
						population.mutate();

						// perform intragenome structural mutations
						population.structural_intra_genome();

						// perform intergenome structural mutations
						population.structural_inter_genome();

						// sample next generation
						let sampled_indices = population.sample_individuals(&mut rng);
						population.next_generation(sampled_indices);
						eprintln!("Finished generation {generation}");
					}
					

					println!("Writing output...");
					let output_fasta = configuration
						.get("output.fasta_file")
						.cloned()
						.unwrap_or_else(|| "final_population.fasta".to_string());

					if let Err(err) = population.write_fasta(&output_fasta) {
						eprintln!("Failed to write final population FASTA files: {err}");
						std::process::exit(1);
					}

					println!("Wrote final population FASTA files with per-genome prefixes based on: {}", output_fasta);
				}
			}
			Err(err) => {
				eprintln!("Failed to read input files: {err}");
				std::process::exit(1);
			}
		}
	} else {
		eprintln!("Failed to read required input files: input.gff_file, input.fasta_file, input.earlgrey_gff_file");
	}
}