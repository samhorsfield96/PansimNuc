mod gff;
mod config;
mod population;
mod mutation;

use clap::Parser;
use gff::read_gff_lines;
use population::Population;
use config::Config;
use std::collections::HashMap;
use crate::mutation::Distribution;
use rand::rngs::StdRng;

#[derive(Parser, Debug)]
#[command(name = "pansimnuc")]
#[command(about = "Forward simulation of base genome at nucleotide level", long_about = None)]
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

	if let (Some(gff_path), Some(fasta_path)) = (configuration.get("input.gff_file"), configuration.get("input.fasta_file")) {
		match read_gff_lines(&gff_path, &fasta_path) {
			Ok(features) => {
				println!("Loaded {} contigs with features", features.len());

				if let (Some(n_individuals_str), Some(n_generation_str)) = (configuration.get("population.n_individuals"), configuration.get("population.n_generations")) {
					
					// generate distributions to draw mutations from
					// placeholder for mutation map, which will be added to each NucElement in the genome
        			let exon_dist = Distribution::new_double_exp(0.0, 1.0, 0.5).expect("Failed to create double exponential distribution for exon features");
        			let intron_dist = Distribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intron features");
        			let intergenic_dist = Distribution::new_uniform(0.0, 1.0).expect("Failed to create uniform distribution for intergenic features");

					let site_mutation_dists = vec![exon_dist, intron_dist, intergenic_dist];

					// generate initial population
					let n_individuals: usize = n_individuals_str.parse::<usize>().expect("n_individuals must be an integer.");
					let mut population = Population::new(features, n_individuals, site_mutation_dists);

					let n_generation: usize = n_generation_str.parse::<usize>().expect("n_generation must be an integer.");
				}
			}
			Err(err) => {
				eprintln!("Failed to read input files: {err}");
				std::process::exit(1);
			}
		}
	} else {
		eprintln!("Failed to read input gff and fasta files");
	}
}