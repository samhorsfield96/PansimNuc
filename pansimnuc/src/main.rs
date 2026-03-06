mod gff;
mod config;

use clap::Parser;
use gff::read_gff_lines;
use config::Config;
use std::collections::HashMap;

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