mod gff;
mod config;

use clap::Parser;
use gff::read_gff_lines;
use config::Config;

#[derive(Parser, Debug)]
#[command(name = "pansimnuc")]
#[command(about = "Extract feature sequences from GFF + FASTA", long_about = None)]
struct Args {
	#[arg(long, help = "Path to GFF3 file")]
	gff: String,

	#[arg(long, help = "Path to FASTA file")]
	fasta: String,

	#[arg(long, help = "Optional path to config file")]
	config: Option<String>,
}

fn main() {
	let args = Args::parse();

	// Load config if provided
	if let Some(config_path) = &args.config {
		match Config::from_file(config_path) {
			Ok(config) => {
				println!("Loaded config from: {}", config_path);
				for section in config.sections.keys() {
					println!("  [{}]", section);
					for key in config.keys_in_section(section) {
						if let Some(value) = config.get(section, &key) {
							println!("    {} = {}", key, value);
						}
					}
				}
			}
			Err(err) => {
				eprintln!("Failed to read config file: {err}");
				std::process::exit(1);
			}
		}
	}

	match read_gff_lines(&args.gff, &args.fasta) {
		Ok(features) => {
			println!("Loaded {} contigs with features", features.len());
		}
		Err(err) => {
			eprintln!("Failed to read input files: {err}");
			std::process::exit(1);
		}
	}

	
}