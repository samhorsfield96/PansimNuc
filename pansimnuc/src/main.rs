mod gff;

use clap::Parser;
use gff::read_gff_lines;

#[derive(Parser, Debug)]
#[command(name = "pansimnuc")]
#[command(about = "Extract feature sequences from GFF + FASTA", long_about = None)]
struct Args {
	#[arg(long, help = "Path to GFF3 file")]
	gff: String,

	#[arg(long, help = "Path to FASTA file")]
	fasta: String,
}

fn main() {
	let args = Args::parse();

	match read_gff_lines(&args.gff, &args.fasta) {
		Ok(features) => {
			println!("Loaded {} features", features.len());
		}
		Err(err) => {
			eprintln!("Failed to read input files: {err}");
			std::process::exit(1);
		}
	}
}