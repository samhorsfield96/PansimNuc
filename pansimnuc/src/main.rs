// TODO allow TEs to have specific multiplative effects governed by a distribution
// TODO plot GFFs with ggGenome to show how the genome evolves over time

mod config;
mod gff;
mod mutation;
mod population;
mod structural;

use crate::mutation::Distribution;
use crate::structural::StructureMutationMap;
use clap::Parser;
use config::Config;
use gff::read_gff_lines;
#[cfg(debug_assertions)]
use gff::write_root_genome_gff;
use itertools::Itertools;
use population::Population;
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::collections::HashMap;

#[derive(Parser, Debug)]
#[command(name = "pansimnuc")]
#[command(about = "Forward simulation of a base at nucleotide level, including structural variation, TE movement and selection", long_about = None)]
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
                for (key, value) in configuration.iter().sorted_by_key(|x| x.0) {
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
        let mut n_threads: usize = n_threads_str
            .parse::<usize>()
            .expect("threads must be an integer.");
        if n_threads < 1 {
            n_threads = 1;
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()
            .expect("Failed to initialise rayon pool");
    }

    // initialise RNG with seed for reproducibility
    let seed_str = configuration
        .get("misc.seed")
        .expect("Require seed for reproducibility. Please specify in config file.");
    let seed: u64 = seed_str.parse::<u64>().expect("seed must be an integer.");
    let mut rng: StdRng = StdRng::seed_from_u64(seed);

    let earlgrey_gff_path: Option<&str> = configuration
        .get("input.earlgrey_gff_file")
        .map(|s| s.as_str());

    if let (Some(gff_path), Some(fasta_path)) = (
        configuration.get("input.gff_file"),
        configuration.get("input.fasta_file"),
    ) {
        match read_gff_lines(&gff_path, &fasta_path, earlgrey_gff_path) {
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

                if let (Some(n_individuals_str), Some(n_generation_str)) = (
                    configuration.get("population.n_individuals"),
                    configuration.get("population.n_generations"),
                ) {
                    let parse_f64 = |key: &str| -> f64 {
                        let value = configuration
                            .get(key)
                            .unwrap_or_else(|| panic!("Missing required config key: {}", key));
                        value
                            .parse::<f64>()
                            .unwrap_or_else(|_| panic!("Config key '{}' must be a float", key))
                    };

                    let parse_usize = |key: &str| -> usize {
                        let value = configuration
                            .get(key)
                            .unwrap_or_else(|| panic!("Missing required config key: {}", key));
                        value
                            .parse::<usize>()
                            .unwrap_or_else(|_| panic!("Config key '{}' must be an integer", key))
                    };

                    let feature_sections = ["exons", "introns", "intergenic", "TE-CUT", "TE-COPY"];

                    // generate distributions to draw mutations from
                    let mut site_mutation_dists: Vec<Distribution> = Vec::new();
                    let mut site_mutation_mus: Vec<Distribution> = Vec::new();
                    let mut structural_dists: Vec<StructureMutationMap> = Vec::new();
					let mut multiplier_dists: Vec<Distribution> = Vec::new();

                    for section in feature_sections {
                        let mutation_rate_key = format!("{}.mutation_rate", section);
                        let duplication_rate_key = format!("{}.duplication_rate", section);
                        let deletion_rate_key = format!("{}.deletion_rate", section);
                        let inversion_rate_key = format!("{}.inversion_rate", section);
                        let max_duplications_key = format!("{}.max_duplications", section);
                        let duplication_insertion_prob_key =
                            format!("{}.duplication_insertion_prob", section);

                        let mutation_rate = parse_f64(&mutation_rate_key);

                        site_mutation_dists.push(
							Distribution::from_selection_config(&configuration, section).unwrap_or_else(|err| {
								panic!(
									"Failed to create selection distribution for section '{}': {}",
									section, err
								)
							}),
						);

                        site_mutation_mus.push(
							Distribution::new_poisson(mutation_rate).unwrap_or_else(|_| {
								panic!(
									"Failed to create mutation-rate distribution for section '{}' from key '{}'",
									section, mutation_rate_key
								)
							}),
						);

                        let max_duplications = configuration
                            .get(&max_duplications_key)
                            .map(|_| parse_usize(&max_duplications_key));

                        structural_dists.push(StructureMutationMap {
                            duplication_rate: parse_f64(&duplication_rate_key),
                            deletion_rate: parse_f64(&deletion_rate_key),
                            inversion_rate: parse_f64(&inversion_rate_key),
                            max_duplications,
                            duplication_insertion_prob: parse_f64(&duplication_insertion_prob_key),
                        });

						// get multiplier just for TE sections
						if section.starts_with("TE") {
							let multiplier_rate_key = format!("{}.multiplier_rate", section);
							let multiplier_scale_key = format!("{}.multiplier_scale", section);
							let multiplier_rate = parse_f64(&multiplier_rate_key);
							let multiplier_scale = parse_f64(&multiplier_scale_key);
							multiplier_dists.push(Distribution::new_gamma(multiplier_rate, multiplier_scale).unwrap());
						} else {
							// set multiplier to 1 for non-TE sections
							multiplier_dists.push(Distribution::new_gamma(1.0, 1.0).unwrap());
						}
                    }

                    // recombination distributions
                    let recombination_rate = parse_f64("population.recombination_rate");
                    let recombination_size_mean = parse_f64("population.recombination_size_mean");
                    let recombination_prob_dist = Distribution::new_poisson(recombination_rate)
                        .expect("Failed to create recombination probability distribution");
                    let recombination_size_dist = Distribution::new_poisson(
                        recombination_size_mean,
                    )
                    .expect("Failed to create recombination distance probability distribution");
                    let recombination_threshold = parse_f64("population.recombination_threshold");

                    let recombination_dists =
                        vec![recombination_prob_dist, recombination_size_dist];

                    // generate initial population
                    let n_individuals: usize = n_individuals_str
                        .parse::<usize>()
                        .expect("n_individuals must be an integer.");
                    println!("Initialising population...");
                    let mut population = Population::new(
                        features,
                        n_individuals,
                        site_mutation_dists,
                        site_mutation_mus,
                        recombination_dists,
                        recombination_threshold,
                        structural_dists,
                        parse_usize("population.max_multiplier_dist"),
						multiplier_dists,
                        &mut rng,
                    );
                    println!("Finished initialising population...");

                    // mutate population
                    let n_generation: usize = n_generation_str
                        .parse::<usize>()
                        .expect("n_generation must be an integer.");

                    for generation in 1..=n_generation {
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

                    if let Some(output_gff) = configuration.get("output.gff_file") {
                        if let Err(err) = population.write_gff(output_gff) {
                            eprintln!("Failed to write final population GFF files: {err}");
                            std::process::exit(1);
                        }
                    }

                    if let Err(err) = population.write_fasta(&output_fasta) {
                        eprintln!("Failed to write final population FASTA files: {err}");
                        std::process::exit(1);
                    }
                }
            }
            Err(err) => {
                eprintln!("Failed to read input files: {err}");
                std::process::exit(1);
            }
        }
    } else {
        eprintln!("Failed to read required input files: input.gff_file, input.fasta_file");
    }
}
