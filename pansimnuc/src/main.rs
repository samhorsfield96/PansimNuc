// TODO plot GFFs with ggGenome to show how the genome evolves over time
// Plot time series data with allelic drift/with selection added
// TODO add selection contraint on genome size to prevent runaway genome expansion

mod config;
mod gff;
mod mutation;
mod population;
mod structural;
mod demography;
mod tracking;
use rayon::prelude::*;

use crate::config::PopulationSplitConfig;
use crate::demography::MetaPopulation;
use crate::mutation::Distribution;
use crate::tracking::write_tracking_header;
use clap::Parser;
use config::Config;
use gff::read_gff_lines;
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
    let mut population_split_config = PopulationSplitConfig::new();
    let mut tracking_regions: Vec<(String, usize, usize)> = Vec::new(); // vector of (contig_id, start, end) tuples for regions to track
    let mut augment_tracking = false;
    let mut genome_size_penalty_per_bp = 0.0;

    // Load config if provided
    let mut verbose = false;
    if let Some(config_path) = &args.config {
        match Config::from_file(config_path) {
            Ok(config) => {
                println!("Loaded config from: {}", config_path);
                // Flatten config into a single HashMap for easy access
                configuration = config.flatten();

                if let Some(verbose_str) = configuration.get("misc.verbose") {
                    verbose = verbose_str
                        .parse::<bool>()
                        .expect("verbose must be a boolean (true/false).");
                }
                if verbose {
                    println!("Verbose mode enabled.");
                }

                if verbose {
                    println!("Configuration values:");
                    for (key, value) in configuration.iter().sorted_by_key(|x| x.0) {
                        println!("  {} = {}", key, value);
                    }
                }

                population_split_config = config.population_split_config().unwrap_or_else(|err| {
                    println!("No population split configuration provided: {err}");
                    PopulationSplitConfig::new()
                });

                tracking_regions = config.tracking_regions().unwrap_or_else(|err| {
                    println!("No tracking regions provided: {err}");
                    Vec::new()
                });

                if let Some(augment_str) = configuration.get("tracking.augmentation") {
                    augment_tracking = augment_str
                        .parse::<bool>()
                        .expect("Tracking augmentation must be a boolean (true/false).");
                }

                if let Some(genome_size_penalty_str) = configuration.get("population.genome_size_penalty_per_bp") {
                    genome_size_penalty_per_bp = genome_size_penalty_str
                        .parse::<f64>()
                        .expect("Genome size penalty must be a float.");
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
            Ok((features, contig_name_to_id)) => {
                println!("Loaded {} contigs with features", features.len());

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

                    let feature_sections = ["exons", "introns", "intergenic", "TE-CUT", "TE-COPY", "tracking"];

                    // generate distributions to draw mutations from
                    let mut site_mutation_dists: Vec<Distribution> = Vec::new();
                    let mut site_mutation_mus_vals: Vec<f64> = Vec::new();
                    let mut site_indel_mus_vals: Vec<f64> = Vec::new();
                    let mut structural_dists: Vec<Vec<Distribution>> = Vec::new();
					let mut multiplier_dists: Vec<Distribution> = Vec::new();

                    for section in feature_sections {
                        let mutation_rate_key = format!("{}.mutation_rate", section);
                        let duplication_rate_key = format!("{}.duplication_rate", section);
                        let deletion_rate_key = format!("{}.deletion_rate", section);
                        let inversion_rate_key = format!("{}.inversion_rate", section);
                        let indel_rate_key = format!("{}.indel_rate", section);

                        site_mutation_dists.push(
							Distribution::from_selection_config(&configuration, section).unwrap_or_else(|err| {
								panic!(
									"Failed to create selection distribution for section '{}': {}",
									section, err
								)
							}),
						);

                        // each position takes two elements of vector
                        site_mutation_mus_vals.push(parse_f64(&mutation_rate_key));
                        site_indel_mus_vals.push(parse_f64(&indel_rate_key));

                        structural_dists.push(vec![
                            Distribution::new_poisson(parse_f64(&duplication_rate_key)).expect("Failed to create duplication distribution"),
                            Distribution::new_poisson(parse_f64(&deletion_rate_key)).expect("Failed to create deletion distribution"),
                            Distribution::new_poisson(parse_f64(&inversion_rate_key)).expect("Failed to create inversion distribution"),
                        ]);

						// get multiplier just for TE sections and tracking
						if section.starts_with("TE") || section == "tracking" {
							let multiplier_rate_key = format!("{}.multiplier_rate", section);
							let multiplier_scale_key = format!("{}.multiplier_scale", section);
							let multiplier_rate = parse_f64(&multiplier_rate_key);
							let multiplier_scale = parse_f64(&multiplier_scale_key);
							multiplier_dists.push(Distribution::new_gamma(multiplier_rate, multiplier_scale).expect("Failed to create multiplier distribution"));
						} else {
							// set multiplier to 1 for non-TE sections
							multiplier_dists.push(Distribution::new_gamma(1.0, 1.0).expect("Failed to create multiplier distribution"));
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

                    let n_generations: usize = n_generation_str
                        .parse::<usize>()
                        .expect("n_generation must be an integer.");
                    
                    println!("Initialising population...");
                    // genome size penalty is 1 - penalty for multiplication
                    genome_size_penalty_per_bp = 1.0 - genome_size_penalty_per_bp;
                    let mut population = Population::new(
                        features,
                        n_individuals,
                        site_mutation_dists,
                        &site_mutation_mus_vals,
                        &site_indel_mus_vals,
                        recombination_dists,
                        recombination_threshold,
                        structural_dists,
                        parse_usize("population.max_multiplier_dist"),
						multiplier_dists,
                        n_generations,
                        &mut rng,
                        verbose,
                        &contig_name_to_id,
                        &tracking_regions,
                        augment_tracking,
                        genome_size_penalty_per_bp,
                    );

                    let mut is_tracking = false;
                    // generate root genome
                    if let Some(outdir) = configuration
                        .get("output.outdir") {
                        let gff_path = format!("{}/.gff", outdir);
                        population.write_gff(&gff_path, true).unwrap_or_else(|err| {
                            panic!("Failed to write root genome GFF file: {err}");
                        });
                        let fasta_path = format!("{}/.fasta", outdir);
                        population.write_fasta(&fasta_path, true).unwrap_or_else(|err| {
                            panic!("Failed to write root genome FASTA file: {err}");
                        });

                        // add tracking regions to population struct so we can track mutations in these regions over time

                        if !tracking_regions.is_empty() {
                            is_tracking = true;
                            let tracking_out_path = format!("{}/tracking.csv", outdir);
                            write_tracking_header(&tracking_out_path);
                        }
                    }

                    // generate metapopulation with different mutation distributions
                    let mut metapopopulation = MetaPopulation::new(
                        population, 
                        population_split_config, 
                        n_generations, 
                        recombination_rate, 
                        recombination_size_mean, 
                        site_mutation_mus_vals,
                    );

                    println!("Finished initialising population...");

                    // run simulation
                    metapopopulation.run_simulation(is_tracking, &configuration);

                    println!("Writing output...");
                    metapopopulation.write_output(&configuration);
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
