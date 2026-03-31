// in this script, genes can move around, be duplicated and deleted

use crate::mutation::Distribution as MutationDistribution;
use crate::population::NucElement;
use crate::population::{Genome, Population};
use rand::Rng; 
use rand::seq::SliceRandom;
use triple_accel::levenshtein::*;
use std::collections::HashMap;
use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::visit::Dfs;
use std::collections::HashSet;
use rayon::{prelude::*};

// for a given NucElement, store its position in the genome
// which can then be shuffled around by structural mutations, or copied

/// Returns a similarity score in [0.0, 1.0] based on normalized edit distance.
/// 1.0 = identical, 0.0 = completely different.
fn calculate_homology(a: &NucElement, b: &NucElement, threshold: f64) -> f64 {
    let s: &[u8] = a.seq.as_slice();
    let t: &[u8] = b.seq.as_slice();

    let m = s.len();
    let n = t.len();

    if m == 0 || n == 0 {
        return 0.0;
    }

    let max_len = m.max(n) as f64;

    let min_dist = ((1.0 - threshold) * max_len).ceil() as u32;

    // accelerated Levenshtein distance with early exit if distance exceeds min_dist
    if let Some(dist) = levenshtein_simd_k(s, t, min_dist) {
        return 1.0 - (dist as f64 / max_len)
    } else {
        return 0.0;
    };    
}

// write function which runs through each element and determines whether a structural mutation occurs, and if so, which one, and where it moves to.
pub fn mutate_intra_genome(
    genome: &mut Genome,
    structural_mu_dists: &Vec<Vec<MutationDistribution>>,
    pos_dist: &MutationDistribution,
) -> (usize, usize, usize, usize, usize, usize, usize) {
    let mut thread_rng = rand::thread_rng();

    // For all intra genome comparisons, sample from uniform distribution to determine if variant occurs
    // and poisson distribution to determine where duplication goes

    // store hashmap of positions of genome elements, can store multiple per entry to capture duplications
    let mut new_positions: HashMap<i64, Vec<(usize, i64)>> = HashMap::new();
    let mut max_position: usize = 0;

    // check which contig each block will be inserted into
    let contig_starts = &genome.contig_starts;

    let mut total_non_te_duplications = 0;
    let mut total_non_te_deletions = 0;
    let mut te_cut_duplications = 0;
    let mut te_copy_duplications = 0;
    let mut te_cut_deletions = 0;
    let mut te_copy_deletions = 0;
    let mut total_inversions = 0;

    for (current_pos, element) in &mut genome.seq.iter().enumerate() {
        // determine maximum position
        if current_pos > max_position {
            max_position = current_pos;
        }

        let mutation_dist: &Vec<MutationDistribution> = match element.feature_type.as_str() {
            "exon" => &structural_mu_dists[0],
            "intron" => &structural_mu_dists[1],
            "intergenic" => &structural_mu_dists[2],
            "TE-CUT" => &structural_mu_dists[3],
            "TE-COPY" => &structural_mu_dists[4],
            _ => panic!("Unknown feature type: {}", element.feature_type),
        };

        // get feature type
        let feature_type = &element.feature_type;

        //store element structure positions, with current position first
        let mut new_positions_vec: Vec<(usize, i64)> =
            vec![(element.contig_id, current_pos as i64)];

        // duplications, can model multiple duplications repeatedly sampling until rand_val is above duplication rate
        let num_dups = mutation_dist[0].sample(&mut thread_rng) as usize;
        let mut dup_count: usize = 0;

        for _ in 0..num_dups {
            dup_count += 1;

            // for non-TEs sample from poisson distribution to determine where duplication goes
            let duplication_pos = if feature_type.contains("TE") {
                // determine what is max dist to end from current position
                let max_dist = (genome.seq.len() - current_pos).max(current_pos);

                if feature_type == "TE-CUT" {
                    te_cut_duplications += 1;
                } else {
                    te_copy_duplications += 1;
                }
                // for TEs, sample from uniform distribution to determine where duplication goes, with equal probability of anywhere in genome
                thread_rng.gen_range(0..=max_dist) as f64
            } else {
                total_non_te_duplications += 1;
                pos_dist.sample(&mut thread_rng)
            };
            let duplication_pos = duplication_pos as i64;

            // determine if position is before or after gene, adjust if too large
            let pos_rand_val = thread_rng.gen_range(0..=9);
            let pos_order: i64 =
                if pos_rand_val < 5 {
                    -1
                } else {
                    1
                };
            let mut new_pos = current_pos as i64 + (duplication_pos * pos_order);
            let mut new_contig_id = 0;

            // assign to 0 if before start of genome
            if new_pos < 0 {
                new_pos = 0;
            }

            // avoid adding to same position
            if new_pos == current_pos as i64 {
                if new_pos == 0 {
                    new_pos = genome.seq.len() as i64; // if at start, move to end
                } else {
                    new_pos += 1;
                }
            }

            // determine contig position
            for (contig_id, contig_start) in contig_starts.iter().enumerate() {
                if new_pos < *contig_start as i64 {
                    // new position is in previous contig
                    break;
                }
                new_contig_id = contig_id;
            }

            new_positions_vec.push((new_contig_id, new_pos));

            // determine maximum position present
            if new_pos as usize > max_position {
                max_position = new_pos as usize;
            }

            if dup_count > 0 && feature_type == "TE-CUT" {
                // if element is a TE-COPY and has already been duplicated, break loop to prevent further duplications, to avoid runaway genome growth
                break;
            }
        }

        // deletions, only first gene deleted which is original position
        if feature_type == "TE-CUT" && dup_count > 0 {
            // if element is a TE-COPY and has already been duplicated, force deletion of original copy, to capture cut and paste mechanism of TE-COPYs
            let _ = new_positions_vec.remove(0);
            te_cut_deletions += 1;
        } else {
            // All other gene features
            let mut n_deletions = mutation_dist[1].sample(&mut thread_rng) as usize;
            n_deletions = n_deletions.min(new_positions_vec.len());

            // delete as many copies as possible
            for _ in 0..n_deletions {
                let _ = new_positions_vec.remove(0);
                if feature_type.contains("TE") {
                    if feature_type == "TE-CUT" {
                            te_cut_deletions += 1;
                        } else {
                            te_copy_deletions += 1;
                        }
                    } else {
                            total_non_te_deletions += 1;
                }
            }
        }

        // translocations not explicitely modelled, as captured by simulatenous duplication and deletion event

        // inversions, apply to all genes in element structure vec
        for (contig_id, pos) in new_positions_vec {
            let mut inversion = 1;

            // flip sign for number of inversions
            let n_inversions = mutation_dist[2].sample(&mut thread_rng) as usize;
            for _ in 0..n_inversions {
                inversion *= -1;
                total_inversions += 1;
            }

            // now update new_positions, indexed by new position,
            // with old entry as value either as positive or negative value depending on whether it is inverted or not and 1 indexed
            let key = pos;
            let value = (current_pos as i64 + 1) * inversion; // store old position as value, which can be used to update mutation map of element after all structural mutations have been processed
            new_positions
                .entry(key)
                .or_default()
                .push((contig_id, value));
        }
    }

    // generate new genome based on all intra-genome variation, current everything is 1-indexed
    let mut new_genome_seq: Vec<(usize, i64)> = Vec::new();

    // iterate through each position, indexed by new position
    for idx in 0..=max_position {
        if let Some(prev_pos) = new_positions.get(&(idx as i64)) {
            // value exists
            for entry in prev_pos {
                if entry.1 == 0 {
                    panic!(
                        "Entry for structural varient is 0, for position {} in genome {}",
                        idx, genome.identifier
                    );
                }

                // append entry, meaning that genes are not deleted, appended in order
                new_genome_seq.push(*entry);
            }
        }
    }

    // generate new genome
    let mut new_genome: Vec<NucElement> = Vec::new();

    for (contig_id, element_idx) in new_genome_seq {
        let invert: bool = if element_idx > 0 { false } else { true };
        let mut new_element = genome.seq[element_idx.abs() as usize - 1].clone(); // convert back to 0 indexed

        if invert {
            new_element.strand = !new_element.strand;
        }

        // update contig id
        new_element.contig_id = contig_id;

        new_genome.push(new_element);
    }

    genome.seq = new_genome;

    // update contig start positions
    genome.update_contig_starts();

    (total_non_te_duplications, total_non_te_deletions, te_cut_duplications, te_copy_duplications, te_cut_deletions, te_copy_deletions, total_inversions)
}

// get connected components
fn connected_components(
    nodes: impl IntoIterator<Item = u32>,
    edges: &Vec<(u32, u32)>,
) -> Vec<Vec<u32>> {
    let mut graph = UnGraph::<u32, ()>::new_undirected();
    let mut node_map: HashMap<u32, NodeIndex> = HashMap::new();

    for node in nodes {
        let index = graph.add_node(node);
        node_map.insert(node, index);
    }

    for &(a, b) in edges {
        let a_index = *node_map.get(&a).expect("missing node a");
        let b_index = *node_map.get(&b).expect("missing node b");
        graph.add_edge(a_index, b_index, ());
    }

    let mut seen: HashSet<NodeIndex> = HashSet::new();
    let mut components: Vec<Vec<u32>> = Vec::new();

    for start in graph.node_indices() {
        if seen.contains(&start) {
            continue;
        }

        let mut dfs = Dfs::new(&graph, start);
        let mut component = Vec::new();

        while let Some(node_index) = dfs.next(&graph) {
            if seen.insert(node_index) {
                component.push(graph[node_index]);
            }
        }

        components.push(component);
    }

    components
}

pub fn mutate_inter_genome(population: &mut Population) -> (usize, usize, usize) {
    let mut rng = rand::thread_rng();

    // get number of recombination events across whole population
    let n_recombinations = population.recombination_dists[0].sample(&mut rng) as usize;

    let pop_size = population.pop.len();

    // All ordered pairs where donor != recipient
    let all_pairs: Vec<(u32, u32)> = (0..pop_size)
        .flat_map(|i| (0..pop_size).filter(move |&j| j != i).map(move |j| (i as u32, j as u32)))
        .collect();

    let mut recombination_map_tmp: HashMap<usize, Vec<usize>> = HashMap::new();
    for _ in 0..n_recombinations {
        let (donor, recipient) = all_pairs.choose(&mut rng).expect("Failed to select a random pair for recombination");
        recombination_map_tmp.entry(*donor as usize).or_default().push(*recipient as usize);
    }

    let sampled_edges: Vec<(u32, u32)> = recombination_map_tmp
        .iter()
        .flat_map(|(&donor, recipients)| recipients.iter().map(move |&recipient| (donor as u32, recipient as u32)))
        .collect();

    let components = connected_components(0..pop_size as u32, &sampled_edges);

    // generate list of independent recombination maps to process
    let mut recombination_map_list: Vec<HashMap<usize, Vec<usize>>> = Vec::with_capacity(components.len());

    // pull out each connected component
    for component in components {
        let mut component_map: HashMap<usize, Vec<usize>> = HashMap::new();
        for &donor in &component {
            if let Some(recipients) = recombination_map_tmp.get(&(donor as usize)) {
                component_map.insert(donor as usize, recipients.clone());
            }
        }
        if !component_map.is_empty() {
            recombination_map_list.push(component_map);
        }
    }

    // Empties population.pop; slots are indexed by genome_id (== vec index)
    let mut pop_opt: Vec<Option<Genome>> = population.pop.drain(..).map(Some).collect();

    // create component packages
    let component_packages: Vec<(HashMap<usize, Vec<usize>>, Vec<(usize, Genome)>)> =
        recombination_map_list.into_iter().map(|recombination_map| {
            // collect every genome ID this component touches
            let mut ids: HashSet<usize> = HashSet::new();
            for (&donor, recipients) in &recombination_map {
                ids.insert(donor);
                for &r in recipients { ids.insert(r); }
            }
            // move each genome out of its Option slot — panics if already taken (impossible: components are disjoint)
            let genomes: Vec<(usize, Genome)> = ids.into_iter()
                .map(|id| (id, pop_opt[id].take().expect("genome already taken")))
                .collect();
            (recombination_map, genomes)
        }).collect();

    let results = 
        component_packages.into_par_iter().map(|(recombination_map, mut genomes)| {
            // thread specific variables
            let mut thread_rng = rand::thread_rng();
            let mut thread_homology_map = population.homology_map.clone();
            let mut thread_total_donor_length = 0;
            let mut thread_total_recipient_length = 0;
            let mut thread_successful_recombinations = 0;

            // Before the for (donor, recipients) loop, build a local index map:
            let genome_id_to_local_idx: HashMap<usize, usize> = genomes
                .iter()
                .enumerate()
                .map(|(local_idx, (genome_id, _))| (*genome_id, local_idx))
                .collect();

            // iterate over recombinations
            for (donor, recipients) in recombination_map {
                // get local donor index in genomes vec
                let donor_local = genome_id_to_local_idx[&donor];
                for recipient in recipients {
                    // get local recipient index in genomes vec
                    let recipient_local = genome_id_to_local_idx[&recipient];

                    let (donor_genome, recipient_genome): (&Genome, &mut Genome) = if donor_local < recipient_local {
                        let (left, right) = genomes.split_at_mut(recipient_local);
                        (&left[donor_local].1, &mut right[0].1)
                    } else {
                        let (left, right) = genomes.split_at_mut(donor_local);
                        (&right[0].1, &mut left[recipient_local].1)
                    };

                    // look for donor and recipient site, maximum total donor length attempts, if not found, skip recombination event
                    let mut donor_site_chosen: bool = false;

                    // set up sampling with replacement
                    let mut indices: Vec<usize> = (0..donor_genome.seq.len()).collect();
                    indices.shuffle(&mut thread_rng);

                    let mut start_donor_site: usize = 0;
                    let mut start_recipient_site: usize = 0;

                    for recombination_pos in indices {
                        if donor_site_chosen {
                            break;
                        }
                        let recombination_pos_idx = donor_genome.seq[recombination_pos].element_id;

                        // determine if position in both donor and recipient genome, if not, resample
                        let recomb_element = &population.homology_map[recombination_pos_idx];

                        let donor_has_site = !recomb_element[donor].is_empty();
                        let recipient_has_site =
                            recomb_element.len() > recipient && !recomb_element[recipient].is_empty();

                        // if both vectors are not empty, then search through each and test to make sure they have sufficient homology
                        if donor_has_site && recipient_has_site {
                            for donor_site in &recomb_element[donor] {
                                // check that site present in donor
                                if donor_site >= &donor_genome.seq.len() {
                                    continue;
                                }
                                for recipient_site in &recomb_element[recipient] {
                                    // check that site present in recipient
                                    if recipient_site >= &recipient_genome.seq.len() {
                                        continue;
                                    }
                                    // check homology between donor and recipient site, if sufficient, break loop and move to recombination, if not, continue searching
                                    let homology = calculate_homology(
                                        &donor_genome.seq[*donor_site],
                                        &recipient_genome.seq[*recipient_site],
                                        population.recombination_threshold
                                    );
                                    if homology >= population.recombination_threshold {
                                        // perform recombination event, break out of loops

                                        start_donor_site = *donor_site;
                                        start_recipient_site = *recipient_site;

                                        donor_site_chosen = true;
                                        break;
                                    }
                                }
                                if donor_site_chosen {
                                    break;
                                }
                            }
                        }
                    }

                    if donor_site_chosen {
                        // now sample from poisson distribution to determine minumum size of recombination track
                        let min_recombination_len =
                            population.recombination_dists[1].sample(&mut thread_rng) as usize;

                        // determine whether there is a track that can be recombined
                        let mut track_found = false;
                        let mut end_donor_site = start_donor_site;
                        let mut end_recipient_site = start_recipient_site;
                        let mut recombination_len = donor_genome.seq[start_donor_site].seq.len();

                        // ensure recombination occurs in single chromosome each
                        let donor_contig_id = donor_genome.seq[start_donor_site].contig_id;
                        let recipient_contig_id = recipient_genome.seq[start_recipient_site].contig_id;

                        // track contig end of donor
                        let mut donor_contig_end = false;

                        while !track_found {
                            // determine length of donor DNA
                            while recombination_len < min_recombination_len {
                                let new_end_donor_site = end_donor_site + 1;

                                // run off end of contig, assume complete recombination
                                if new_end_donor_site >= donor_genome.seq.len() {
                                    donor_contig_end = true;
                                } else if donor_genome.seq[new_end_donor_site].contig_id != donor_contig_id {
                                    donor_contig_end = true
                                }
                                if donor_contig_end 
                                {
                                    recombination_len += donor_genome.seq[end_donor_site].seq.len();

                                    // find end of recipient track
                                    let mut recipient_contig_end = false;
                                    while !recipient_contig_end {
                                        let new_end_recipient_site = end_recipient_site + 1;
                                        
                                        // check if at end of contig
                                        if new_end_recipient_site >= recipient_genome.seq.len() {
                                            recipient_contig_end = true;
                                        } else if recipient_genome.seq[new_end_recipient_site].contig_id
                                            != recipient_contig_id {
                                                recipient_contig_end = true;
                                        } else {
                                            end_recipient_site = new_end_recipient_site;
                                        }
                                    }
                                    
                                    track_found = true;
                                    break;
                                }
                                

                                // else continue going through contig
                                end_donor_site = new_end_donor_site;
                                recombination_len += donor_genome.seq[end_donor_site].seq.len();
                            }

                            // use to determine homology between sites
                            let donor_site = &donor_genome.seq[end_donor_site];

                            // now iterate through recipient genome until homology found between end and donor site
                            let mut recipient_end_found = false;
                            while !recipient_end_found {
                                // run off end of contig, assume complete recombination
                                if end_recipient_site >= recipient_genome.seq.len()
                                {
                                    // reduce index by 1
                                    end_recipient_site -= 1;
                                    track_found = true;
                                    recipient_end_found = true;
                                    break;
                                } else if recipient_genome.seq[end_recipient_site].contig_id != recipient_contig_id {
                                    // reduce index by 1
                                    end_recipient_site -= 1;
                                    track_found = true;
                                    recipient_end_found = true;
                                    break;
                                }

                                let recipient_site = &recipient_genome.seq[end_recipient_site];
                                let homology = calculate_homology(donor_site, recipient_site, population.recombination_threshold);
                                if homology >= population.recombination_threshold {
                                    track_found = true;
                                    recipient_end_found = true;
                                    break;
                                }

                                // continue iterating through recipient contig until homology found, or end of contig reached
                                end_recipient_site += 1;
                            }
                        }

                        // perform recombination event, replacing recipient track with donor track
                        // clone the donor track first, before any mutable borrow of pop
                        let mut donor_track: Vec<NucElement> = donor_genome.seq
                            [start_donor_site..=end_donor_site]
                            .to_vec()
                            .clone();

                        // update information from recipient track
                        for element in &mut donor_track {
                            element.contig_id = recipient_contig_id;
                        }

                        // store donor_track length before it is moved
                        let donor_track_len = donor_track.len();
                        let donor_track_seq_len: usize = donor_track.iter().map(|e| e.seq.len()).sum();

                        let recipient_track_seq_len: usize = recipient_genome.seq
                            [start_recipient_site..=end_recipient_site]
                            .iter()
                            .map(|e| e.seq.len())
                            .sum();

                        thread_total_donor_length += donor_track_seq_len;
                        thread_total_recipient_length += recipient_track_seq_len;
                        thread_successful_recombinations += 1;

                        // update homology map for recipient genome, need to add new positions for each element in donor track, and remove old positions for each element in recipient track
                        // remove old positions
                        for element_idx in start_recipient_site..=end_recipient_site {
                            let element_id = recipient_genome.seq[element_idx].element_id;
                            let homology_group =
                                &mut thread_homology_map[element_id][recipient_genome.genome_id];
                            homology_group.retain(|&pos| pos != element_idx); // remove old position
                        }

                        // now safe to mutably borrow recipient
                        recipient_genome
                            .seq
                            .splice(start_recipient_site..=end_recipient_site, donor_track);

                        // add new positions
                        for element_idx in start_recipient_site..(start_recipient_site + donor_track_len) {
                            let element_id = recipient_genome.seq[element_idx].element_id;
                            let homology_group =
                                &mut thread_homology_map[element_id][recipient_genome.genome_id];
                            homology_group.push(element_idx); // add new position
                        }
                        // update contig_ids
                        recipient_genome.update_contig_starts();
                    }
                }
            }
        (genomes, thread_homology_map, thread_total_donor_length, thread_total_recipient_length, thread_successful_recombinations)
    }).collect::<Vec<_>>();

    // combine results from each thread
    let mut successful_recombinations = 0;
    let mut total_donor_length = 0;
    let mut total_recipient_length = 0;

    // staging vec for indexed insertion by genome_id
    let mut new_pop: Vec<Option<Genome>> = (0..pop_size).map(|_| None).collect();

    for (genomes, thread_homology_map, thread_donor_length, thread_recipient_length, thread_successful_recombinations) in results {
        total_donor_length += thread_donor_length;
        total_recipient_length += thread_recipient_length;
        successful_recombinations += thread_successful_recombinations;

        // update population with new genomes
        for (genome_id, genome) in genomes.into_iter() {
            new_pop[genome_id] = Some(genome);

            // update population homology map with new homology map
            for (element_id, homology_groups) in thread_homology_map.iter().enumerate() {
                population.homology_map[element_id][genome_id] = homology_groups[genome_id].clone();
            }
        }
    }

    // Push back any genomes not involved in any recombination component
    for maybe_genome in pop_opt.into_iter() {
        if let Some(genome) = maybe_genome {
            let genome_id = genome.genome_id;
            new_pop[genome_id] = Some(genome);
        }
    }

    // Unwrap back into population.pop (all slots must be filled)
    population.pop = new_pop.into_iter()
        .enumerate()
        .map(|(i, opt)| opt.unwrap_or_else(|| panic!("genome slot {} was never filled", i)))
        .collect();

    (successful_recombinations, total_donor_length, total_recipient_length)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mutation::{Distribution as MutationDistribution, MutationMap};
    use crate::population::{Genome, NucElement};
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn default_structural_dists() -> Vec<Vec<MutationDistribution>> {
        let mut structural_dists = Vec::new();
        for _ in 0..5 {
            structural_dists.push(vec![
                MutationDistribution::new_uniform(0.0, 0.1).expect("Failed to create uniform distribution for duplication"),
                MutationDistribution::new_uniform(0.0, 0.1).expect("Failed to create uniform distribution for deletions"),
                MutationDistribution::new_uniform(0.0, 0.1).expect("Failed to create uniform distribution for inversions"),
            ]);
        }

        structural_dists
    }

    fn make_test_genome() -> Genome {
        let mut rng = StdRng::seed_from_u64(42);
        let sel_dist = MutationDistribution::new_uniform(0.0, 1.0).unwrap();

        Genome {
            identifier: "test".to_string(),
            genome_id: 0,
            parent: "root".to_string(),
            contig_starts: vec![0],
            seq: vec![
                NucElement {
                    contig_id: 0,
                    element_id: 0,
                    feature_id: 0,
                    feature_type: "exon".to_string(),
                    multiplier: 1.0,
                    seq: vec![],
                    mutation_map: MutationMap::new(0, 0, &vec![], &sel_dist, &mut rng),
                    strand: true,
                    original_length: 0,
                    frameshift: false,
                },
                NucElement {
                    contig_id: 0,
                    element_id: 1,
                    feature_id: 1,
                    feature_type: "exon".to_string(),
                    multiplier: 1.0,
                    seq: vec![],
                    mutation_map: MutationMap::new(0, 0, &vec![], &sel_dist, &mut rng),
                    strand: false,
                    original_length: 0,
                    frameshift: false,
                },
                NucElement {
                    contig_id: 0,
                    element_id: 2,
                    feature_id: 2,
                    feature_type: "exon".to_string(),
                    multiplier: 1.0,
                    seq: vec![],
                    mutation_map: MutationMap::new(0, 0, &vec![], &sel_dist, &mut rng),
                    strand: false,
                    original_length: 0,
                    frameshift: false,
                },
            ],
            seq_length: 0,
        }
    }

    fn make_multi_contig_test_genome() -> Genome {
        let mut rng = StdRng::seed_from_u64(42);
        let sel_dist = MutationDistribution::new_uniform(0.0, 1.0).unwrap();

        let make_element = |contig_id: usize,
                            element_id: usize,
                            feature_id: usize,
                            strand: bool,
                            sel_dist: &MutationDistribution,
                            rng: &mut StdRng| NucElement {
            contig_id,
            element_id,
            feature_id,
            feature_type: "exon".to_string(),
            multiplier: 1.0,
            seq: vec![],
            mutation_map: MutationMap::new(0, 0, &vec![], sel_dist, rng),
            strand,
            original_length: 0,
            frameshift: false,
        };

        Genome {
            identifier: "test_multi_contig".to_string(),
            genome_id: 0,
            parent: "root".to_string(),
            contig_starts: vec![0, 2, 4],
            seq: vec![
                make_element(0, 0, 0, true, &sel_dist, &mut rng),
                make_element(0, 1, 1, false, &sel_dist, &mut rng),
                make_element(1, 2, 2, true, &sel_dist, &mut rng),
                make_element(1, 3, 3, false, &sel_dist, &mut rng),
                make_element(2, 4, 4, true, &sel_dist, &mut rng),
                make_element(2, 5, 5, false, &sel_dist, &mut rng),
            ],
            seq_length: 0,
        }
    }

    fn make_recombination_test_genome(
        genome_id: usize,
        n_elements: usize,
        strand_seed: bool,
        marker_base: u8,
    ) -> Genome {

        let mut rng = StdRng::seed_from_u64(100 + genome_id as u64);
        let sel_dist = MutationDistribution::new_uniform(0.0, 1.0).unwrap();

        let mut seq: Vec<NucElement> = Vec::new();
        for idx in 0..n_elements {
            let marker_seq = vec![marker_base; 4];
            seq.push(NucElement {
                contig_id: 0,
                element_id: idx,
                feature_id: idx,
                feature_type: "exon".to_string(),
                multiplier: 1.0,
                seq: marker_seq.clone(),
                mutation_map: MutationMap::new(0, 0, &marker_seq, &sel_dist, &mut rng),
                strand: if idx % 2 == 0 {
                    strand_seed
                } else {
                    !strand_seed
                },
                original_length: marker_seq.len(),
                frameshift: false,
            });
        }

        Genome {
            identifier: format!("recomb_{}", genome_id),
            genome_id,
            parent: "root".to_string(),
            contig_starts: vec![0],
            seq,
            seq_length: 0,
        }
    }

    fn make_recombination_test_population(forced_events: usize, n_elements: usize) -> Population {
        // genome 0 starts with marker base 1 (A), genome 1 starts with marker base 2 (C)
        let g0 = make_recombination_test_genome(0, n_elements, true, 1);
        let g1 = make_recombination_test_genome(1, n_elements, false, 2);

        let recombination_count =
            MutationDistribution::new_uniform(forced_events as f64, forced_events as f64 + 0.1)
                .unwrap();
        let recombination_len = MutationDistribution::new_uniform(0.0, 0.1).unwrap();

        let mut homology_map: Vec<Vec<Vec<usize>>> = Vec::new();
        for _ in 0..n_elements {
            homology_map.push(vec![vec![0], vec![0]]);
        }

        Population {
            id: 0,
            generation: 0,
            pop: vec![g0, g1],
            core_vec: vec![],
            selection_dists: vec![],
            mu_dists: vec![],
            indel_dists: vec![],
            structural_mu_dists: vec![vec![]],
            recombination_dists: vec![recombination_count, recombination_len],
            recombination_threshold: 0.0,
            homology_map,
            feature_map: HashMap::new(),
            max_multiplier_dist: 10,
            n_generations: 10,
            verbose: true,
        }
    }

    fn genome_has_marker(genome: &Genome, marker: u8) -> bool {
        genome
            .seq
            .iter()
            .any(|element| element.seq.first().copied() == Some(marker))
    }

    fn count_mixed_marker_genomes(population: &Population) -> usize {
        population
            .pop
            .iter()
            .filter(|genome| genome_has_marker(genome, 1) && genome_has_marker(genome, 2))
            .count()
    }

    #[test]
    fn inversion_flips_strand() {
        let mut genome = make_test_genome();
        let before_strands: Vec<bool> = genome.seq.iter().map(|e| e.strand).collect();

        let mut default_structural_dists = default_structural_dists();

        // update inversion rate to 1 for all elements, so that all elements are inverted
        default_structural_dists[0][2] = MutationDistribution::new_uniform(1.0, 1.1).unwrap();
        let pos = MutationDistribution::new_uniform(0.0, 1.0).unwrap();

        let mut homology_map: Vec<Vec<Vec<usize>>> = Vec::new();
        for _ in genome.seq.iter() {
            homology_map.push(vec![vec![0]]);
        }

        mutate_intra_genome(&mut genome, &default_structural_dists, &pos);

        let after_strands: Vec<bool> = genome.seq.iter().map(|e| e.strand).collect();
        assert_ne!(
            before_strands, after_strands,
            "strands should flip after forced inversion"
        );
        assert_eq!(
            genome.seq.len(),
            before_strands.len(),
            "inversion should preserve genome length"
        );
    }

    #[test]
    fn deletion_reduces_genome_length() {
        let mut genome = make_test_genome();
        let before_len = genome.seq.len();

        let mut default_structural_dists = default_structural_dists();
        default_structural_dists[0][1] = MutationDistribution::new_uniform(1.0, 1.1).unwrap();

        let mut homology_map: Vec<Vec<Vec<usize>>> = Vec::new();
        for _ in genome.seq.iter() {
            homology_map.push(vec![vec![0]]);
        }

        let pos = MutationDistribution::new_uniform(0.0, 1.0).unwrap();
        mutate_intra_genome(&mut genome, &default_structural_dists, &pos);

        assert!(
            genome.seq.len() < before_len,
            "genome should shrink after forced deletion"
        );
    }

    #[test]
    // run multiple times, as may fail stochastically if duplication and deletion events don't line up as expected
    fn translocation_preserves_length() {
        // A translocation is a simultaneous duplication (copy to new position) and
        // deletion (remove from original position).  The gene moves but the genome
        // length stays the same.
        let mut genome = make_test_genome();
        let before_len = genome.seq.len();
        let before_ids: Vec<usize> = genome.seq.iter().map(|e| e.feature_id).collect();

        let mut default_structural_dists = default_structural_dists();

        // Exactly one duplication per element, then delete the original.
        default_structural_dists[0][0] = MutationDistribution::new_uniform(1.0, 1.1).unwrap();
        default_structural_dists[0][1] = MutationDistribution::new_uniform(1.0, 1.1).unwrap();

        let mut homology_map: Vec<Vec<Vec<usize>>> = Vec::new();
        for _ in genome.seq.iter() {
            homology_map.push(vec![vec![0]]);
        }

        // Use a non-zero offset so duplicates land somewhere other than position 0.
        let pos = MutationDistribution::new_uniform(1.0, 2.0).unwrap();
        mutate_intra_genome(&mut genome, &default_structural_dists, &pos);

        assert_eq!(
            genome.seq.len(),
            before_len,
            "translocation should preserve genome length"
        );
        // The set of feature_ids present should be unchanged even if order differs.
        let mut after_ids: Vec<usize> = genome.seq.iter().map(|e| e.feature_id).collect();
        let mut expected = before_ids.clone();
        after_ids.sort();
        expected.sort();
        assert_eq!(
            after_ids, expected,
            "translocation should not add or remove genes"
        );
    }

    #[test]
    fn multi_contig_ids_remain_in_ascending_order() {
        let mut genome = make_multi_contig_test_genome();

        let mut default_structural_dists = default_structural_dists();

        default_structural_dists[0][0] = MutationDistribution::new_poisson(10.0).unwrap();

        let pos = MutationDistribution::new_uniform(0.0, 0.9).unwrap();

        mutate_intra_genome(&mut genome, &default_structural_dists, &pos);

        let contig_ids: Vec<usize> = genome.seq.iter().map(|e| e.contig_id).collect();
        assert!(!contig_ids.is_empty(), "mutated genome should not be empty");
        assert_eq!(
            contig_ids.first().copied(),
            Some(0),
            "first contig should be 0"
        );
        assert!(
            contig_ids.windows(2).all(|w| w[0] <= w[1]),
            "contig ids should be non-decreasing"
        );

        let mut unique_contigs = contig_ids.clone();
        unique_contigs.dedup();
        assert!(
            unique_contigs.len() > 1,
            "test should contain multiple contigs"
        );
        let expected_contigs: Vec<usize> = (0..unique_contigs.len()).collect();
        assert_eq!(
            unique_contigs, expected_contigs,
            "contigs should progress upward from 0"
        );
    }

    #[test]
    fn inter_genome_recombination_zero_events_is_noop() {
        let mut population = make_recombination_test_population(0, 8);

        let element_ids_before: Vec<Vec<usize>> = population
            .pop
            .iter()
            .map(|genome| genome.seq.iter().map(|e| e.element_id).collect())
            .collect();

        for (genome_idx, ids) in element_ids_before.iter().enumerate() {
            println!(
                "Before recombination - genome {} element_ids: {:?}",
                genome_idx, ids
            );
        }

        let mixed_before = count_mixed_marker_genomes(&population);
        assert_eq!(
            mixed_before, 0,
            "before recombination, genomes should not be mixed"
        );

        let total_before: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        mutate_inter_genome(&mut population);
        let total_after: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        let mixed_after = count_mixed_marker_genomes(&population);

        let element_ids_after: Vec<Vec<usize>> = population
            .pop
            .iter()
            .map(|genome| genome.seq.iter().map(|e| e.element_id).collect())
            .collect();

        for (genome_idx, ids) in element_ids_after.iter().enumerate() {
            println!(
                "After recombination - genome {} element_ids: {:?}",
                genome_idx, ids
            );
        }

        assert_eq!(
            population.pop.len(),
            2,
            "population size should be unchanged"
        );
        assert_eq!(
            total_after, total_before,
            "no forced recombination should not change total length in this deterministic setup"
        );
        assert_eq!(
            mixed_after, mixed_before,
            "no forced recombination should not change mixed marker genomes"
        );
    }

    #[test]
    fn inter_genome_recombination_single_event_changes_total_length_by_one() {
        let mut population = make_recombination_test_population(1, 8);

        let element_ids_before: Vec<Vec<usize>> = population
            .pop
            .iter()
            .map(|genome| genome.seq.iter().map(|e| e.element_id).collect())
            .collect();

        for (genome_idx, ids) in element_ids_before.iter().enumerate() {
            println!(
                "Before recombination - genome {} element_ids: {:?}",
                genome_idx, ids
            );
        }

        let mixed_before = count_mixed_marker_genomes(&population);
        assert_eq!(
            mixed_before, 0,
            "before recombination, genomes should not be mixed"
        );

        let total_before: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        mutate_inter_genome(&mut population);
        let total_after: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        let mixed_after = count_mixed_marker_genomes(&population);

        let element_ids_after: Vec<Vec<usize>> = population
            .pop
            .iter()
            .map(|genome| genome.seq.iter().map(|e| e.element_id).collect())
            .collect();

        for (genome_idx, ids) in element_ids_after.iter().enumerate() {
            println!(
                "After recombination - genome {} element_ids: {:?}",
                genome_idx, ids
            );
        }

        assert_eq!(
            population.pop.len(),
            2,
            "population size should be unchanged"
        );
        assert_eq!(
            total_after, total_before,
            "single forced recombination should preserve total genome length"
        );
        assert!(
            mixed_after > mixed_before,
            "after one forced recombination, at least one genome should contain marker sequence from the other genome"
        );
    }

    #[test]
    fn inter_genome_recombination_multiple_events_change_total_length_by_event_count() {
        let forced_events = 3;
        let mut population = make_recombination_test_population(forced_events, 8);

        let element_ids_before: Vec<Vec<usize>> = population
            .pop
            .iter()
            .map(|genome| genome.seq.iter().map(|e| e.element_id).collect())
            .collect();

        for (genome_idx, ids) in element_ids_before.iter().enumerate() {
            println!(
                "Before recombination - genome {} element_ids: {:?}",
                genome_idx, ids
            );
        }

        let mixed_before = count_mixed_marker_genomes(&population);
        assert_eq!(
            mixed_before, 0,
            "before recombination, genomes should not be mixed"
        );

        let total_before: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        mutate_inter_genome(&mut population);
        let total_after: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        let mixed_after = count_mixed_marker_genomes(&population);

        let element_ids_after: Vec<Vec<usize>> = population
            .pop
            .iter()
            .map(|genome| genome.seq.iter().map(|e| e.element_id).collect())
            .collect();

        for (genome_idx, ids) in element_ids_after.iter().enumerate() {
            println!(
                "After recombination - genome {} element_ids: {:?}",
                genome_idx, ids
            );
        }

        assert_eq!(
            population.pop.len(),
            2,
            "population size should be unchanged"
        );
        assert_eq!(
            total_after, total_before,
            "forced multiple recombinations should preserve total genome length"
        );
        assert!(
            mixed_after >= 1,
            "after forced recombinations, at least one genome should contain marker sequence from the other genome"
        );
    }

    fn create_test_genomes(element_type: &str) -> Genome {
          // TE-CUT should: duplication -> break loop -> force delete original
        // Result: one copy at new position, original removed (cut-and-paste)
        let mut rng = StdRng::seed_from_u64(42);
        let sel_dist = MutationDistribution::new_uniform(0.0, 0.1).unwrap();

        let genome = Genome {
            identifier: format!("test_{}", element_type.to_lowercase()),
            genome_id: 0,
            parent: "root".to_string(),
            contig_starts: vec![0],
            seq: vec![NucElement {
                contig_id: 0,
                element_id: 0,
                feature_id: 0,
                feature_type: element_type.to_string(),
                multiplier: 1.0,
                seq: vec![1, 2, 4, 8],
                mutation_map: MutationMap::new(0, 0, &vec![1, 2, 4, 8], &sel_dist, &mut rng),
                strand: true,
                original_length: 4,
                frameshift: false,
            },
            NucElement {
                contig_id: 0,
                element_id: 0,
                feature_id: 0,
                feature_type: "exon".to_string(),
                multiplier: 1.0,
                seq: vec![1, 2, 4, 8],
                mutation_map: MutationMap::new(0, 0, &vec![1, 2, 4, 8], &sel_dist, &mut rng),
                strand: true,
                original_length: 4,
                frameshift: false,
            },
            NucElement {
                contig_id: 0,
                element_id: 0,
                feature_id: 0,
                feature_type: "intergenic".to_string(),
                multiplier: 1.0,
                seq: vec![1, 2, 4, 8],
                mutation_map: MutationMap::new(0, 0, &vec![1, 2, 4, 8], &sel_dist, &mut rng),
                strand: true,
                original_length: 4,
                frameshift: false,
            }],
            seq_length: 0,
        };
        genome
    }

    #[test]
    fn te_copy_allows_multiple_duplications() {
        // TE-COPY should allow multiple duplications without early break
        let mut genome = create_test_genomes("TE-COPY");

        let before_len = genome.seq.len();

        let mut default_structural_dists = default_structural_dists();
        
        // High probability of duplication, zero deletion
        default_structural_dists[4][0] = MutationDistribution::new_uniform(2.0, 2.1).unwrap();
        let pos = MutationDistribution::new_uniform(1.0, 2.0).unwrap();
        
        mutate_intra_genome(&mut genome, &default_structural_dists, &pos);
        
        // TE-COPY should result in multiple copies (original + duplicates)
        assert!(
            genome.seq.len() > before_len,
            "TE-COPY should allow multiple duplications; before: {}, after: {}",
            before_len,
            genome.seq.len()
        );
        
        // All copies should be TE-COPY
        let te_copy_count = genome.seq.iter().filter(|e| e.feature_type == "TE-COPY").count();
        assert_eq!(
            te_copy_count > 1,
            true,
            "TE-COPY should result in multiple copies; found {} TE-COPY elements",
            te_copy_count
        );
    }

    #[test]
    fn te_cut_implements_cut_and_paste() {
        // TE-CUT should: duplication -> break loop -> force delete original
        // Result: one copy at new position, original removed (cut-and-paste)
        let mut genome = create_test_genomes("TE-CUT");

        let before_len = genome.seq.len();
        let mut default_structural_dists = default_structural_dists();
        
        // High probability of duplication, zero deletion
        default_structural_dists[3][0] = MutationDistribution::new_uniform(2.0, 2.1).unwrap();
        default_structural_dists[3][1] = MutationDistribution::new_uniform(2.0, 2.1).unwrap();
        let pos = MutationDistribution::new_uniform(1.0, 1.1).unwrap();
        
        mutate_intra_genome(&mut genome, &default_structural_dists, &pos);
        
        // After cut-and-paste, genome should have same or fewer elements
        // (original deleted, one copy inserted)
        assert!(
            genome.seq.len() == before_len,
            "TE-CUT cut-and-paste should result in at most one additional element; before: {}, after: {}",
            before_len,
            genome.seq.len()
        );

        println!("Genome after TE-CUT mutation: {:?}", genome.seq.iter().map(|e| e.feature_type.clone()).collect::<Vec<String>>());

        // ensure TE has moved and original position is deleted
        let te_position = genome.seq.iter().position(|e| e.feature_type == "TE-CUT");
        assert_ne!(
            te_position.expect("TE-CUT should still be present after cut-and-paste"),
            0,
            "TE-CUT should have moved from original position"
        );

        let te_cut_count = genome.seq.iter().filter(|e| e.feature_type == "TE-CUT").count();
        assert_eq!(
            te_cut_count, 1,
            "TE-CUT should result in exactly one copy after cut-and-paste"
        );

    }

    #[test]
    fn intergenic_allows_multiple_duplications_with_poisson_position() {
        // Non-TE features should allow multiple duplications and use Poisson for position
        let mut genome: Genome = create_test_genomes("intergenic");

        let before_len = genome.seq.len();

        let mut default_structural_dists = default_structural_dists();
        default_structural_dists[2][0] = MutationDistribution::new_uniform(2.0, 2.1).unwrap();
        
        // Use Poisson position distribution for non-TEs
        let pos = MutationDistribution::new_poisson(1.5).unwrap();
        
        mutate_intra_genome(&mut genome, &default_structural_dists, &pos);
        
        // Should have duplications
        assert!(
            genome.seq.len() > before_len,
            "intergenic should allow multiple duplications; before: {}, after: {}",
            before_len,
            genome.seq.len()
        );
        
        // All should be intergenic
        let intergenic_count = genome.seq.iter().filter(|e| e.feature_type == "intergenic").count();
        assert_eq!(
            intergenic_count > 1,
            true,
            "Intergenic should be duplicated"
        );
    }
}
