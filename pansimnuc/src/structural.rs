// in this script, genes can move around, be duplicated and deleted

use crate::mutation::Distribution as MutationDistribution;
use crate::population::NucElement;
use crate::population::HomologyPositions;
use crate::population::{Genome, Population};
use rand::Rng; 
use rand::seq::SliceRandom;
use triple_accel::levenshtein::*;
use std::borrow::Cow;
use std::collections::HashMap;
use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::visit::Dfs;
use std::collections::HashSet;
use rayon::{prelude::*};

// for a given NucElement, store its position in the genome
// which can then be shuffled around by structural mutations, or copied

/// Returns a similarity score in [0.0, 1.0] based on normalized edit distance.
/// 1.0 = identical, 0.0 = completely different.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base {
            1 => 8,   // A -> T
            2 => 4,   // C -> G
            4 => 2,   // G -> C
            8 => 1,   // T -> A
            16 => 16, // N -> N
            _ => panic!("Allele code must be one-hot (1, 2, 4, 8, 16); got {}", base),
        })
        .collect()
}

fn calculate_homology(a: &NucElement, b: &NucElement, threshold: f64) -> f64 {
    let s: &[u8] = a.seq.as_slice();
    let t: Cow<[u8]> = if a.strand == b.strand {
        Cow::Borrowed(b.seq.as_slice())
    } else {
        Cow::Owned(reverse_complement(b.seq.as_slice()))
    };

    let m = s.len();
    let n = t.len();

    if m == 0 || n == 0 {
        return 0.0;
    }

    let max_len = m.max(n) as f64;

    let min_dist = ((1.0 - threshold) * max_len).ceil() as u32;

    // accelerated Levenshtein distance with early exit if distance exceeds min_dist
    if let Some(dist) = levenshtein_simd_k(s, t.as_ref(), min_dist) {
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
    augment_tracking: bool,
) -> (usize, usize, usize, usize, usize, usize, usize) {
    let mut thread_rng = rand::thread_rng();

    // For all intra genome comparisons, sample from uniform distribution to determine if variant occurs
    // and poisson distribution to determine where duplication goes

    // store hashmap of positions of genome elements, can store multiple per entry to capture duplications
    let mut new_positions: HashMap<i64, Vec<(usize, i64)>> = HashMap::new();

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
        let mut mutation_dist: &Vec<MutationDistribution> = match element.feature_type.as_ref() {
            "exon" => &structural_mu_dists[0],
            "intron" => &structural_mu_dists[1],
            "intergenic" => &structural_mu_dists[2],
            "TE-CUT" => &structural_mu_dists[3],
            "TE-COPY" => &structural_mu_dists[4],
            _ => panic!("Unknown feature type: {}", element.feature_type),
        };

        // override for tracked elements
        if element.tracked && augment_tracking {
            mutation_dist = &structural_mu_dists[5];
        }

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

            let genome_len = genome.seq.len() as i64;
            let mut new_pos = if feature_type.contains("TE") {
                if feature_type.as_ref() == "TE-CUT" {
                    te_cut_duplications += 1;
                } else {
                    te_copy_duplications += 1;
                }
                // For TEs, sample a uniform absolute insertion position across the genome.
                thread_rng.gen_range(0..genome_len)
            } else {
                total_non_te_duplications += 1;
                // For non-TEs, sample a displacement around the current position and
                // wrap around contig boundaries to avoid start/end clamping bias.
                let displacement = pos_dist.sample(&mut thread_rng) as i64;
                let pos_order: i64 = if thread_rng.gen_bool(0.5) { -1 } else { 1 };
                current_pos as i64 + (displacement * pos_order)
            };

            // Wrap into valid genome index range [0, genome_len - 1].
            new_pos = new_pos.rem_euclid(genome_len);
            let mut new_contig_id = 0;

            // determine contig position
            for (contig_id, contig_start) in contig_starts.iter().enumerate() {
                if new_pos < *contig_start as i64 {
                    // new position is in previous contig
                    break;
                }
                new_contig_id = contig_id;
            }

            new_positions_vec.push((new_contig_id, new_pos));

            if dup_count > 0 && feature_type.as_ref() == "TE-CUT" {
                // if element is a TE-CUT and has already been duplicated, break loop to mimic cut and paste mechanism
                break;
            }
        }

        // deletions, only first gene deleted which is original position
        if feature_type.as_ref() == "TE-CUT" && dup_count > 0 {
            // if element is a TE-CUT and has already been duplicated, force deletion of original copy, to capture cut and paste mechanism of TE-CUTs
            let _ = new_positions_vec.remove(0);
            //te_cut_deletions += 1;
        } else {
            // All other gene features
            let mut n_deletions = mutation_dist[1].sample(&mut thread_rng) as usize;
            n_deletions = n_deletions.min(new_positions_vec.len());

            // delete as many copies as possible
            for _ in 0..n_deletions {
                let _ = new_positions_vec.remove(0);
                if feature_type.contains("TE") {
                    if feature_type.as_ref() == "TE-CUT" {
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

    // sort actual keys instead of iterating a sparse 0..=max_position range,
    // which could be enormous when duplications land far from the current position
    let mut sorted_positions: Vec<i64> = new_positions.keys().copied().collect();
    sorted_positions.sort_unstable();

    for idx in sorted_positions {
        for entry in &new_positions[&idx] {
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

    // generate new genome
    // count how many times each source index is referenced so we can move
    // single-reference elements instead of cloning them, avoiding a full
    // duplicate of genome.seq in memory for the common (no-duplication) case
    let mut ref_counts: Vec<usize> = vec![0usize; genome.seq.len()];
    for &(_, element_idx) in &new_genome_seq {
        ref_counts[element_idx.abs() as usize - 1] += 1;
    }

    let mut source: Vec<Option<NucElement>> = genome.seq.drain(..).map(Some).collect();
    let mut new_genome: Vec<NucElement> = Vec::with_capacity(new_genome_seq.len());

    for (contig_id, element_idx) in new_genome_seq {
        let invert = element_idx < 0;
        let src_idx = element_idx.abs() as usize - 1;

        ref_counts[src_idx] -= 1;
        let mut new_element = if ref_counts[src_idx] == 0 {
            // last (or only) reference — move instead of clone
            source[src_idx].take().unwrap()
        } else {
            source[src_idx].as_ref().unwrap().clone()
        };

        if invert {
            new_element.strand = !new_element.strand;
        }
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

pub fn mutate_inter_genome(population: &mut Population, bidirectional: bool) -> (usize, usize, usize) {
    let mut rng = rand::thread_rng();

    // get number of recombination events across whole population
    let n_recombinations = population.recombination_dists[0].sample(&mut rng) as usize;

    let pop_size = population.pop.len();

    if n_recombinations == 0 || pop_size < 2 {
        return (0, 0, 0);
    }

    let mut recombination_map_tmp: HashMap<usize, Vec<usize>> = HashMap::new();
    let mut active_nodes: HashSet<u32> = HashSet::new();
    for _ in 0..n_recombinations {
        // Sample donor/recipient directly to avoid materializing O(pop^2) pair lists.
        let donor = rng.gen_range(0..pop_size);
        let mut recipient = rng.gen_range(0..(pop_size - 1));
        if recipient >= donor {
            recipient += 1;
        }

        recombination_map_tmp
            .entry(donor)
            .or_default()
            .push(recipient);
        active_nodes.insert(donor as u32);
        active_nodes.insert(recipient as u32);
    }

    let sampled_edges: Vec<(u32, u32)> = recombination_map_tmp
        .iter()
        .flat_map(|(&donor, recipients)| recipients.iter().map(move |&recipient| (donor as u32, recipient as u32)))
        .collect();

    if active_nodes.is_empty() {
        return (0, 0, 0);
    }

    let components = connected_components(active_nodes.into_iter(), &sampled_edges);

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
            let mut thread_homology_updates: HashMap<(usize, usize), HomologyPositions> = HashMap::new();
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

                    let (donor_genome, recipient_genome): (&mut Genome, &mut Genome) = if donor_local < recipient_local {
                        let (left, right) = genomes.split_at_mut(recipient_local);
                        (&mut left[donor_local].1, &mut right[0].1)
                    } else {
                        let (left, right) = genomes.split_at_mut(donor_local);
                        (&mut right[0].1, &mut left[recipient_local].1)
                    };

                    // look for donor and recipient site, maximum total donor length attempts, if not found, skip recombination event
                    let mut donor_site_chosen: bool = false;

                    // now sample from poisson distribution to determine minumum size of recombination track
                    let min_recombination_len =
                            population.recombination_dists[1].sample(&mut thread_rng) as usize;
                    
                    // set up sampling with replacement
                    let mut indices: Vec<usize> = (0..donor_genome.seq.len()).collect();
                    indices.shuffle(&mut thread_rng);

                    let mut start_donor_site: usize = 0;
                    let mut start_recipient_site: usize = 0;

                    for recombination_pos in indices {
                        if donor_site_chosen {
                            break;
                        }
                        let element = &donor_genome.seq[recombination_pos];
                        let recombination_pos_idx = element.element_id;

                        // determine if recombination position is usable.
                        let element_contig_length = donor_genome.contig_lengths[element.contig_id];
                        let element_pos = element.feature_pos;

                        // skip if too short for recombination event, only if element is not first in contig, otherwise just recombine the whole chromosome
                        if element_contig_length - element_pos < min_recombination_len && element_pos > 0 {
                            continue;
                        }

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

                        // hold recipient track in for bidirectional recombination
                        let mut recipient_track: Vec<NucElement> = vec![];
                        let mut recipient_track_seq_len: usize = 0;


                        if bidirectional {
                            recipient_track = recipient_genome.seq
                                [start_recipient_site..=end_recipient_site]
                                .to_vec()
                                .clone();
                            recipient_track_seq_len = recipient_track.iter()
                                .map(|e| e.seq.len())
                                .sum();
                        } else {
                            recipient_track_seq_len = recipient_genome.seq
                                [start_recipient_site..=end_recipient_site]
                                .iter()
                                .map(|e| e.seq.len())
                                .sum();
                        }

                        // store donor_track length before it is moved
                        let donor_track_len = donor_track.len();
                        let donor_track_seq_len: usize = donor_track.iter().map(|e| e.seq.len()).sum();

                        // determine recipient track length
                        let recipient_track_len = recipient_track.len();

                        thread_total_donor_length += donor_track_seq_len;
                        thread_total_recipient_length += recipient_track_seq_len;
                        thread_successful_recombinations += 1;

                        // update homology map for recipient genome, need to add new positions for each element in donor track, and remove old positions for each element in recipient track
                        // remove old positions in recipient site
                        for element_idx in start_recipient_site..=end_recipient_site {
                            let element_id = recipient_genome.seq[element_idx].element_id;
                            let homology_group = thread_homology_updates
                                .entry((element_id, recipient_genome.genome_id))
                                .or_insert_with(|| population.homology_map[element_id][recipient_genome.genome_id].clone());
                            homology_group.retain(|pos| *pos != element_idx); // remove old position
                        }

                        // now safe to mutably borrow recipient and update recipient
                        recipient_genome
                            .seq
                            .splice(start_recipient_site..=end_recipient_site, donor_track);

                        // add new positions
                        for element_idx in start_recipient_site..(start_recipient_site + donor_track_len) {
                            let element_id = recipient_genome.seq[element_idx].element_id;
                            let homology_group = thread_homology_updates
                                .entry((element_id, recipient_genome.genome_id))
                                .or_insert_with(|| population.homology_map[element_id][recipient_genome.genome_id].clone());
                            homology_group.push(element_idx); // add new position
                        }
                        // update contig_ids
                        recipient_genome.update_contig_starts();

                        // do the same for donor track
                        if bidirectional {
                            // remove old positions in donor site
                            for element_idx in start_donor_site..=end_donor_site {
                                let element_id = donor_genome.seq[element_idx].element_id;
                                let homology_group = thread_homology_updates
                                    .entry((element_id, donor_genome.genome_id))
                                    .or_insert_with(|| population.homology_map[element_id][donor_genome.genome_id].clone());
                                homology_group.retain(|pos| *pos != element_idx); // remove old position
                            }

                            // update donor genome
                            donor_genome
                                .seq
                                .splice(start_donor_site..=end_donor_site, recipient_track);

                            // add new positions
                            for element_idx in start_donor_site..(start_donor_site + recipient_track_len) {
                                let element_id = donor_genome.seq[element_idx].element_id;
                                let homology_group = thread_homology_updates
                                    .entry((element_id, donor_genome.genome_id))
                                    .or_insert_with(|| population.homology_map[element_id][donor_genome.genome_id].clone());
                                homology_group.push(element_idx); // add new position
                            }
                            // update contig_ids
                            donor_genome.update_contig_starts();
                        }
                    }
                }
            }
        (genomes, thread_homology_updates, thread_total_donor_length, thread_total_recipient_length, thread_successful_recombinations)
    }).collect::<Vec<_>>();

    // combine results from each thread
    let mut successful_recombinations = 0;
    let mut total_donor_length = 0;
    let mut total_recipient_length = 0;

    // staging vec for indexed insertion by genome_id
    let mut new_pop: Vec<Option<Genome>> = (0..pop_size).map(|_| None).collect();

    for (genomes, thread_homology_updates, thread_donor_length, thread_recipient_length, thread_successful_recombinations) in results {
        total_donor_length += thread_donor_length;
        total_recipient_length += thread_recipient_length;
        successful_recombinations += thread_successful_recombinations;

        // update population with new genomes
        for (genome_id, genome) in genomes.into_iter() {
            new_pop[genome_id] = Some(genome);
        }

        // only write back homology groups that changed in this component
        for ((element_id, genome_id), positions) in thread_homology_updates.into_iter() {
            population.homology_map[element_id][genome_id] = positions;
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
    use std::sync::Arc;

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
            contig_lengths: vec![0],
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
            seq: vec![
                NucElement {
                    contig_id: 0,
                    element_id: 0,
                    feature_id: 0,
                    feature_pos: 0,
                    feature_type: Arc::from("exon"),
                    multiplier: 1.0,
                    seq: Arc::new(vec![]),
                    mutation_map: Arc::new(MutationMap::new(0, 0, &vec![], &sel_dist, &mut rng)),
                    strand: true,
                    inverted: false,
                    original_length: 0,
                    frameshift: false,
                    tracked: false,
                    selection_coeff: 0.0,
                },
                NucElement {
                    contig_id: 0,
                    element_id: 1,
                    feature_id: 1,
                    feature_pos: 1,
                    feature_type: Arc::from("exon"),
                    multiplier: 1.0,
                    seq: Arc::new(vec![]),
                    mutation_map: Arc::new(MutationMap::new(0, 0, &vec![], &sel_dist, &mut rng)),
                    strand: false,
                    inverted: false,
                    original_length: 0,
                    frameshift: false,
                    tracked: false,
                    selection_coeff: 0.0,
                },
                NucElement {
                    contig_id: 0,
                    element_id: 2,
                    feature_id: 2,
                    feature_pos: 2,
                    feature_type: Arc::from("exon"),
                    multiplier: 1.0,
                    seq: Arc::new(vec![]),
                    mutation_map: Arc::new(MutationMap::new(0, 0, &vec![], &sel_dist, &mut rng)),
                    strand: false,
                    inverted: false,
                    original_length: 0,
                    frameshift: false,
                    tracked: false,
                    selection_coeff: 0.0,
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
                            feature_pos: usize,
                            strand: bool,
                            sel_dist: &MutationDistribution,
                            rng: &mut StdRng| NucElement {
            contig_id,
            element_id,
            feature_id,
            feature_pos,
            feature_type: Arc::from("exon"),
            multiplier: 1.0,
            seq: Arc::new(vec![]),
            mutation_map: Arc::new(MutationMap::new(0, 0, &vec![], sel_dist, rng)),
            strand,
            original_length: 0,
            frameshift: false,
            tracked: false,
            selection_coeff: 0.0,
            inverted: false,
        };

        Genome {
            identifier: "test_multi_contig".to_string(),
            genome_id: 0,
            parent: "root".to_string(),
            contig_starts: vec![0, 2, 4],
            contig_lengths: vec![0, 0, 0],
            seq: vec![
                make_element(0, 0, 0, 0, true, &sel_dist, &mut rng),
                make_element(0, 1, 1, 1, false, &sel_dist, &mut rng),
                make_element(1, 2, 2, 2, true, &sel_dist, &mut rng),
                make_element(1, 3, 3, 3, false, &sel_dist, &mut rng),
                make_element(2, 4, 4, 4, true, &sel_dist, &mut rng),
                make_element(2, 5, 5, 5, false, &sel_dist, &mut rng),
            ],
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
                feature_pos: idx,
                feature_type: Arc::from("exon"),
                multiplier: 1.0,
                seq: Arc::new(marker_seq.clone()),
                mutation_map: Arc::new(MutationMap::new(0, 0, &marker_seq, &sel_dist, &mut rng)),
                strand: strand_seed,
                original_length: marker_seq.len(),
                frameshift: false,
                tracked: false,
                selection_coeff: 0.0,
                inverted: false,
            });
        }

        let mut genome = Genome {
            identifier: format!("recomb_{}", genome_id),
            genome_id,
            parent: "root".to_string(),
            contig_starts: vec![0],
            contig_lengths: vec![0],
            seq,
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
        };
        genome.update_contig_starts();
        genome
    }

    fn print_genome(population: &Population, genome_idx: usize) -> String {
        let genome: &Genome = &population.pop[genome_idx];
        genome
            .seq
            .iter()
            .flat_map(|element| {
                if element.strand {
                    element.seq.as_ref().clone()
                } else {
                    reverse_complement(element.seq.as_slice())
                }
            })
            .map(|base| base.to_string())
            .collect::<Vec<String>>()
            .join("")
    }

    fn make_recombination_test_population(forced_events: usize, n_elements: usize) -> Population {
        // genome 0 starts with marker base 1 (A), genome 1 starts with marker base 2 (C)
        let g0 = make_recombination_test_genome(0, n_elements, true, 1);
        let g1 = make_recombination_test_genome(1, n_elements, true, 2);

        let recombination_count =
            MutationDistribution::new_uniform(forced_events as f64, forced_events as f64 + 0.1)
                .unwrap();
        let recombination_len = MutationDistribution::new_uniform(0.0, 0.1).unwrap();

        let mut homology_map: Vec<Vec<HomologyPositions>> = Vec::new();
        for idx in 0..n_elements {
            // Map each element_id to its actual position in each genome so
            // recombination start sites can vary across the genome.
            homology_map.push(vec![smallvec::smallvec![idx], smallvec::smallvec![idx]]);
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
            augment_tracking: false,
            genome_size_penalty_per_bp: 0.01,
            optimal_genome_size: 1000,
            compress_output: false,
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

    fn make_homology_test_element(seq: Vec<u8>, strand: bool) -> NucElement {
        let mut rng = StdRng::seed_from_u64(999);
        let sel_dist = MutationDistribution::new_uniform(0.0, 1.0).unwrap();

        NucElement {
            contig_id: 0,
            element_id: 0,
            feature_id: 0,
            feature_pos: 0,
            feature_type: Arc::from("exon"),
            multiplier: 1.0,
            seq: Arc::new(seq.clone()),
            mutation_map: Arc::new(MutationMap::new(0, 0, &seq, &sel_dist, &mut rng)),
            strand,
            original_length: seq.len(),
            frameshift: false,
            tracked: false,
            selection_coeff: 0.0,
            inverted: false,
        }
    }

    #[test]
    fn reverse_complement_produces_expected_sequence() {
        let seq = vec![1, 2, 4, 8, 16];
        let rc = reverse_complement(&seq);
        assert_eq!(
            rc,
            vec![16, 1, 2, 4, 8],
            "reverse complement should reverse sequence and swap one-hot nucleotide codes"
        );
    }

    #[test]
    fn homology_uses_reverse_complement_for_opposite_strands() {
        let forward = vec![1, 2, 4, 8, 1];
        let reversed_complement = reverse_complement(&forward);

        let query = make_homology_test_element(forward, true);

        // Same strand: do not reverse complement, so this should not be a perfect match.
        let subject_same_strand = make_homology_test_element(reversed_complement.clone(), true);
        let homology_without_rc = calculate_homology(&query, &subject_same_strand, 0.99);
        assert!(
            homology_without_rc < 1.0,
            "same-strand comparison should not reverse complement and should be imperfect for RC-only sequence"
        );

        // Opposite strand: reverse complement should be applied, yielding a perfect match.
        let subject_opposite_strand = make_homology_test_element(reversed_complement, false);
        let homology_with_rc = calculate_homology(&query, &subject_opposite_strand, 0.99);
        assert!(
            homology_with_rc >= 0.99,
            "opposite-strand comparison should reverse complement and recover a near-perfect match"
        );
    }

    #[test]
    fn inversion_flips_strand() {
        let mut genome = make_test_genome();
        let before_strands: Vec<bool> = genome.seq.iter().map(|e| e.strand).collect();

        let mut default_structural_dists = default_structural_dists();

        // update inversion rate to 1 for all elements, so that all elements are inverted
        default_structural_dists[0][2] = MutationDistribution::new_uniform(1.0, 1.1).unwrap();
        let pos = MutationDistribution::new_uniform(0.0, 1.0).unwrap();

        let mut homology_map: Vec<Vec<HomologyPositions>> = Vec::new();
        for _ in genome.seq.iter() {
            homology_map.push(vec![smallvec::smallvec![0]]);
        }

        mutate_intra_genome(&mut genome, &default_structural_dists, &pos, false);

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

        let mut homology_map: Vec<Vec<HomologyPositions>> = Vec::new();
        for _ in genome.seq.iter() {
            homology_map.push(vec![smallvec::smallvec![0]]);
        }

        let pos = MutationDistribution::new_uniform(0.0, 1.0).unwrap();
        mutate_intra_genome(&mut genome, &default_structural_dists, &pos, false);

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

        let mut homology_map: Vec<Vec<HomologyPositions>> = Vec::new();
        for _ in genome.seq.iter() {
            homology_map.push(vec![smallvec::smallvec![0]]);
        }

        // Use a non-zero offset so duplicates land somewhere other than position 0.
        let pos = MutationDistribution::new_uniform(1.0, 2.0).unwrap();
        mutate_intra_genome(&mut genome, &default_structural_dists, &pos, false);

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

        mutate_intra_genome(&mut genome, &default_structural_dists, &pos, false);

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

        println!("Genome 0 pre-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 pre-recomb: {}", print_genome(&population, 1));

        let total_before: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        mutate_inter_genome(&mut population, false);
        let total_after: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        let mixed_after = count_mixed_marker_genomes(&population);

        println!("Genome 0 post-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 post-recomb: {}", print_genome(&population, 1));

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

        println!("Genome 0 pre-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 pre-recomb: {}", print_genome(&population, 1));

        let total_before: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        mutate_inter_genome(&mut population, false);
        let total_after: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        let mixed_after = count_mixed_marker_genomes(&population);

        println!("Genome 0 post-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 post-recomb: {}", print_genome(&population, 1));

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
            mixed_after == 1,
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

        println!("Genome 0 pre-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 pre-recomb: {}", print_genome(&population, 1));

        let total_before: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        mutate_inter_genome(&mut population, false);
        let total_after: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        let mixed_after = count_mixed_marker_genomes(&population);

        println!("Genome 0 post-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 post-recomb: {}", print_genome(&population, 1));

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
            "after forced recombinations, at one genome should contain marker sequence from the other genome"
        );
    }

    #[test]
    fn inter_genome_recombination_single_event_is_bidirectional() {
        let mut population = make_recombination_test_population(1, 8);

        let mixed_before = count_mixed_marker_genomes(&population);
        assert_eq!(
            mixed_before, 0,
            "before recombination, genomes should not be mixed"
        );

        assert!(
            !genome_has_marker(&population.pop[0], 2),
            "before bidirectional recombination, genome 0 should not have marker sequence from genome 1"
        );
        assert!(
            !genome_has_marker(&population.pop[1], 1),
            "before ith bidirectional recombination, genome 1 should not have marker sequence from genome 0"
        );

        println!("Genome 0 pre-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 pre-recomb: {}", print_genome(&population, 1));

        let total_before: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        let (successful_recombinations, _, _) = mutate_inter_genome(&mut population, true);
        let total_after: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        let mixed_after = count_mixed_marker_genomes(&population);

        println!("Genome 0 post-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 post-recomb: {}", print_genome(&population, 1));

        assert_eq!(
            population.pop.len(),
            2,
            "population size should be unchanged"
        );
        assert_eq!(
            total_after, total_before,
            "single forced bidirectional recombination should preserve total genome length"
        );
        assert!(
            successful_recombinations >= 1,
            "at least one recombination should succeed in this deterministic setup"
        );
        assert!(
            genome_has_marker(&population.pop[0], 2),
            "with bidirectional recombination, genome 0 should gain marker sequence from genome 1"
        );
        assert!(
            genome_has_marker(&population.pop[1], 1),
            "with bidirectional recombination, genome 1 should gain marker sequence from genome 0"
        );
        assert!(
            mixed_after == 2,
            "after forced recombinations, both genomes should contain marker sequence from the other genome"
        );
    }

    #[test]
    fn inter_genome_recombination_multiple_events_are_bidirectional() {
        let forced_events = 3;
        let mut population = make_recombination_test_population(forced_events, 8);

        let mixed_before = count_mixed_marker_genomes(&population);
        assert_eq!(
            mixed_before, 0,
            "before recombination, genomes should not be mixed"
        );

        assert!(
            !genome_has_marker(&population.pop[0], 2),
            "before bidirectional recombination, genome 0 should not have marker sequence from genome 1"
        );
        assert!(
            !genome_has_marker(&population.pop[1], 1),
            "before ith bidirectional recombination, genome 1 should not have marker sequence from genome 0"
        );

        println!("Genome 0 pre-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 pre-recomb: {}", print_genome(&population, 1));

        let total_before: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        let (successful_recombinations, _, _) = mutate_inter_genome(&mut population, true);
        let total_after: usize = population.pop.iter().map(|g| g.seq.len()).sum();
        let mixed_after = count_mixed_marker_genomes(&population);

        println!("Genome 0 post-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 post-recomb: {}", print_genome(&population, 1));

        assert_eq!(
            population.pop.len(),
            2,
            "population size should be unchanged"
        );
        assert_eq!(
            total_after, total_before,
            "forced multiple bidirectional recombinations should preserve total genome length"
        );
        assert!(
            successful_recombinations >= 1,
            "at least one recombination should succeed in this deterministic setup"
        );
        assert!(
            genome_has_marker(&population.pop[0], 2),
            "with bidirectional recombination, genome 0 should gain marker sequence from genome 1"
        );
        assert!(
            genome_has_marker(&population.pop[1], 1),
            "with bidirectional recombination, genome 1 should gain marker sequence from genome 0"
        );
        assert!(
            mixed_after == 2,
            "after forced recombinations, both genomes should contain marker sequence from the other genome"
        );
    }

    #[test]
    fn inter_genome_multiple_events_recombine_more_sites_than_single_event() {
        let n_replicates = 30;
        let n_elements = 8;

        let mut single_total_foreign_sites = 0usize;
        let mut multiple_total_foreign_sites = 0usize;

        for _ in 0..n_replicates {
            let mut single_pop = make_recombination_test_population(1, n_elements);
            mutate_inter_genome(&mut single_pop, false);

            let single_foreign_sites: usize = single_pop
                .pop
                .iter()
                .enumerate()
                .map(|(genome_idx, genome)| {
                    let foreign_marker = if genome_idx == 0 { 2 } else { 1 };
                    genome
                        .seq
                        .iter()
                        .filter(|element| element.seq.first().copied() == Some(foreign_marker))
                        .count()
                })
                .sum();
            single_total_foreign_sites += single_foreign_sites;

            let mut multiple_pop = make_recombination_test_population(5, n_elements);
            mutate_inter_genome(&mut multiple_pop, false);

            let multiple_foreign_sites: usize = multiple_pop
                .pop
                .iter()
                .enumerate()
                .map(|(genome_idx, genome)| {
                    let foreign_marker = if genome_idx == 0 { 2 } else { 1 };
                    genome
                        .seq
                        .iter()
                        .filter(|element| element.seq.first().copied() == Some(foreign_marker))
                        .count()
                })
                .sum();
            multiple_total_foreign_sites += multiple_foreign_sites;
        }

        assert!(
            multiple_total_foreign_sites > single_total_foreign_sites,
            "across replicates, multiple forced recombinations should yield more recombined sites than single-event runs (single total: {}, multiple total: {})",
            single_total_foreign_sites,
            multiple_total_foreign_sites
        );
    }

    #[test]
    fn recombination_length_relative_to_contig_controls_track_coverage() {
        let n_elements = 8;

        // Small minimum recombination length relative to contig length should
        // produce a partial-track swap (not the whole contig).
        let mut short_len_population = make_recombination_test_population(1, n_elements);
        short_len_population.recombination_dists[1] =
            MutationDistribution::new_uniform(1.0, 1.1).unwrap();

        let (short_successful_recombinations, _, _) =
            mutate_inter_genome(&mut short_len_population, false);
        assert!(
            short_successful_recombinations >= 1,
            "expected at least one successful recombination with short minimum length"
        );

        let short_foreign_sites: usize = short_len_population
            .pop
            .iter()
            .enumerate()
            .map(|(genome_idx, genome)| {
                let foreign_marker = if genome_idx == 0 { 2 } else { 1 };
                genome
                    .seq
                    .iter()
                    .filter(|element| element.seq.first().copied() == Some(foreign_marker))
                    .count()
            })
            .sum();

        assert!(
            short_foreign_sites > 0 && short_foreign_sites < n_elements,
            "short minimum recombination length should recombine part, not all, of a contig (foreign sites: {}, contig elements: {})",
            short_foreign_sites,
            n_elements
        );

        // Large minimum recombination length relative to contig length should
        // force a whole-contig swap.
        let mut long_len_population = make_recombination_test_population(1, n_elements);
        let contig_length_bp: usize = long_len_population.pop[0]
            .seq
            .iter()
            .map(|element| element.seq.len())
            .sum();
        let huge_min_track = (contig_length_bp * 10) as f64;
        long_len_population.recombination_dists[1] =
            MutationDistribution::new_uniform(huge_min_track, huge_min_track + 0.1).unwrap();

        let (long_successful_recombinations, _, _) =
            mutate_inter_genome(&mut long_len_population, false);
        assert!(
            long_successful_recombinations >= 1,
            "expected at least one successful recombination with large minimum length"
        );

        let long_foreign_sites: usize = long_len_population
            .pop
            .iter()
            .enumerate()
            .map(|(genome_idx, genome)| {
                let foreign_marker = if genome_idx == 0 { 2 } else { 1 };
                genome
                    .seq
                    .iter()
                    .filter(|element| element.seq.first().copied() == Some(foreign_marker))
                    .count()
            })
            .sum();

        assert!(long_foreign_sites > short_foreign_sites, 
            "Should be more recombined sites in longer track recombination.");

        assert_eq!(
            long_foreign_sites,
            n_elements,
            "large minimum recombination length should recombine the entire contig (foreign sites: {}, contig elements: {})",
            long_foreign_sites,
            n_elements
        );
    }

    #[test]
    fn whole_genome_inversion_prevents_recombination() {
        let forced_events = 5;
        let n_elements = 8;

        let mut g0 = make_recombination_test_genome(0, n_elements, true, 1);
        let mut g1 = make_recombination_test_genome(1, n_elements, false, 1);

        let recombination_count =
            MutationDistribution::new_uniform(forced_events as f64, forced_events as f64 + 0.1)
                .unwrap();
        let recombination_len = MutationDistribution::new_uniform(0.0, 0.1).unwrap();

        let mut homology_map: Vec<Vec<HomologyPositions>> = Vec::new();
        for idx in 0..n_elements {
            homology_map.push(vec![smallvec::smallvec![idx], smallvec::smallvec![idx]]);
        }

        let mut population = Population {
            id: 0,
            generation: 0,
            pop: vec![g0, g1],
            core_vec: vec![],
            selection_dists: vec![],
            mu_dists: vec![],
            indel_dists: vec![],
            structural_mu_dists: vec![vec![]],
            recombination_dists: vec![recombination_count, recombination_len],
            // Require exact sequence identity at candidate sites.
            recombination_threshold: 1.0,
            homology_map,
            feature_map: HashMap::new(),
            max_multiplier_dist: 10,
            n_generations: 10,
            verbose: true,
            augment_tracking: false,
            genome_size_penalty_per_bp: 0.01,
            optimal_genome_size: 1000,
            compress_output: false,
        };

        println!("Genome 0 pre-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 pre-recomb: {}", print_genome(&population, 1));

        let (successful_recombinations, _, _) = mutate_inter_genome(&mut population, false);

        println!("Genome 0 post-recomb: {}", print_genome(&population, 0));
        println!("Genome 1 post-recomb: {}", print_genome(&population, 1));

        assert_eq!(
            successful_recombinations, 0,
            "whole-genome inversion with strict homology should prevent recombination"
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
            contig_lengths: vec![0],
            seq: vec![NucElement {
                contig_id: 0,
                element_id: 0,
                feature_id: 0,
                feature_pos: 0,
                feature_type: Arc::from(element_type),
                multiplier: 1.0,
                seq: Arc::new(vec![1, 2, 4, 8]),
                mutation_map: Arc::new(MutationMap::new(0, 0, &vec![1, 2, 4, 8], &sel_dist, &mut rng)),
                strand: true,
                inverted: false,
                original_length: 4,
                frameshift: false,
                tracked: false,
                selection_coeff: 0.0,
            },
            NucElement {
                contig_id: 0,
                element_id: 0,
                feature_id: 0,
                feature_pos: 0,
                feature_type: Arc::from("exon"),
                multiplier: 1.0,
                seq: Arc::new(vec![1, 2, 4, 8]),
                mutation_map: Arc::new(MutationMap::new(0, 0, &vec![1, 2, 4, 8], &sel_dist, &mut rng)),
                strand: true,
                inverted: false,
                original_length: 4,
                frameshift: false,
                tracked: false,
                selection_coeff: 0.0,
            },
            NucElement {
                contig_id: 0,
                element_id: 0,
                feature_id: 0,
                feature_pos: 0,
                feature_type: Arc::from("intergenic"),
                multiplier: 1.0,
                seq: Arc::new(vec![1, 2, 4, 8]),
                mutation_map: Arc::new(MutationMap::new(0, 0, &vec![1, 2, 4, 8], &sel_dist, &mut rng)),
                strand: true,
                original_length: 4,
                frameshift: false,
                tracked: false,
                inverted: false,
                selection_coeff: 0.0,
            }],
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
        
        mutate_intra_genome(&mut genome, &default_structural_dists, &pos, false);
        
        // TE-COPY should result in multiple copies (original + duplicates)
        assert!(
            genome.seq.len() > before_len,
            "TE-COPY should allow multiple duplications; before: {}, after: {}",
            before_len,
            genome.seq.len()
        );
        
        // All copies should be TE-COPY
        let te_copy_count = genome.seq.iter().filter(|e| e.feature_type.as_ref() == "TE-COPY").count();
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
        
        // try randomly permuting sequence to maximum attempts to pass all tests
        let n_attempts = 100;

        for _ in 0..n_attempts {
            mutate_intra_genome(&mut genome, &default_structural_dists, &pos, false);

            let te_position = genome.seq.iter().position(|e| e.feature_type.as_ref() == "TE-CUT");
            if te_position.expect("TE-CUT should still be present after cut-and-paste") != 0 {
                break;
            }
        }

        // After cut-and-paste, genome should have same or fewer elements
        // (original deleted, one copy inserted)
        assert!(
            genome.seq.len() == before_len,
            "TE-CUT cut-and-paste should result in at most one additional element; before: {}, after: {}",
            before_len,
            genome.seq.len()
        );

        println!("Genome after TE-CUT mutation: {:?}", genome.seq.iter().map(|e| e.feature_type.to_string()).collect::<Vec<String>>());

        // ensure TE has moved and original position is deleted
        let te_position = genome.seq.iter().position(|e| e.feature_type.as_ref() == "TE-CUT");
        assert_ne!(
            te_position.expect("TE-CUT should still be present after cut-and-paste"),
            0,
            "TE-CUT should have moved from original position"
        );

        let te_cut_count = genome.seq.iter().filter(|e| e.feature_type.as_ref() == "TE-CUT").count();
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
        
        mutate_intra_genome(&mut genome, &default_structural_dists, &pos, false);
        
        // Should have duplications
        assert!(
            genome.seq.len() > before_len,
            "intergenic should allow multiple duplications; before: {}, after: {}",
            before_len,
            genome.seq.len()
        );
        
        // All should be intergenic
        let intergenic_count = genome.seq.iter().filter(|e| e.feature_type.as_ref() == "intergenic").count();
        assert_eq!(
            intergenic_count > 1,
            true,
            "Intergenic should be duplicated"
        );
    }
}
