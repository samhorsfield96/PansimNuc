use noodles_gff::{self as gff};
use noodles_gff::feature::record::Strand;
use noodles_fasta as fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};

#[derive(Clone)]
struct TeInterval {
    start: usize,
    end: usize,
    feature_type: String,
    strand: bool,
}

pub struct FeaturePos {
    pub contig_id: usize,
    pub feature_id: usize,
    pub feature_type: String,
    pub start: usize,
    pub end: usize,
    pub strand: bool, // true for +, false for -
    pub seq: Vec<u8>
}

fn encode_dna(seq: &str) -> Vec<u8> {
    seq.bytes()
        .map(|b| match b {
            b'A' => 1,
            b'C' => 2,
            b'G' => 4,
            b'T' => 8,
            _ => 16, // N or any other non-ACGT character
        })
        .collect()
}

fn classify_te_feature_type(raw_type: &str) -> Option<String> {
    let upper = raw_type.to_ascii_uppercase();

    if upper.contains("UNCLASSIFIED") {
        return None;
    }

    if upper.contains("DNA") || upper.contains("MITE") {
        Some("TE-CUT".to_string())
    } else {
        Some("TE-COPY".to_string())
    }
}

fn get_contig_order_from_gff(gff_path: &str) -> io::Result<Vec<String>> {
    let file_gff = File::open(gff_path)?;
    let mut gff_reader = gff::io::Reader::new(BufReader::new(file_gff));
    let mut contigs: Vec<String> = Vec::new();

    for result in gff_reader.record_bufs() {
        let record: noodles_gff::feature::RecordBuf = result?;
        let seqname = record.reference_sequence_name().to_string();
        if contigs.last().map_or(true, |last| last != &seqname) {
            contigs.push(seqname);
        }
    }

    Ok(contigs)
}

fn parse_earlgrey_intervals(
    earlgrey_gff_path: &str,
    contig_map: &HashMap<String, usize>,
) -> io::Result<HashMap<usize, Vec<TeInterval>>> {
    let file_gff = File::open(earlgrey_gff_path)?;
    let mut gff_reader = gff::io::Reader::new(BufReader::new(file_gff));
    let mut intervals_by_contig: HashMap<usize, Vec<TeInterval>> = HashMap::new();

    for result in gff_reader.record_bufs() {
        let record: noodles_gff::feature::RecordBuf = result?;
        let seqname = record.reference_sequence_name().to_string();
        let Some(&contig_id) = contig_map.get(&seqname) else {
            continue;
        };

        let Some(feature_type) = classify_te_feature_type(&record.ty().to_string()) else {
            continue;
        };

        let start = usize::from(record.start()).saturating_sub(1);
        let end = usize::from(record.end());
        if start >= end {
            continue;
        }

        intervals_by_contig
            .entry(contig_id)
            .or_default()
            .push(TeInterval {
                start,
                end,
                feature_type,
                strand: record.strand() == Strand::Forward,
            });
    }

    for intervals in intervals_by_contig.values_mut() {
        intervals.sort_by_key(|interval| interval.start);
    }

    Ok(intervals_by_contig)
}

fn push_feature_segment(
    out: &mut Vec<FeaturePos>,
    contig_id: usize,
    feature_id: usize,
    feature_type: &str,
    start: usize,
    end: usize,
    strand: bool,
    seq: &str,
) {
    if start >= end || end > seq.len() {
        return;
    }

    out.push(FeaturePos {
        contig_id,
        feature_id,
        feature_type: feature_type.to_string(),
        start,
        end,
        strand,
        seq: encode_dna(&seq[start..end]),
    });
}

fn overlay_te_intervals(
    features: &mut Vec<FeaturePos>,
    intervals: &[TeInterval],
    contig_id: usize,
    contig_seq: &str,
) {
    for interval in intervals {
        let mut updated: Vec<FeaturePos> = Vec::new();

        for feature in &*features {
            let overlap_start = feature.start.max(interval.start);
            let overlap_end = feature.end.min(interval.end);

            if overlap_start >= overlap_end {
                updated.push(FeaturePos {
                    contig_id: feature.contig_id,
                    feature_id: feature.feature_id,
                    feature_type: feature.feature_type.clone(),
                    start: feature.start,
                    end: feature.end,
                    strand: feature.strand,
                    seq: feature.seq.clone(),
                });
                continue;
            }

            push_feature_segment(
                &mut updated,
                contig_id,
                feature.feature_id,
                &feature.feature_type,
                feature.start,
                overlap_start,
                feature.strand,
                contig_seq,
            );

            push_feature_segment(
                &mut updated,
                contig_id,
                0,
                &interval.feature_type,
                overlap_start,
                overlap_end,
                interval.strand,
                contig_seq,
            );

            push_feature_segment(
                &mut updated,
                contig_id,
                feature.feature_id,
                &feature.feature_type,
                overlap_end,
                feature.end,
                feature.strand,
                contig_seq,
            );
        }

        *features = updated;
    }
}

pub fn extract_feature_positions(file_gff: File) -> io::Result<Vec<Vec<FeaturePos>>> {
    let mut gff_reader = gff::io::Reader::new(BufReader::new(file_gff));

    // keep track of current feature ID, dictacted by gene and its upstream region
    let mut current_feature_id: usize = 1;
    
    // hold features
    let mut features: Vec<Vec<FeaturePos>> = Vec::new();
    let mut contig_id: i32 = -1;
    let mut prev_seqname: String = String::new();

    for result in gff_reader.record_bufs() {
        let record: noodles_gff::feature::RecordBuf = result?;
        let feature_start: usize = usize::from(record.start()).saturating_sub(1);
        let feature_end: usize = usize::from(record.end());
        let seqname: String = record.reference_sequence_name().to_string();
        let feature_type: String = record.ty().to_string();

        if prev_seqname != seqname {
            contig_id += 1;
            prev_seqname = seqname.clone();
            features.push(Vec::new());
        }

        // get last record if present
        if let Some(last_feature) = features[contig_id as usize].last_mut() {

            // determine if there is overlap
            let last_feature_end = last_feature.end.clone();
            let last_feature_type = last_feature.feature_type.clone();
            let last_feature_id = last_feature.feature_id.clone();

            // if at new gene that is non-overlapping, previous region must be intergenic
            if feature_type == "gene" && last_feature_type == "exon" {
                current_feature_id += 1;

                // add intergenic region between last exon and current gene, if non-overlapping
                if feature_start >= last_feature_end {
                    features[contig_id as usize]
                        .push(FeaturePos {
                            contig_id: contig_id as usize,
                            feature_id: 0,
                            feature_type: "intergenic".to_string(),
                            start: last_feature_end,
                            end: feature_start,
                            strand: record.strand() == Strand::Forward,
                            seq: vec![0]
                        });
                }
            } 
            // non-overlapping exon after intergenic region
            else if feature_type == "exon" && feature_start >= last_feature_end && current_feature_id != last_feature_id {
                // if no intergenic region, need to add due to gene overlap
                if last_feature_type != "intergenic" {
                    features[contig_id as usize]
                        .push(FeaturePos {
                            contig_id: contig_id as usize,
                            feature_id: 0,
                            feature_type: "intergenic".to_string(),
                            start: last_feature_end,
                            end: feature_start,
                            strand: record.strand() == Strand::Forward,
                            seq: vec![0]
                        });
                } else {
                    // if intergenic region, need to update end coordinate to current exon start
                    last_feature.end = feature_start;
                }

                features[contig_id as usize]
                    .push(FeaturePos {
                        contig_id: contig_id as usize,
                        feature_id: current_feature_id,
                        feature_type: feature_type.clone(),
                        start: feature_start,
                        end: feature_end,
                        strand: record.strand() == Strand::Forward,
                        seq: vec![0]
                });
            } else if feature_type == "exon" && feature_start >= last_feature_end && last_feature_type == "exon" {

                // add intron feature between last exon and current exon
                features[contig_id as usize]
                    .push(FeaturePos {
                        contig_id: contig_id as usize,
                        feature_id: current_feature_id,
                        feature_type: "intron".to_string(),
                        start: last_feature_end,
                        end: feature_start,
                        strand: record.strand() == Strand::Forward,
                        seq: vec![0]
                    });

                // only add exons as features
                features[contig_id as usize]
                    .push(FeaturePos {
                        contig_id: contig_id as usize,
                        feature_id: current_feature_id,
                        feature_type: feature_type.clone(),
                        start: feature_start,
                        end: feature_end,
                        strand: record.strand() == Strand::Forward,
                        seq: vec![0]
                });    
            } 
        } else { 
            // if no last feature, add intergenic region from start of contig to first feature
            if feature_start > 0 {
                features[contig_id as usize]
                    .push(FeaturePos {
                        contig_id: contig_id as usize,
                        feature_id: 0,
                        feature_type: "intergenic".to_string(),
                        start: 0,
                        end: feature_start,
                        strand: true,
                        seq: vec![0]
                    });
            } else  {
                // unless if first feature starts at 0 add feature
                if feature_type == "gene" {
                    current_feature_id += 1;
                }
                else if feature_type == "exon" {
                    features[contig_id as usize]
                        .push(FeaturePos {
                            contig_id: contig_id as usize,
                            feature_id: current_feature_id,
                            feature_type: feature_type,
                            start: feature_start,
                            end: feature_end,
                            strand: record.strand() == Strand::Forward,
                            seq: vec![0]
                        });
                }
            }
        }
    }

    Ok(features)
}

pub fn read_gff_lines(
    gff_path: &str,
    fasta_path: &str,
    earlgrey_gff_path: Option<&str>,
) -> io::Result<Vec<Vec<FeaturePos>>> {
    let file_gff = File::open(gff_path)?;
    let file_fasta = File::open(fasta_path)?;

    let mut features = extract_feature_positions(file_gff)?;

    // Load FASTA records into memory keyed by contig name.
    let mut fasta_reader = fasta::io::Reader::new(BufReader::new(file_fasta));

    let mut genome: Vec<String> = Vec::new();

    for result in fasta_reader.records() {
        let record = result?;
        genome.push(String::from_utf8_lossy(record.sequence().as_ref()).into_owned());
    }

    for (contig_id, results) in features.iter_mut().enumerate() {
        if let Some(seq) = genome.get(contig_id) {
            let mut last_feature_end: usize = 0;
            
            for result in &mut **results {
            
                if result.start >= result.end || result.end > seq.len() {
                    continue;
                }

                let subseq = &seq[result.start..result.end];

                result.seq = encode_dna(subseq);

                last_feature_end = result.end;
            }
            
            // add final intergenic region, if contig empty adds full contig
            let len_seq: usize  = seq.len();
            let feature_start = last_feature_end;
            let feature_end = len_seq;
            let subseq = encode_dna(&seq[feature_start..feature_end]);

            results.push(FeaturePos {
                contig_id: contig_id,
                feature_id: 0,
                feature_type: "intergenic".to_string(),
                start: feature_start,
                end: feature_end,
                strand: true,
                seq: subseq
            });
        }
    }

    if let Some(earlgrey_path) = earlgrey_gff_path {
        let contig_order = get_contig_order_from_gff(gff_path)?;
        let contig_map: HashMap<String, usize> = contig_order
            .iter()
            .enumerate()
            .map(|(idx, name)| (name.clone(), idx))
            .collect();

        let intervals_by_contig = parse_earlgrey_intervals(earlgrey_path, &contig_map)?;

        for (contig_id, results) in features.iter_mut().enumerate() {
            let Some(seq) = genome.get(contig_id) else {
                continue;
            };

            let Some(intervals) = intervals_by_contig.get(&contig_id) else {
                continue;
            };

            overlay_te_intervals(results, intervals, contig_id, seq);
        }
    }

    Ok(features)
}

#[cfg(debug_assertions)]
pub fn write_root_genome_gff(features: &[Vec<FeaturePos>], output_path: &str) -> io::Result<()> {
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "##gff-version 3")?;

    let mut record_id: usize = 1;
    for (contig_idx, contig_features) in features.iter().enumerate() {
        for feature in contig_features {
            if feature.start >= feature.end {
                continue;
            }

            let seq_id = format!("contig_{}", contig_idx + 1);
            let start_1based = feature.start + 1;
            let end_1based = feature.end;
            let strand = if feature.strand { "+" } else { "-" };
            let attributes = format!(
                "ID=root_feature_{};feature_id={};feature_type={}",
                record_id, feature.feature_id, feature.feature_type
            );

            writeln!(
                writer,
                "{}\tPansimNuc\t{}\t{}\t{}\t.\t{}\t.\t{}",
                seq_id,
                feature.feature_type,
                start_1based,
                end_1based,
                strand,
                attributes
            )?;

            record_id += 1;
        }
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn create_temp_file(prefix: &str, suffix: &str, content: &str) -> String {
        let temp_dir = std::env::temp_dir();
        let temp_path = temp_dir.join(format!("{}_{}{}", prefix, std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos(), suffix));
        let mut file = File::create(&temp_path).expect("Failed to create temp file");
        file.write_all(content.as_bytes()).expect("Failed to write temp file");
        drop(file);
        temp_path.to_string_lossy().to_string()
    }

    #[test]
    fn test_read_multi_contig_gff() {
        let gff_content = "##gff-version 3
contig1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1
contig1\t.\texon\t100\t200\t.\t+\t.\tID=exon1
contig2\t.\tgene\t50\t150\t.\t+\t.\tID=gene2
contig2\t.\texon\t50\t150\t.\t+\t.\tID=exon2";

        let fasta_content = ">contig1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>contig2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let gff_file = create_temp_file("test", ".gff", gff_content);
        let fasta_file = create_temp_file("test", ".fasta", fasta_content);

        let result = read_gff_lines(&gff_file, &fasta_file, None);
        assert!(result.is_ok());

        let features = result.unwrap();
        assert_eq!(features.len(), 2);

        // Check coordinates for contig1 exon
        let contig1_features = &features[0];
        let exon = contig1_features.iter().find(|f| f.feature_type == "exon").unwrap();
        assert_eq!(exon.start, 99);  // GFF 100 becomes 0-indexed 99
        assert_eq!(exon.end, 200);

        // Check coordinates for contig2 exon
        let contig2_features = &features[1];
        let exon = contig2_features.iter().find(|f| f.feature_type == "exon").unwrap();
        assert_eq!(exon.start, 49);  // GFF 50 becomes 0-indexed 49
        assert_eq!(exon.end, 150);

        let _ = std::fs::remove_file(&gff_file);
        let _ = std::fs::remove_file(&fasta_file);
    }

    #[test]
    fn test_extract_multi_contig_features() {
        let gff_content = "##gff-version 3
contig1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1
contig1\t.\texon\t100\t200\t.\t+\t.\tID=exon1
contig2\t.\tgene\t50\t150\t.\t+\t.\tID=gene2
contig2\t.\texon\t50\t150\t.\t+\t.\tID=exon2";

        let gff_file = create_temp_file("test", ".gff", gff_content);
        let file = File::open(&gff_file).expect("Failed to open test file");

        let result = extract_feature_positions(file);
        assert!(result.is_ok());

        let features = result.unwrap();

        // Check coordinates for contig1 exon
        let contig1_features = &features[0];
        let exon = contig1_features.iter().find(|f| f.feature_type == "exon").unwrap();
        assert_eq!(exon.start, 99);  // GFF 100 becomes 0-indexed 99
        assert_eq!(exon.end, 200);

        // Check coordinates for contig2 exon
        let contig2_features = &features[1];
        let exon = contig2_features.iter().find(|f| f.feature_type == "exon").unwrap();
        assert_eq!(exon.start, 49);  // GFF 50 becomes 0-indexed 49
        assert_eq!(exon.end, 150);

        let _ = std::fs::remove_file(&gff_file);
    }

    #[test]
    fn test_earlgrey_overlay_replaces_overlaps_and_skips_unclassified() {
        let gff_content = "##gff-version 3
contig1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1
contig1\t.\texon\t100\t200\t.\t+\t.\tID=exon1";

        let mut fasta_content = String::from(">contig1\n");
        fasta_content.push_str(&"A".repeat(260));

        let earlgrey_content = "##gff-version 3
contig1\t.\tDNA\t120\t130\t.\t+\t.\tID=te_cut
contig1\t.\tLINE\t10\t20\t.\t-\t.\tID=te_copy
contig1\t.\tUnclassified\t140\t145\t.\t+\t.\tID=skip_me";

        let gff_file = create_temp_file("test", ".gff", gff_content);
        let fasta_file = create_temp_file("test", ".fasta", &fasta_content);
        let earlgrey_file = create_temp_file("test", ".earlgrey.gff", earlgrey_content);

        let result = read_gff_lines(&gff_file, &fasta_file, Some(&earlgrey_file));
        assert!(result.is_ok());

        let features = result.unwrap();
        assert_eq!(features.len(), 1);
        let contig_features = &features[0];

        // DNA -> TE-CUT and overlaps exon coordinates.
        let te_cut = contig_features
            .iter()
            .find(|f| f.feature_type == "TE-CUT" && f.start == 119 && f.end == 130)
            .expect("Expected TE-CUT segment at overlapped exon coordinates");
        assert!(te_cut.strand);

        // LINE -> TE-COPY and replaces intergenic segment.
        let te_copy = contig_features
            .iter()
            .find(|f| f.feature_type == "TE-COPY" && f.start == 9 && f.end == 20)
            .expect("Expected TE-COPY segment in intergenic coordinates");
        assert!(!te_copy.strand);

        // Exon should be split around TE-CUT overlap.
        assert!(contig_features
            .iter()
            .any(|f| f.feature_type == "exon" && f.start == 99 && f.end == 119));
        assert!(contig_features
            .iter()
            .any(|f| f.feature_type == "exon" && f.start == 130 && f.end == 200));

        // Unclassified entries must not be added.
        assert!(!contig_features
            .iter()
            .any(|f| f.feature_type.to_ascii_uppercase().contains("UNCLASSIFIED")));

        let _ = std::fs::remove_file(&gff_file);
        let _ = std::fs::remove_file(&fasta_file);
        let _ = std::fs::remove_file(&earlgrey_file);
    }
}
