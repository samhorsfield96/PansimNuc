use noodles_gff::{self as gff};
use noodles_gff::feature::record::Strand;
use noodles_fasta as fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader};

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

pub fn extract_feature_positions(file_gff: File) -> io::Result<Vec<Vec<FeaturePos>>> {
    let mut gff_reader = gff::io::Reader::new(BufReader::new(file_gff));

    // keep track of current feature ID, dictacted by gene and its upstream region
    let mut current_feature_id: usize = 0;

    // keep track of previous feature end
    let mut last_feature_end: usize = 0;
    
    // hold features
    let mut features: Vec<Vec<FeaturePos>> = Vec::new();
    let mut contig_id: i32 = -1;
    let mut prev_seqname: String = String::new();

    for result in gff_reader.record_bufs() {
        let record: noodles_gff::feature::RecordBuf = result?;
        let feature_start: usize = usize::from(record.start()).saturating_sub(1);
        let feature_end: usize = usize::from(record.end());
        let seqname: String = record.reference_sequence_name().to_string();

        if prev_seqname != seqname {
            contig_id += 1;
            prev_seqname = seqname.clone();
            features.push(Vec::new());
            last_feature_end = 0;
        }

        let feature_type: String = record.ty().to_string();

        // if last feature was a gene, then region must be intergenic
        if feature_type == "Gene" {
            // increment feature ID
            current_feature_id += 1;

            // in case gene is start of contig
            if feature_start > last_feature_end {
                features[contig_id as usize]
                    .push(FeaturePos {
                        contig_id: contig_id as usize,
                        feature_id: current_feature_id,
                        feature_type: "intergenic".to_string(),
                        start: last_feature_end,
                        end: feature_start,
                        strand: true,
                        seq: vec![0]
                    });
            }

            // update last_feature_end
            last_feature_end = feature_end;
        } 
        // if next feature is exon, need to check than current end is not identical otherwise still at start of gene
        else if feature_type == "exon" {
            
            // if next feature is exon, need to check than current end is not identical otherwise still at start of gene
            if feature_start > last_feature_end {
                features[contig_id as usize]
                    .push(FeaturePos {
                        contig_id: contig_id as usize,
                        feature_id: current_feature_id,
                        feature_type: "intron".to_string(),
                        start: last_feature_end,
                        end: feature_start,
                        strand: true,
                        seq: vec![0]
                    });
            }

            // only add exons as features
            if feature_type == "exon" {
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

            // update last_feature_end
            last_feature_end = feature_end;
        }
        
    }

    Ok(features)
}

pub fn read_gff_lines(gff_path: &str, fasta_path: &str) -> io::Result<Vec<Vec<FeaturePos>>> {
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
            let mut last_feature_id:usize = 0;
            
            for result in &mut **results {
            
                if result.start >= result.end || result.end > seq.len() {
                    continue;
                }

                let subseq = &seq[result.start..result.end];

                result.seq = encode_dna(subseq);

                last_feature_end = result.end;
                last_feature_id = result.feature_id;
            }
            
            // add final intergenic region, if contig empty adds full contig
            let len_seq: usize  = seq.len();
            let feature_start = last_feature_end;
            let feature_end = len_seq;
            let subseq = encode_dna(&seq[feature_start..feature_end]);

            results.push(FeaturePos {
                contig_id: contig_id,
                feature_id: last_feature_id + 1,
                feature_type: "intergenic".to_string(),
                start: feature_start,
                end: feature_end,
                strand: true,
                seq: subseq
            });
        }
    }

    Ok(features)
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
contig1\t.\tGene\t100\t200\t.\t+\t.\tID=gene1
contig1\t.\texon\t100\t200\t.\t+\t.\tID=exon1
contig2\t.\tGene\t50\t150\t.\t+\t.\tID=gene2
contig2\t.\texon\t50\t150\t.\t+\t.\tID=exon2";

        let fasta_content = ">contig1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>contig2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

        let gff_file = create_temp_file("test", ".gff", gff_content);
        let fasta_file = create_temp_file("test", ".fasta", fasta_content);

        let result = read_gff_lines(&gff_file, &fasta_file);
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
contig1\t.\tGene\t100\t200\t.\t+\t.\tID=gene1
contig1\t.\texon\t100\t200\t.\t+\t.\tID=exon1
contig2\t.\tGene\t50\t150\t.\t+\t.\tID=gene2
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
}
