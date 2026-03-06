use noodles_gff::{self as gff};
use noodles_gff::feature::record::Strand;
use noodles_fasta as fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader};

pub struct FeaturePos {
    pub seqname: String,
    pub feature_id: usize,
    pub feature_type: String,
    pub start: usize,
    pub end: usize,
    pub strand: bool, // true for +, false for -
    pub seq: String
}

pub fn extract_feature_positions(file_gff: File) -> io::Result<HashMap<String, Vec<FeaturePos>>> {
    let mut gff_reader = gff::io::Reader::new(BufReader::new(file_gff));

    // keep track of current feature ID, dictacted by gene and its upstream region
    let mut current_feature_id: usize = 0;

    // keep track of previous feature end
    let mut last_feature_end: usize = 0;
    
    // hold features
    let mut features: HashMap<String, Vec<FeaturePos>> = HashMap::new();

    for result in gff_reader.record_bufs() {
        let record: noodles_gff::feature::RecordBuf = result?;

        let feature_start: usize = usize::from(record.start()).saturating_sub(1);
        let feature_end: usize = usize::from(record.end());
        let seqname: String = record.reference_sequence_name().to_string();
        let feature_type: String = record.ty().to_string();

        // if last feature was a gene, then region must be intergenic
        if feature_type == "Gene" {
            // increment feature ID
            current_feature_id += 1;

            // in case gene is start of contig
            if feature_start > last_feature_end {
                features.entry(seqname.clone())
                    .or_default()
                    .push(FeaturePos {
                        seqname: seqname.clone(),
                        feature_id: current_feature_id,
                        feature_type: "intergenic".to_string(),
                        start: last_feature_end,
                        end: feature_start,
                        strand: true,
                        seq: "".to_string()
                    });
            }

            // update last_feature_end
            last_feature_end = feature_end;
        } 
        // if next feature is exon, need to check than current end is not identical otherwise still at start of gene
        else if feature_type == "exon" {
            
            // if next feature is exon, need to check than current end is not identical otherwise still at start of gene
            if feature_start > last_feature_end {
                features.entry(seqname.clone())
                    .or_default()
                    .push(FeaturePos {
                        seqname: seqname.clone(),
                        feature_id: current_feature_id,
                        feature_type: "intron".to_string(),
                        start: last_feature_end,
                        end: feature_start,
                        strand: true,
                        seq: "".to_string()
                    });
            }

            // only add exons as features
            if feature_type == "exon" {
                features.entry(seqname.clone())
                    .or_default()
                    .push(FeaturePos {
                        seqname: seqname.clone(),
                        feature_id: current_feature_id,
                        feature_type: feature_type.clone(),
                        start: feature_start,
                        end: feature_end,
                        strand: record.strand() == Strand::Forward,
                        seq: "".to_string()
                    });
            }

            // update last_feature_end
            last_feature_end = feature_end;
        }
        
    }

    Ok(features)
}

pub fn read_gff_lines(gff_path: &str, fasta_path: &str) -> io::Result<HashMap<String, Vec<FeaturePos>>> {
    let file_gff = File::open(gff_path)?;
    let file_fasta = File::open(fasta_path)?;

    let mut features = extract_feature_positions(file_gff)?;

    // Load FASTA records into memory keyed by contig name.
    let mut fasta_reader = fasta::io::Reader::new(BufReader::new(file_fasta));

    let mut genome: HashMap<String, String> = HashMap::new();

    for result in fasta_reader.records() {
        let record = result?;
        genome.insert(
            String::from_utf8_lossy(record.name()).into_owned(),
            String::from_utf8_lossy(record.sequence().as_ref()).into_owned(),
        );
    }

    for (seqname, results) in &mut features {
        if let Some(seq) = genome.get(seqname) {
            let mut last_feature_end: usize = 0;
            let mut last_feature_id:usize = 0;
            
            for result in &mut **results {
            
                if result.start >= result.end || result.end > seq.len() {
                    continue;
                }

                let subseq = &seq[result.start..result.end];

                result.seq = subseq.to_string();

                last_feature_end = result.end;
                last_feature_id = result.feature_id;
            }
            
            // add final intergenic region, if contig empty adds full contig
            let len_seq: usize  = seq.len();
            let feature_start = last_feature_end;
            let feature_end = len_seq;
            let subseq = &seq[feature_start..feature_end];

            results.push(FeaturePos {
                seqname: seqname.clone(),
                feature_id: last_feature_id + 1,
                feature_type: "intergenic".to_string(),
                start: feature_start,
                end: feature_end,
                strand: true,
                seq: subseq.to_string()
            });
        }
    }

    Ok(features)
}
