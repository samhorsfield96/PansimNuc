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

pub fn extract_feature_positions(file_gff: File) -> io::Result<Vec<FeaturePos>> {
    let mut gff_reader = gff::io::Reader::new(BufReader::new(file_gff));

    // keep track of current feature ID, dictacted by gene and its upstream region
    let mut current_feature_id: usize = 0;

    // keep track of previous feature end
    let mut last_feature_end: usize = 0;
    
    // hold features
    let mut features: Vec<FeaturePos> = Vec::new();

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
                features.push(FeaturePos {
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
                features.push(FeaturePos {
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
                features.push(FeaturePos {
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

pub fn read_gff_lines(gff_path: &str, fasta_path: &str) -> io::Result<Vec<FeaturePos>> {
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

    for result in &mut features {

        if let Some(seq) = genome.get(&result.seqname) {
            if result.start >= result.end || result.end > seq.len() {
                continue;
            }

            let subseq = &seq[result.start..result.end];

            result.seq = subseq.to_string();
        }
    }

    if cfg!(debug_assertions) {
        println!("{}\n{}\n{}\n{}\n{}\n{}\n{}", 
            features[0].seqname,
            features[0].feature_id,
            features[0].feature_type,
            features[0].start,
            features[0].end,
            features[0].strand,
            features[0].seq);
        println!("{}\n{}\n{}\n{}\n{}\n{}\n{}", 
            features[1].seqname,
            features[1].feature_id,
            features[1].feature_type,
            features[1].start,
            features[1].end,
            features[1].strand,
            features[1].seq);
        println!("{}\n{}\n{}\n{}\n{}\n{}\n{}", 
            features[2].seqname,
            features[2].feature_id,
            features[2].feature_type,
            features[2].start,
            features[2].end,
            features[2].strand,
            features[2].seq);
        println!("{}\n{}\n{}\n{}\n{}\n{}\n{}", 
            features[3].seqname,
            features[3].feature_id,
            features[3].feature_type,
            features[3].start,
            features[3].end,
            features[3].strand,
            features[3].seq);
    }

    Ok(features)
}
