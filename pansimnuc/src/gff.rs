use noodles_gff as gff;
use noodles_fasta as fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader};

pub struct FeatureSeq {
    pub seqname: String,
    pub feature_type: String,
    pub sequence: Vec<u8>,
}

pub fn read_gff_lines(gff_path: &str, fasta_path: &str) -> io::Result<Vec<FeatureSeq>> {
    let file_gff = File::open(gff_path)?;
    let file_fasta = File::open(fasta_path)?;

    // Load FASTA records into memory keyed by contig name.
    let mut fasta_reader = fasta::io::Reader::new(BufReader::new(file_fasta));

    let mut genome: HashMap<String, Vec<u8>> = HashMap::new();

    for result in fasta_reader.records() {
        let record = result?;
        genome.insert(
            String::from_utf8_lossy(record.name()).into_owned(),
            record.sequence().as_ref().to_vec(),
        );
    }

    // Parse GFF features and extract subsequences from loaded FASTA.
    let mut gff_reader = gff::io::Reader::new(BufReader::new(file_gff));

    let mut features = Vec::new();

    for result in gff_reader.record_bufs() {
        let record = result?;

        let seqname = record.reference_sequence_name().to_string();
        let start = usize::from(record.start()).saturating_sub(1);
        let end = usize::from(record.end());

        if let Some(seq) = genome.get(&seqname) {
            if start >= end || end > seq.len() {
                continue;
            }

            let subseq = seq[start..end].to_vec();

            features.push(FeatureSeq {
                seqname,
                feature_type: record.ty().to_string(),
                sequence: subseq,
            });
        }
    }

    Ok(features)
}
