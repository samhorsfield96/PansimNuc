#!/usr/bin/env python3
"""Extract GFF or FASTA records matching a specific chromosome/sequence identifier."""

import sys
import argparse


def extract_gff(filepath, chrom):
    """Print GFF lines where column 1 matches chrom, preserving ## headers."""
    with open(filepath) as f:
        for line in f:
            if line.startswith("##") and not line.startswith("##sequence-region"):
                print(line, end="")
            elif line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 2 and parts[1] == chrom:
                    print(line, end="")
            elif not line.startswith("#"):
                fields = line.split("\t")
                if fields and fields[0] == chrom:
                    print(line, end="")


def extract_fasta(filepath, chrom):
    """Print FASTA records whose header identifier (first word after >) matches chrom."""
    inside = False
    with open(filepath) as f:
        for line in f:
            if line.startswith(">"):
                header_id = line[1:].split()[0]
                inside = (header_id == chrom)
                if inside:
                    print(line, end="")
            elif inside:
                print(line, end="")


def detect_format(filepath):
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("##gff") or line.startswith("#"):
                return "gff"
            if line.startswith(">"):
                return "fasta"
            # GFF line: tab-separated with 9 fields
            if "\t" in line and len(line.split("\t")) == 9:
                return "gff"
            break
    return "unknown"


def main():
    parser = argparse.ArgumentParser(
        description="Extract GFF or FASTA records for a specific chromosome/sequence ID."
    )
    parser.add_argument("file", help="Input GFF or FASTA file")
    parser.add_argument("chrom", help="Chromosome/sequence identifier to extract")
    parser.add_argument(
        "--format", choices=["gff", "fasta"], help="Force input format (auto-detected by default)"
    )
    args = parser.parse_args()

    fmt = args.format or detect_format(args.file)
    if fmt == "gff":
        extract_gff(args.file, args.chrom)
    elif fmt == "fasta":
        extract_fasta(args.file, args.chrom)
    else:
        print(f"Could not detect format of {args.file}. Use --format gff or --format fasta.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()