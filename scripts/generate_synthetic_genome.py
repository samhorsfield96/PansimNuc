#!/usr/bin/env python3
"""Generate synthetic FASTA + paired coding/Earl Grey GFF files from a simple config.

Config format (CSV-like, one feature per line):
    feature,start,end[,strand][,attributes]
or:
    contig,feature,start,end[,strand][,attributes]

Examples:
    gene,1,200,+,ID=gene1
    exon,1,20,+
    intron,21,180,+
    TE-COPY,300,650,-,ID=my_te_1
    chr2,LINE,1000,1500,-,ID=te_chr2_1

Rules:
- Blank lines and lines starting with '#' are ignored.
- start/end are 1-based inclusive coordinates.
- If start > end, values are swapped automatically.
- strand defaults to '+' when omitted.
- attributes are optional and should be semicolon-separated key=value pairs.
"""

from __future__ import annotations

import argparse
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List


TE_PREFIX_GUESSES = (
    "te",
    "ltr",
    "line",
    "sine",
    "mite",
    "dna/",
    "rc/",
    "tir",
)


@dataclass
class FeatureRecord:
    contig: str
    feature: str
    start: int
    end: int
    strand: str
    attributes: Dict[str, str]
    line_no: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_attributes(raw: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in raw.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key.strip()] = value.strip()
        else:
            attrs[item] = ""
    return attrs


def parse_config(config_path: Path, default_contig: str) -> List[FeatureRecord]:
    records: List[FeatureRecord] = []
    with config_path.open("r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = [part.strip() for part in line.split(",")]
            if len(parts) < 3:
                raise ValueError(
                    f"Invalid config at line {line_no}: expected at least 3 columns "
                    f"(feature,start,end), got: {raw_line.rstrip()}"
                )

            contig = default_contig
            feature = ""
            coords_offset = 0

            has_leading_contig = False
            if len(parts) >= 4:
                try:
                    int(parts[2])
                    int(parts[3])
                    has_leading_contig = True
                except ValueError:
                    has_leading_contig = False

            if has_leading_contig:
                contig = parts[0]
                feature = parts[1]
                coords_offset = 2
            else:
                feature = parts[0]
                coords_offset = 1

            if not contig:
                raise ValueError(f"Invalid config at line {line_no}: contig name is empty")
            if not feature:
                raise ValueError(f"Invalid config at line {line_no}: feature name is empty")

            try:
                start = int(parts[coords_offset])
                end = int(parts[coords_offset + 1])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid config at line {line_no}: start/end must be integers"
                ) from exc

            if start < 1 or end < 1:
                raise ValueError(f"Invalid config at line {line_no}: coordinates must be >= 1")

            if start > end:
                start, end = end, start

            strand = "+"
            attributes: Dict[str, str] = {}

            strand_idx = coords_offset + 2
            attrs_idx = coords_offset + 3
            if len(parts) > strand_idx and parts[strand_idx] in {"+", "-"}:
                strand = parts[strand_idx]
                if len(parts) > attrs_idx:
                    attributes = parse_attributes(",".join(parts[attrs_idx:]))
            elif len(parts) > strand_idx:
                attributes = parse_attributes(",".join(parts[strand_idx:]))

            records.append(
                FeatureRecord(
                    contig=contig,
                    feature=feature,
                    start=start,
                    end=end,
                    strand=strand,
                    attributes=attributes,
                    line_no=line_no,
                )
            )

    if not records:
        raise ValueError("Config file has no feature records")

    return records


def is_te_feature(feature: str) -> bool:
    lowered = feature.strip().lower()
    return lowered.startswith(TE_PREFIX_GUESSES)


def random_dna(length: int, seed: int | None) -> str:
    rng = random.Random(seed)
    bases = "ACGT"
    return "".join(rng.choice(bases) for _ in range(length))


def wrap_sequence(sequence: str, width: int = 80) -> List[str]:
    return [sequence[i : i + width] for i in range(0, len(sequence), width)]


def build_sequences(contig_lengths: Dict[str, int], seed: int | None) -> Dict[str, str]:
    rng = random.Random(seed)
    bases = "ACGT"
    sequences: Dict[str, str] = {}
    for contig in sorted(contig_lengths):
        length = contig_lengths[contig]
        sequences[contig] = "".join(rng.choice(bases) for _ in range(length))
    return sequences


def write_fasta(path: Path, sequences: Dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for contig in sorted(sequences):
            sequence = sequences[contig]
            handle.write(f">{contig} synthetic_dna length={len(sequence)}\n")
            for line in wrap_sequence(sequence):
                handle.write(f"{line}\n")


def format_attrs(attrs: Dict[str, str]) -> str:
    if not attrs:
        return "."

    chunks: List[str] = []
    for key, value in attrs.items():
        if value == "":
            chunks.append(key)
        else:
            chunks.append(f"{key}={value}")
    return ";".join(chunks)


def write_coding_gff(path: Path, contig_lengths: Dict[str, int], records: List[FeatureRecord]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for contig in sorted(contig_lengths):
            handle.write(f"##sequence-region {contig} 1 {contig_lengths[contig]}\n")
        for idx, record in enumerate(records, start=1):
            attrs = {
                "ID": record.attributes.get("ID", f"{record.feature}_{idx}"),
                "Name": record.attributes.get("Name", record.feature),
                "source_line": str(record.line_no),
            }
            for key, value in record.attributes.items():
                if key not in attrs:
                    attrs[key] = value

            cols = [
                record.contig,
                "Synthetic",
                record.feature,
                str(record.start),
                str(record.end),
                ".",
                record.strand,
                ".",
                format_attrs(attrs),
            ]
            handle.write("\t".join(cols) + "\n")


def write_earlgrey_gff(
    path: Path,
    te_records: List[FeatureRecord],
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for idx, record in enumerate(te_records, start=1):
            attrs = {
                "Tstart": "1",
                "Tend": str(record.length),
                "ID": record.attributes.get("ID", f"{idx}"),
                "shortTE": record.attributes.get("shortTE", "F"),
            }
            for key, value in record.attributes.items():
                if key not in attrs:
                    attrs[key] = value

            cols = [
                record.contig,
                "RepeatMasker",
                record.feature,
                str(record.start),
                str(record.end),
                str(record.length),
                record.strand,
                ".",
                format_attrs(attrs),
            ]
            handle.write("\t".join(cols) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate synthetic FASTA + coding GFF + Earl Grey TE GFF from a "
            "simple component-order config."
        )
    )
    parser.add_argument("--config", required=True, type=Path, help="Path to the component config")
    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
        help="Output prefix, e.g. testing/synth_chr19 (creates .fasta, _coding.gff, _earlgrey.gff)",
    )
    parser.add_argument("--seed", type=int, default=42, help="Optional random seed for reproducible DNA")

    args = parser.parse_args()

    config_path = args.config
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")

    records = parse_config(config_path, default_contig="1")

    te_records: List[FeatureRecord] = []
    coding_records: List[FeatureRecord] = []
    for record in records:
        if is_te_feature(record.feature):
            te_records.append(record)
        else:
            coding_records.append(record)

    contig_lengths: Dict[str, int] = {}
    for record in records:
        current = contig_lengths.get(record.contig, 0)
        if record.end > current:
            contig_lengths[record.contig] = record.end

    sequences = build_sequences(contig_lengths, args.seed)

    output_prefix = args.output_prefix
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    fasta_path = output_prefix.with_suffix(".fasta")
    coding_gff_path = output_prefix.parent / f"{output_prefix.name}_coding.gff"
    earlgrey_gff_path = output_prefix.parent / f"{output_prefix.name}_earlgrey.gff"

    write_fasta(fasta_path, sequences)
    write_coding_gff(coding_gff_path, contig_lengths, coding_records)
    write_earlgrey_gff(earlgrey_gff_path, te_records)

    print(f"Wrote FASTA: {fasta_path}")
    print(f"Wrote coding GFF: {coding_gff_path} ({len(coding_records)} features)")
    print(f"Wrote Earl Grey GFF: {earlgrey_gff_path} ({len(te_records)} TE features)")


if __name__ == "__main__":
    main()
