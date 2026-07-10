#!/usr/bin/env python3
"""Generate synthetic FASTA + paired coding/Earl Grey GFF files from a simple config.

Config format (CSV-like, one feature per line):
    feature,start,end[,strand][,attributes]

Examples:
    gene,1,200,+,ID=gene1
    exon,1,20,+
    intron,21,180,+
    TE-COPY,300,650,-,ID=my_te_1

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


def parse_config(config_path: Path) -> List[FeatureRecord]:
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

            feature = parts[0]
            if not feature:
                raise ValueError(f"Invalid config at line {line_no}: feature name is empty")

            try:
                start = int(parts[1])
                end = int(parts[2])
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

            if len(parts) >= 4 and parts[3] in {"+", "-"}:
                strand = parts[3]
                if len(parts) >= 5:
                    attributes = parse_attributes(",".join(parts[4:]))
            elif len(parts) >= 4:
                attributes = parse_attributes(",".join(parts[3:]))

            records.append(
                FeatureRecord(
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


def is_te_feature(feature: str, explicit_te_features: set[str]) -> bool:
    lowered = feature.strip().lower()
    if lowered in explicit_te_features:
        return True
    return lowered.startswith(TE_PREFIX_GUESSES)


def random_dna(length: int, seed: int | None) -> str:
    rng = random.Random(seed)
    bases = "ACGT"
    return "".join(rng.choice(bases) for _ in range(length))


def wrap_sequence(sequence: str, width: int = 80) -> List[str]:
    return [sequence[i : i + width] for i in range(0, len(sequence), width)]


def write_fasta(path: Path, contig: str, sequence: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
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


def write_coding_gff(path: Path, contig: str, contig_len: int, records: List[FeatureRecord]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        handle.write(f"##sequence-region {contig} 1 {contig_len}\n")
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
                contig,
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
    contig: str,
    te_records: List[FeatureRecord],
    family_prefix: str,
) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for idx, record in enumerate(te_records, start=1):
            attrs = {
                "Tstart": "1",
                "Tend": str(record.length),
                "ID": record.attributes.get("ID", f"{family_prefix}_{idx}"),
                "shortTE": record.attributes.get("shortTE", "F"),
            }
            for key, value in record.attributes.items():
                if key not in attrs:
                    attrs[key] = value

            cols = [
                contig,
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
    parser.add_argument("--contig", default="contig_0", help="Contig/chromosome name for outputs")
    parser.add_argument("--seed", type=int, default=None, help="Optional random seed for reproducible DNA")
    parser.add_argument(
        "--te-features",
        default="TE-CUT,TE-COPY",
        help=(
            "Comma-separated feature names to force-classify as TE in Earl Grey output. "
            "Case-insensitive."
        ),
    )
    parser.add_argument(
        "--family-prefix",
        default="SYNTE_FAMILY",
        help="Prefix used when auto-generating Earl Grey TE IDs",
    )

    args = parser.parse_args()

    config_path = args.config
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")

    records = parse_config(config_path)
    explicit_te_features = {
        item.strip().lower() for item in args.te_features.split(",") if item.strip()
    }

    te_records: List[FeatureRecord] = []
    coding_records: List[FeatureRecord] = []
    for record in records:
        if is_te_feature(record.feature, explicit_te_features):
            te_records.append(record)
        else:
            coding_records.append(record)

    contig_len = max(record.end for record in records)
    sequence = random_dna(contig_len, args.seed)

    output_prefix = args.output_prefix
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    fasta_path = output_prefix.with_suffix(".fasta")
    coding_gff_path = output_prefix.parent / f"{output_prefix.name}_coding.gff"
    earlgrey_gff_path = output_prefix.parent / f"{output_prefix.name}_earlgrey.gff"

    write_fasta(fasta_path, args.contig, sequence)
    write_coding_gff(coding_gff_path, args.contig, contig_len, coding_records)
    write_earlgrey_gff(earlgrey_gff_path, args.contig, te_records, args.family_prefix)

    print(f"Wrote FASTA: {fasta_path}")
    print(f"Wrote coding GFF: {coding_gff_path} ({len(coding_records)} features)")
    print(f"Wrote Earl Grey GFF: {earlgrey_gff_path} ({len(te_records)} TE features)")


if __name__ == "__main__":
    main()
