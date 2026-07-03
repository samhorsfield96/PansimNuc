#!/usr/bin/env python3
"""Grid (combinatorial) parameter sweep generator for PansimNuc.

Generates a CSV summarising all parameter combinations, plus an individual
.conf file for every row so that each Slurm array job can start immediately.

Usage
-----
    python grid_parameter_sweep.py \\
        --config example.conf \\
        --output grid_samples.csv \\
        --configs-dir grid_configs/ \\
        --outdir-base grid_sweep/ \\
        --param exons.mutation_rate:1e-10,1e-9,1e-8 \\
        --param population.n_generations:100,500,1000 \\
        --param exons.selection_distribution:uniform,normal

Parameter spec format
---------------------
    SECTION.KEY:VAL1,VAL2,VAL3,...   (set key only in the named section)
    KEY:VAL1,VAL2,VAL3,...           (set key in every section that contains it)

SECTION and KEY must match the section headers and option names in the
PansimNuc .conf file exactly (case-sensitive).  Omitting the SECTION prefix
is useful for parameters like ``mutation_rate`` that appear in multiple
sections and should be swept together.

All combinations of all parameter values are generated (full Cartesian product).
"""

from __future__ import annotations

import argparse
import csv
import itertools
import re
import sys
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_param_spec(spec: str) -> dict[str, Any]:
    """Parse one --param specification into a structured dict."""
    if ":" not in spec:
        raise ValueError(
            f"Invalid parameter spec: {spec!r}. "
            "Expected SECTION.KEY:VAL1,VAL2,..."
        )

    key_path, values_str = spec.split(":", 1)

    if "." in key_path:
        section, key = key_path.split(".", 1)
    else:
        section = None  # broadcast: apply to every section that has this key
        key = key_path
    values = [v.strip() for v in values_str.split(",") if v.strip()]

    if not values:
        raise ValueError(f"No values specified for {key_path!r}")

    return {
        "section": section,
        "key": key,
        "values": values,
        "name": key_path,
    }


# ---------------------------------------------------------------------------
# Config file generation
# ---------------------------------------------------------------------------

def read_base_lines(config_path: Path) -> list[str]:
    with open(config_path) as fh:
        return fh.readlines()


def existing_params(base_lines: list[str]) -> dict[str, set[str]]:
    """Return {section: {key, ...}} for all key=value lines in the base config."""
    result: dict[str, set[str]] = {}
    current_section: str | None = None
    for line in base_lines:
        m = re.match(r"^\[([^\]]+)\]", line.strip())
        if m:
            current_section = m.group(1).strip()
            result.setdefault(current_section, set())
        elif current_section and "=" in line and not line.strip().startswith("#"):
            key = line.strip().split("=", 1)[0].strip()
            result[current_section].add(key)
    return result


def write_config(
    base_lines: list[str],
    param_overrides: dict[tuple[str, str], str],
    run_outdir: str,
    out_path: Path,
) -> None:
    """Write a modified config file, substituting param values and outdir."""
    current_section: str | None = None
    result: list[str] = []

    for line in base_lines:
        stripped = line.strip()

        # Section header
        m = re.match(r"^\[([^\]]+)\]", stripped)
        if m:
            current_section = m.group(1).strip()
            result.append(line)
            continue

        # Blank lines and comments pass through unchanged
        if not stripped or stripped.startswith("#"):
            result.append(line)
            continue

        # Key=value line
        if "=" in stripped:
            key = stripped.split("=", 1)[0].strip()

            if current_section == "output" and key == "outdir":
                result.append(f"outdir={run_outdir}\n")
                continue

            if current_section and (current_section, key) in param_overrides:
                result.append(f"{key}={param_overrides[(current_section, key)]}\n")
                continue

            if current_section and (None, key) in param_overrides:
                result.append(f"{key}={param_overrides[(None, key)]}\n")
                continue

        result.append(line)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.writelines(result)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate a full-grid combinatorial parameter sweep for PansimNuc.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--config", "-c",
        required=True,
        metavar="FILE",
        help="Base PansimNuc .conf file used as the template for all runs.",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        metavar="CSV",
        help="Path for the output CSV file (run_id + parameter values).",
    )
    parser.add_argument(
        "--configs-dir",
        required=True,
        metavar="DIR",
        help="Directory where per-run .conf files (run_0.conf, run_1.conf, …) are written.",
    )
    parser.add_argument(
        "--outdir-base",
        required=True,
        metavar="DIR",
        help=(
            "Base directory for PansimNuc output. Each run writes to "
            "<outdir-base>/run_<id>/. Must match OUTDIR_BASE in the Slurm script."
        ),
    )
    parser.add_argument(
        "--param", "-p",
        action="append",
        dest="params",
        default=[],
        metavar="SECTION.KEY:VAL1,VAL2,...",
        help=(
            "Parameter to vary with explicit comma-separated values. "
            "Can be repeated. Numeric and string (e.g. distribution names) values are both supported."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.params:
        parser.error("Specify at least one --param SECTION.KEY:VAL1,VAL2,...")

    # Parse and validate parameter specs
    params = []
    for raw in args.params:
        try:
            params.append(parse_param_spec(raw))
        except ValueError as exc:
            parser.error(str(exc))

    base_config = Path(args.config)
    if not base_config.exists():
        parser.error(f"Config file not found: {base_config}")

    base_lines = read_base_lines(base_config)
    known = existing_params(base_lines)

    for param in params:
        sec, key = param["section"], param["key"]
        if sec is None:
            # broadcast: check that the key exists in at least one section
            if not any(key in keys for keys in known.values()):
                print(
                    f"Warning: key '{key}' not found in any section of base config.",
                    file=sys.stderr,
                )
        elif sec not in known:
            print(
                f"Warning: section [{sec}] not found in base config.",
                file=sys.stderr,
            )
        elif key not in known[sec]:
            print(
                f"Warning: key '{key}' not found in [{sec}] of base config.",
                file=sys.stderr,
            )

    # Cartesian product of all parameter value lists
    value_lists = [p["values"] for p in params]
    combinations = list(itertools.product(*value_lists))
    n_runs = len(combinations)

    print(
        f"Parameters: {[p['name'] for p in params]}\n"
        f"Value counts: {[len(p['values']) for p in params]}\n"
        f"Total combinations: {n_runs}"
    )

    # Write CSV
    csv_path = Path(args.output)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["run_id"] + [p["name"] for p in params])
        for run_id, combo in enumerate(combinations):
            writer.writerow([run_id] + list(combo))

    # Write per-run config files
    configs_dir = Path(args.configs_dir)
    outdir_base = Path(args.outdir_base)

    for run_id, combo in enumerate(combinations):
        overrides = {
            (p["section"], p["key"]): v
            for p, v in zip(params, combo)
        }
        run_outdir = str(outdir_base / f"run_{run_id}")
        write_config(
            base_lines,
            overrides,
            run_outdir,
            configs_dir / f"run_{run_id}.conf",
        )

    print(
        f"Generated {n_runs} runs.\n"
        f"  CSV     : {csv_path}\n"
        f"  Configs : {configs_dir}/run_0.conf … run_{n_runs - 1}.conf"
    )


if __name__ == "__main__":
    main()
