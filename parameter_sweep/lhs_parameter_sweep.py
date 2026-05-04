#!/usr/bin/env python3
"""Latin Hypercube Sampling parameter sweep generator for PansimNuc.

Generates a CSV summarising all parameter combinations, plus an individual
.conf file for every row so that each Slurm array job can start immediately.

Requirements: numpy, scipy (both available via conda-forge / pixi).

Usage
-----
    python lhs_parameter_sweep.py \\
        --config baseline.conf \\
        --samples 100 \\
        --output lhs_samples.csv \\
        --configs-dir lhs_configs/ \\
        --outdir-base lhs_sweep/ \\
        --param exons.mutation_rate:1e-10:1e-8:log \\
        --param population.n_generations:100:1000 \\
        --param population.recombination_rate:1e-10:1e-6:log \\
        --seed 42

Parameter spec format
---------------------
    SECTION.KEY:LOW:HIGH         – sample uniformly in [LOW, HIGH]
    SECTION.KEY:LOW:HIGH:log     – sample uniformly in log10-space

SECTION and KEY must match the section headers and option names in the
PansimNuc .conf file exactly (case-sensitive).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats.qmc import LatinHypercube


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_param_spec(spec: str) -> dict[str, Any]:
    """Parse one --param specification into a structured dict."""
    parts = spec.split(":")
    if len(parts) not in (3, 4):
        raise ValueError(
            f"Invalid parameter spec: {spec!r}. "
            "Expected SECTION.KEY:LOW:HIGH or SECTION.KEY:LOW:HIGH:log"
        )

    key_path = parts[0]
    if "." not in key_path:
        raise ValueError(
            f"Parameter {key_path!r} must be in SECTION.KEY format "
            "(e.g. exons.mutation_rate)"
        )

    section, key = key_path.split(".", 1)
    low = float(parts[1])
    high = float(parts[2])
    log_scale = len(parts) == 4 and parts[3].strip().lower() == "log"

    if low >= high:
        raise ValueError(f"LOW must be strictly less than HIGH for {key_path!r}")
    if log_scale and (low <= 0 and high >= 0 or low >= 0 and high <= 0):
        raise ValueError(
            f"Log-scale sampling requires both bounds to be strictly positive or strictly negative for {key_path!r}"
        )

    return {
        "section": section,
        "key": key,
        "low": low,
        "high": high,
        "log_scale": log_scale,
        "name": key_path,
    }


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------

def sample_lhs(
    params: list[dict],
    n_samples: int,
    seed: int | None = None,
) -> np.ndarray:
    """Return an (n_samples, n_params) array of LHS-sampled values."""
    sampler = LatinHypercube(d=len(params), seed=seed)
    unit = sampler.random(n=n_samples)  # values in [0, 1)

    scaled = np.empty_like(unit)
    for i, param in enumerate(params):
        lo, hi = param["low"], param["high"]
        if param["log_scale"]:
            neg_lo = lo < 0

            log_lo, log_hi = np.log10(abs(lo)), np.log10(abs(hi))
            scaled[:, i] = 10.0 ** (unit[:, i] * (log_hi - log_lo) + log_lo)
            
            # account for negative bounds by flipping the sign after scaling
            if neg_lo:
                scaled[:, i] *= -1
        else:
            scaled[:, i] = unit[:, i] * (hi - lo) + lo

    return scaled


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


def format_value(value: float, log_scale: bool) -> str:
    """Format a float for writing into a config file."""
    if log_scale or abs(value) < 1e-3 or abs(value) >= 1e6:
        return f"{value:.6e}"
    return f"{value:.6g}"


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

        result.append(line)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.writelines(result)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate a Latin Hypercube parameter sweep for PansimNuc.",
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
        "--samples", "-n",
        type=int,
        required=True,
        metavar="N",
        help="Number of LHS samples (= number of parameter combinations).",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        metavar="CSV",
        help="Path for the output CSV file (run_id + sampled parameter values).",
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
        metavar="SECTION.KEY:LOW:HIGH[:log]",
        help=(
            "Parameter to vary. Can be repeated. "
            "Append ':log' to sample in log10-space (recommended for rates)."
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Integer random seed for reproducibility.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.params:
        parser.error("Specify at least one --param SECTION.KEY:LOW:HIGH[:log].")

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
        if sec not in known:
            print(
                f"Warning: section [{sec}] not found in base config.",
                file=sys.stderr,
            )
        elif key not in known[sec]:
            print(
                f"Warning: key '{key}' not found in [{sec}] of base config.",
                file=sys.stderr,
            )

    # Draw LHS samples
    samples = sample_lhs(params, args.samples, seed=args.seed)

    # Write CSV
    csv_path = Path(args.output)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["run_id"] + [p["name"] for p in params])
        for run_id, row in enumerate(samples):
            writer.writerow([run_id] + [format_value(v, params[i]["log_scale"]) for i, v in enumerate(row)])

    # Write per-run config files
    configs_dir = Path(args.configs_dir)
    outdir_base = Path(args.outdir_base)

    for run_id, row in enumerate(samples):
        overrides = {
            (p["section"], p["key"]): format_value(v, p["log_scale"])
            for p, v in zip(params, row)
        }
        run_outdir = str(outdir_base / f"run_{run_id}")
        write_config(
            base_lines,
            overrides,
            run_outdir,
            configs_dir / f"run_{run_id}.conf",
        )

    print(
        f"Generated {args.samples} LHS samples.\n"
        f"  CSV     : {csv_path}\n"
        f"  Configs : {configs_dir}/run_0.conf … run_{args.samples - 1}.conf"
    )


if __name__ == "__main__":
    main()
