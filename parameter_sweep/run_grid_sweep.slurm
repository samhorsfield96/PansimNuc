#!/bin/bash
# =============================================================================
# Slurm array job – PansimNuc Latin Hypercube parameter sweep
#
# Prerequisite: run lhs_parameter_sweep.py first to generate
#   1. A CSV of sampled parameter values  (e.g. lhs_samples.csv)
#   2. Per-run config files               (e.g. lhs_configs/run_0.conf …)
#
# Example:
#   python lhs_parameter_sweep.py \
#       --config baseline.conf \
#       --samples 100 \
#       --output lhs_samples.csv \
#       --configs-dir lhs_configs/ \
#       --outdir-base lhs_sweep/ \
#       --param exons.mutation_rate:1e-10:1e-8:log \
#       --param population.n_generations:100:1000 \
#       --seed 42
#
# Submit (adjust --array upper bound to N_SAMPLES-1):
#   sbatch --array=0-99 run_lhs_sweep.slurm
#
# Or set the array range here and submit without the flag:
#   sbatch run_lhs_sweep.slurm
# =============================================================================

#SBATCH --job-name=pansimnuc_grid
#SBATCH --output=/scratch/sam_simulations/logs/%A_%a.out
#SBATCH --error=/scratch/sam_simulations/logs/%A_%a.err
##SBATCH --array=0-242%10             # <-- Set to 0-(N_SAMPLES-1)
##SBATCH --array=0-2915%10             # <-- Set to 0-(N_SAMPLES-1)
##SBATCH --array=0-728%10             # <-- Set to 0-(N_SAMPLES-1)
#SBATCH --array=0-1151%10             # <-- Set to 0-(N_SAMPLES-1)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4        # Should match threads= in your .conf
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --partition=normal.12h   # Uncomment and set if required
##SBATCH --mail-user=samuel.horsfield@unine.ch       # Uncomment and set if required
##SBATCH --mail-type=END,FAIL

# =============================================================================
# User configuration – edit these paths before submitting
# =============================================================================
EXECUTABLE="/home/sam/Software/PansimNuc/pansimnuc/target/release/PansimNuc"

#change for each simulation
#RELATIVE_SOURCE_PATH="mu_vs_selection_vs_recombination_grid"
#RELATIVE_SOURCE_PATH="mu_vs_selection_vs_recombination_vs_TE_grid"
#RELATIVE_SOURCE_PATH="mu_vs_selection_vs_recombination_vs_migration_grid"
RELATIVE_SOURCE_PATH="AGR_grid"

# Directory containing run_0.conf, run_1.conf, … (--configs-dir from lhs script)
CONFIGS_DIR="/data/sam/analysis/PansimNuc/grid_parameter_sweep/${RELATIVE_SOURCE_PATH}"

# Base output directory for PansimNuc results (--outdir-base from lhs script)
# Each run writes to $OUTDIR_BASE/run_<id>/
OUTDIR_BASE="/scratch/sam_simulations/grid_parameter_sweep/${RELATIVE_SOURCE_PATH}"

mkdir -p $OUTDIR_BASE
# =============================================================================

set -euo pipefail

RUN_ID=${SLURM_ARRAY_TASK_ID}
RUN_CONFIG="${CONFIGS_DIR}/run_${RUN_ID}.conf"
RUN_OUTDIR="${OUTDIR_BASE}/run_${RUN_ID}"

# Sanity checks
if [[ ! -x "$EXECUTABLE" ]]; then
    echo "ERROR: Executable not found or not executable: $EXECUTABLE" >&2
    exit 1
fi

if [[ ! -f "$RUN_CONFIG" ]]; then
    echo "ERROR: Config file not found: $RUN_CONFIG" >&2
    echo "       Did you run lhs_parameter_sweep.py first?" >&2
    exit 1
fi

mkdir -p "$RUN_OUTDIR"

echo "========================================================"
echo "Job ID   : ${SLURM_JOB_ID:-local}_${RUN_ID}"
echo "Run ID   : ${RUN_ID}"
echo "Config   : ${RUN_CONFIG}"
echo "Output   : ${RUN_OUTDIR}"
echo "Started  : $(date -Is)"
echo "========================================================"

"$EXECUTABLE" --config "$RUN_CONFIG" \
    > "${RUN_OUTDIR}/pansimnuc.log" 2>&1

EXIT_CODE=$?

echo "Finished : $(date -Is)"
echo "Exit code: ${EXIT_CODE}"

exit ${EXIT_CODE}
