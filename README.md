# PansimNuc
A nucleotide-level pangenome simulator.

## Summary

PansimNuc takes a root genome and simulates mutation (SNPs and indels), selection, recombination, gene and transposable element (TE) mobility and gain and loss, and demography using a Wright-Fisher simulation framework.

All simulation is done at the nucelotide level, providing annotated sequences of all members of the population throughout the simulation.

## Installation

First clone this repository: 

```
git clone https://github.com/samhorsfield96/PansimNuc.git && cd PansimNuc
```

Then set up for environment...

### with mamba/conda

Install rust (>=1.89.0) using [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), and activate the environment:

```
micromamba create -n PansimNuc
micromamba activate PansimNuc
micromamba install conda-forge::rust
```

### with pixi

Install [pixi](https://pixi.prefix.dev/latest/installation/) and activate the environment:

```
curl -fsSL https://pixi.sh/install.sh | sh
pixi shell
```

### Building the executable

If the environment set up correctly, you should see this in your terminal:

```
(PansimNuc) PansimNuc %
```

With the environment set up, build PansimNuc:

```
cd pansimnuc && cargo build --release && cd ../..
```

The executable will be in `PansimNuc/pansimnuc/target/release`. You can add this to your path to allow running at anytime:

```
export PATH="absolute/path/to/PansimNuc/pansimnuc/target/release:$PATH"
```

To make this permanent, add the above line to your `~/.bashrc` or `~/.zshrc` file.

To run the exectuable, use the command:

```
PansimNuc/pansimnuc/target/release/PansimNuc --config path/to/config.conf
```

See below for setting up your `config.conf` file.

## Configuration

PansimNuc is configured via a `.conf` file passed with `--config`. The file uses INI-style sections (`[section]`) with `key=value` pairs. Lines beginning with `#` are comments. See `example.conf` and files in `testing` directory for examples with setup.

---

## Input / Output parameters

### `[input]`

| Parameter | Description |
|---|---|
| `gff_file` | Path to the input GFF annotation file |
| `fasta_file` | Path to the input FASTA genome file |
| `earlgrey_gff_file` | (Optional) Path to an EarlGrey TE annotation GFF file |

### `[output]`

| Parameter | Description |
|---|---|
| `outdir` | Directory where all output files are written |

---

## Population-level parameters

### `[population]`

| Parameter | Description |
|---|---|
| `n_individuals` | Number of haploid individuals in the population |
| `n_generations` | Number of generations to simulate |
| `recombination_rate` | Number of recombination events per base per generation |
| `recombination_size_mean` | Mean recombination tract length (bp); tract lengths are Poisson-distributed |
| `recombination_threshold` | Minimum sequence homology (0–1) required for two sequences to recombine |
| `max_multiplier_dist` | Maximum distance (bp) between elements for TE multiplier/recombination pairing |
| `population_splits` | Comma-separated list of population sizes after each split event |
| `generation_splits` | Comma-separated list of generations at which population split events occur (must match length of `population_splits`) |
| `migration_rate` | Per-genome probability of migration between sub-populations each generation |
| `genome_size_penalty_per_bp` | (Optional) Fitness penalty per bp of genome size; applied as a multiplicative factor `(1 - penalty)^genome_size` |

---

## Per-region mutation, indel, selection and structural variation parameters

The following parameters apply independently to each functional region section (`[exons]`, `[introns]`, `[intergenic]`, `[TE-CUT]`, `[TE-COPY]`, and `[tracking]`).

### Nucleotide mutation and indels

| Parameter | Description |
|---|---|
| `mutation_rate` | Mean number of nucleotide substitutions per site per genome per generation (Poisson-distributed) |
| `indel_rate` | Mean number of insertion/deletion events per site per genome per generation (Poisson-distributed) |

### Selection

Each site's selection coefficient has a multiplicative effect of `1 + X` on fitness, where `X >= -1`. The distribution of selection coefficients is set with `selection_distribution`.

| Parameter | Description |
|---|---|
| `selection_distribution` | Distribution to draw selection coefficients from. One of: `normal`, `uniform`, `exp`, `double_exp`, `poisson`, `gamma` |

Required parameters by distribution:

| Distribution | Required parameters |
|---|---|
| `normal` | `selection_mean`, `selection_std_dev` |
| `uniform` | `selection_low`, `selection_high` |
| `exp` | `selection_lambda` |
| `double_exp` | `selection_lambda1` (rate of negative selection), `selection_lambda2` (rate of positive selection), `selection_cutoff` (proportion of positively selected mutations, 0–1) |
| `poisson` | `selection_lambda` |
| `gamma` | `selection_lambda` (scale), `selection_shape` |

### Intra-genomic structural variation

Rates are given in events per element per genome per generation. By default, counts are Poisson-distributed; if a `*_variance` key is also provided, a negative binomial distribution is used instead (allowing overdispersion).

| Parameter | Description |
|---|---|
| `duplication_rate` | Mean rate of element duplications |
| `duplication_variance` | (Optional) Variance of duplication counts; enables negative binomial sampling |
| `deletion_rate` | Mean rate of element deletions |
| `deletion_variance` | (Optional) Variance of deletion counts; enables negative binomial sampling |
| `inversion_rate` | Mean rate of element inversions |
| `inversion_variance` | (Optional) Variance of inversion counts; enables negative binomial sampling |

### TE copy-number multiplier (TE regions and `[tracking]` only)

The multiplier controls the relative transposition activity and is drawn from a gamma distribution.

| Parameter | Description |
|---|---|
| `multiplier_rate` | Rate (mean) parameter of the gamma distribution for the TE copy-number multiplier |
| `multiplier_scale` | Scale parameter of the gamma distribution for the TE copy-number multiplier |

---

## Tracking parameters

The `[tracking]` section defines one or more genomic regions to monitor over time. It accepts all per-region parameters above, plus:

| Parameter | Description |
|---|---|
| `contig` | Comma-separated list of contig identifiers to track |
| `start` | Comma-separated list of region start coordinates (must match length of `contig`) |
| `end` | Comma-separated list of region end coordinates (must match length of `contig`) |
| `augmentation` | If `true`, augmented tracking is enabled (boolean) |

WARNING: specifying too large of a region in your root genome will result in very large output files.

---

## Miscellaneous parameters

### `[misc]`

| Parameter | Description |
|---|---|
| `threads` | Number of threads to use (parallelised with Rayon) |
| `seed` | Integer seed for the random number generator (required for reproducibility) |
| `verbose` | If `true`, print detailed progress information (boolean) |
| `print_DFE` | If `true`, write 1000 samples from each selection distribution to `<outdir>/selection_samples.csv` (boolean) |
| `print_all_generations` | If `true`, write GFF and FASTA output at every generation rather than only at the end (boolean) |

---

## Functional region sections

PansimNuc models the following distinct genomic feature types. Each has its own independent set of mutation, indel, selection and structural variation parameters as described above.

| Section | Description |
|---|---|
| `[exons]` | Protein-coding exonic regions |
| `[introns]` | Intronic (non-coding, spliced) regions |
| `[intergenic]` | Intergenic (non-coding, non-spliced) regions |
| `[TE-CUT]` | Cut-and-paste transposable elements. When duplicated, the element is immediately deleted from its original location (transposition only) |
| `[TE-COPY]` | Copy-and-paste transposable elements. Duplications retain the element at the original location |
| `[tracking]` | User-defined regions to track allele frequencies and mutations over time. Supports multiple comma-separated regions via `contig`, `start`, and `end` |

## Output

PansimNuc outputs:

| File | Description |
|---|---|
| `[pop_X_gen_Y_genome_Z.gff]` | GFF file containing function annotations of all elements for genome Z of population X, generation Y |
| `[pop_X_gen_Y_genome_Z.fasta]` | Fasta file containing sequence for genome Z of population X, generation Y |
| `[root.gff]` | GFF file containing function annotations of all elements in the root genome for the simulation |
| `[root.fasta]` | Fasta file containing sequence for root genome for the simulation |
| `[selection_samples.csv]` | Sample of selection coefficients for all elements used in the simulation |
| `[tracking.csv]` | Tracking information for elements within specifie regions in `tracking` section |


### Explanation of the GFF files

Each line in the GFF output follows standard GFF3 column layout, with PansimNuc-specific attributes in the final column:

```
contig_0	PansimNuc	TE-COPY	1	838	.	-	.	genome_id=4;element_id=0;...
```

**Standard GFF3 columns (tab-separated):**

| Column | Example | Description |
|---|---|---|
| seqname | `contig_0` | Contig/chromosome name |
| source | `PansimNuc` | Always `PansimNuc` |
| feature | `TE-COPY` | Feature type (see [Functional region sections](#functional-region-sections)) |
| start | `1` | 1-based start coordinate of the element on the contig |
| end | `838` | 1-based end coordinate of the element on the contig |
| score | `.` | Not used (always `.`) |
| strand | `-` | Strand: `+` for forward, `-` for reverse |
| frame | `.` | Not used (always `.`) |
| attributes | `genome_id=4;...` | Semicolon-separated key=value attributes (see below) |

**Attributes:**

| Attribute | Description |
|---|---|
| `genome_id` | Integer ID of the genome within the population |
| `element_id` | Integer ID of this element within the genome |
| `feature_type` | Feature type, matching column 3 |
| `feature_id` | Integer ID of this element in the root genome (shared across all genomes for the same ancestral element) |
| `contig_id` | Integer ID of the contig on which this element resides |
| `parent` | Hyphen-separated list of generation indices tracing the lineage of this element back to the root |
| `multiplier` | Gamma-distributed copy-number multiplier (TE elements and tracking regions only) |
| `sequence_length` | Current length of the element sequence (bp); may differ from `end - start + 1` if indels have occurred |
| `log_genome_selection_coefficient` | Log-sum of per-site selection coefficients across the entire genome |
| `log_element_selection_coefficient` | Log-sum of per-site selection coefficients for this element only |