# PansimNuc
A nucleotide-level pangenome simulator.

## Nucleotide mutation rates

Each genome feature has it's own average nucleotide mutation rate, where the number of mutations is sampled from a Poisson distribution.

This is described by `mutation_rate`, and is measured in events per site per genome per generation.

## Nucleotide selection distributions

Each genome feature in the config can define its selection distribution with `selection_distribution`.
Supported values are `normal`, `uniform`, `exp`, `double_exp`, and `poisson`.

Each selection coefficient has a multiplicative effect of 1 + X on fitness, where X >= -1. 

Required parameters by distribution:

- `normal`: `selection_mean`, `selection_std_dev`
- `uniform`: `selection_low`, `selection_high`
- `exp`: `selection_lambda`
- `double_exp`: `selection_lambda1` (strenght of negative selection), `selection_lambda2` (strength of positive selection), `selection_cutoff` (proportion of positively selected genes)
- `poisson`: `selection_lambda`
- `gamma`: `selection_lamba`, `selection_shape`

## Intra-genomic variation

Each genome feature has it's own probability for duplications, deletions and inversions. Rates for each event are given in events per element per genome per generation.

Note that `TE-CUT` are only transposed, meaning that if they are duplicated, they are immediately deleted at their old location.

## Inter-genomic variation

Recombination is governed by the rate of recombinations, mean recombination length and sequence homology threshold. The rate of recombination in the population is governed by the number of SNPs, with `recombination_rate` giving the number of bases recombined per SNP. The `recombination_size_mean` dictates how long the recombination tracts are on average. The `recombination_threshold` dictates the lower limit of sequence homology for two sequences to recombine.
