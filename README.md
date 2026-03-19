# PansimNuc
A nucleotide-level pangenome simulator.

## Selection distributions

Each feature section in the config can define its selection distribution with `selection_distribution`.
Supported values are `normal`, `uniform`, `exp`, `double_exp`, and `poisson`.

Required parameters by distribution:

- `normal`: `selection_mean`, `selection_std_dev`
- `uniform`: `selection_low`, `selection_high`
- `exp`: `selection_lambda`
- `double_exp`: `selection_lambda1`, `selection_lambda2`, `selection_cutoff`
- `poisson`: `selection_lambda`
- `gamma`: `selection_lamba`, `selection_shape`

Legacy configs remain supported: `selection_coefficient=<lambda>` is treated as exponential shorthand.
