# mu vs selection vs recombination
python grid_parameter_sweep.py \
    --config /home/sam/Software/PansimNuc/publication_configs/mu_vs_selection_vs_recombination_base.conf \
    --output /home/sam/Software/PansimNuc/parameter_sweep/mu_vs_selection_vs_recombination_grid.csv \
    --configs-dir /data/sam/analysis/PansimNuc/grid_parameter_sweep/mu_vs_selection_vs_recombination_grid \
    --outdir-base /scratch/sam_simulations/grid_parameter_sweep/mu_vs_selection_vs_recombination_grid \
    --param input.gff_file:/home/sam/Software/PansimNuc/testing/Zymo_chr19.gff \
    --param input.fasta_file:/home/sam/Software/PansimNuc/testing/Zymo_chr19.fasta \
    --param input.earlgrey_gff_file:/home/sam/Software/PansimNuc/testing/Zymo_EarlGrey_chr19.gff \
    --param population.n_individuals:50 \
    --param population.n_generations:500 \
    --param mutation_rate:1e-10,1e-8,1e-6 \
    --param selection_distribution:double_exp \
    --param selection_lambda1:1e14,1e7,1e1 \
    --param selection_lambda2:1e14,1e7,1e1 \
    --param selection_cutoff:0.1,0.5,0.9 \
    --param misc.print_all_generations:true \
    --param population.recombination_rate:1e-8,1e-6,1e-4

# TE dynamics
python grid_parameter_sweep.py \
    --config /home/sam/Software/PansimNuc/publication_configs/mu_vs_selection_vs_recombination_base.conf \
    --output /home/sam/Software/PansimNuc/parameter_sweep/mu_vs_selection_vs_recombination_vs_TE_grid.csv \
    --configs-dir /data/sam/analysis/PansimNuc/grid_parameter_sweep/mu_vs_selection_vs_recombination_vs_TE_grid \
    --outdir-base /scratch/sam_simulations/grid_parameter_sweep/mu_vs_selection_vs_recombination_vs_TE_grid \
    --param input.gff_file:/home/sam/Software/PansimNuc/testing/Zymo_chr19.gff \
    --param input.fasta_file:/home/sam/Software/PansimNuc/testing/Zymo_chr19.fasta \
    --param input.earlgrey_gff_file:/home/sam/Software/PansimNuc/testing/Zymo_EarlGrey_chr19.gff \
    --param population.n_individuals:50 \
    --param population.n_generations:500 \
    --param mutation_rate:1e-8 \
    --param selection_distribution:double_exp \
    --param selection_lambda1:1e14,1e1 \
    --param selection_lambda2:1e14,1e1 \
    --param selection_cutoff:0.1,0.5,0.9 \
    --param misc.print_all_generations:false \
    --param population.recombination_rate:1e-8,1e-6,1e-4 \
    --param TE-CUT.duplication_rate:1e-30,1e-5,1e-3 \
    --param TE-CUT.deletion:1e-30,1e-5,1e-3 \
    --param TE-CUT.duplication_rate:1e-30,1e-5,1e-3 \
    --param TE-COPY.deletion:1e-30,1e-5,1e-3

# multiple populations
python grid_parameter_sweep.py \
    --config /home/sam/Software/PansimNuc/publication_configs/mu_vs_selection_vs_recombination_base.conf \
    --output /home/sam/Software/PansimNuc/parameter_sweep/mu_vs_selection_vs_recombination_vs_migration_grid.csv \
    --configs-dir /data/sam/analysis/PansimNuc/grid_parameter_sweep/mu_vs_selection_vs_recombination_vs_migration_grid \
    --outdir-base /scratch/sam_simulations/grid_parameter_sweep/mu_vs_selection_vs_recombination_vs_migration_grid \
    --param input.gff_file:/home/sam/Software/PansimNuc/testing/Zymo_chr19.gff \
    --param input.fasta_file:/home/sam/Software/PansimNuc/testing/Zymo_chr19.fasta \
    --param input.earlgrey_gff_file:/home/sam/Software/PansimNuc/testing/Zymo_EarlGrey_chr19.gff \
    --param population.n_individuals:50 \
    --param population.n_generations:500 \
    --param mutation_rate:1e-10,1e-8,1e-6 \
    --param selection_distribution:double_exp \
    --param selection_lambda1:1e14,1e7,1e1 \
    --param selection_lambda2:1e14,1e7,1e1 \
    --param selection_cutoff:0.1,0.5,0.9 \
    --param misc.print_all_generations:true \
    --param population.recombination_rate:1e-8,1e-6,1e-4 \
    --param population.population_splits:2 \
    --param population.generation_splits:250 \
    --param population.migration_rate:0.0,0.01,0.1