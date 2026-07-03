executable=../pansimnuc/target/release/PansimNuc
config_dir=../test_configs
outdir=/data/sam/analysis/PansimNuc/parameter_sweep
SNAKEMAKE_DIR="/home/sam/Software/PansimNucWF"
SNAKMAKE_CONFIG="/home/sam/Software/PansimNucWF/config/config.yaml"

mkdir -p $outdir

#baseline
# for config in "$config_dir"/baseline*.conf; do
#     run_name=$(basename "$config" .conf)
#     outdir_run=$outdir/$run_name
#     mkdir -p "$outdir_run"
#     echo "Running $run_name..."
#     $executable --config "$config" > "$outdir_run/$run_name.log"
#     Rscript "$DFE_script" "$outdir_run/selection_samples.csv" "$outdir_run/DFE_plot"
#     Rscript "$allele_freq_script" "$outdir_run/tracking.csv" "$outdir_run/allele_freq_plot" 3
#     Rscript "$TE_script" "$outdir_run" "$outdir_run/TE_copy_numbers"
#     Rscript "$sfs_script" "$outdir_run" "$outdir_run/sfs_plot"
#     Rscript "$haplotype_script" "$outdir_run/tracking.csv" "$outdir_run/haplotype_analysis" 5
#     Rscript "$plotting_script" "$outdir_run/root_out.gff" "$outdir_run" --out "$outdir_run/sv_plot.pdf" --width 16 --height 16 --types exon,intron,intergenic,TE-COPY,TE-CUT --link-types exon,TE-COPY,TE-CUT --gap 50000
# done

# #specific runs
for config in "$config_dir"/baseline_uniform_selection_one_split_gen_25*.conf; do
    run_name=$(basename "$config" .conf)
    outdir_run=$outdir/$run_name
    mkdir -p "$outdir_run"
    echo "Running $run_name..."
    #$executable --config "$config" > "$outdir_run/$run_name.log"

    echo "Running $run_name with PansimNucWF..."
    cd /home/sam/Software/PansimNucWF && pixi run snakemake -s ${SNAKEMAKE_DIR}/Snakefile --directory ${outdir_run} --configfile ${SNAKMAKE_CONFIG} --cores 4 --use-conda --conda-prefix ${SNAKEMAKE_DIR}/.snakemake/conda --config reference=${outdir_run}/root.fasta input_dir=${outdir_run} output_dir=${outdir_run} simulated="True" --unlock
    cd /home/sam/Software/PansimNucWF && pixi run snakemake -s ${SNAKEMAKE_DIR}/Snakefile --directory ${outdir_run} --configfile ${SNAKMAKE_CONFIG} --cores 4 --use-conda --conda-prefix ${SNAKEMAKE_DIR}/.snakemake/conda --config reference=${outdir_run}/root.fasta input_dir=${outdir_run} output_dir=${outdir_run} simulated="True" > "${outdir_run}/pansimnucWF.log" 2>&1

    #Rscript "$DFE_script" "$outdir_run/selection_samples.csv" "$outdir_run/DFE_plot"
    #Rscript "$allele_freq_script" "$outdir_run/tracking.csv" "$outdir_run/allele_freq_plot" 3
    #Rscript "$TE_script" "$outdir_run" "$outdir_run/TE_copy_numbers"
    #Rscript "$sfs_script" "$outdir_run" "$outdir_run/sfs_plot"
    #Rscript "$haplotype_script" "$outdir_run/tracking.csv" "$outdir_run/haplotype_analysis" 5
    #Rscript "$plotting_script" "$outdir_run/root.gff" "$outdir_run" --out "$outdir_run/sv_plot.pdf" --width 16 --height 16 --types exon,intron,intergenic,TE-COPY,TE-CUT --link-types exon,TE-COPY,TE-CUT --gap 50000
done

#all
# for config in "$config_dir"/*.conf; do
#     run_name=$(basename "$config" .conf)
#     outdir_run=$outdir/$run_name
#     mkdir -p "$outdir_run"
#     echo "Running $run_name..."
#     $executable --config "$config" > "$outdir_run/$run_name.log"
#     Rscript "$DFE_script" "$outdir_run/selection_samples.csv" "$outdir_run/DFE_plot"
#     Rscript "$allele_freq_script" "$outdir_run/tracking.csv" "$outdir_run/allele_freq_plot" 3
#     Rscript "$TE_script" "$outdir_run" "$outdir_run/TE_copy_numbers"
#     Rscript "$sfs_script" "$outdir_run" "$outdir_run/sfs_plot"
#     Rscript "$haplotype_script" "$outdir_run/tracking.csv" "$outdir_run/haplotype_analysis" 5
#     Rscript "$plotting_script" "$outdir_run/root.gff" "$outdir_run" --out "$outdir_run/sv_plot.pdf" --width 16 --height 16 --types exon,intron,intergenic,TE-COPY,TE-CUT --link-types exon,TE-COPY,TE-CUT --gap 50000
# done