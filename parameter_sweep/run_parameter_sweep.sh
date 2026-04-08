executable=../pansimnuc/target/release/PansimNuc
config_dir=../test_configs
outdir=../parameter_sweep
plotting_script=../sv_plot.R
allele_freq_script=../plot_allele_frequencies.R

mkdir -p $outdir

#baseline
for config in "$config_dir"/baseline*.conf; do
    run_name=$(basename "$config" .conf)
    outdir_run=$outdir/$run_name
    mkdir -p "$outdir_run"
    echo "Running $run_name..."
    $executable --config "$config" > "$outdir_run/$run_name.log"
    #Rscript "$plotting_script" "$outdir_run/root_out.gff" "$outdir_run" --out "$outdir_run/sv_plot.pdf" --width 16 --height 16 --types exon,intron,intergenic,TE-COPY,TE-CUT --link-types exon,TE-COPY,TE-CUT --gap 50000
    Rscript "$allele_freq_script" "$outdir_run/tracking.csv" "$outdir_run/allele_freq_plot" 3
done

# exon selection parameter sweep
for config in "$config_dir"/exon_*_selection*.conf; do
    run_name=$(basename "$config" .conf)
    outdir_run=$outdir/$run_name
    mkdir -p "$outdir_run"
    echo "Running $run_name..."
    $executable --config "$config" > "$outdir_run/$run_name.log"
    #Rscript "$plotting_script" "$outdir_run/root_out.gff" "$outdir_run" --out "$outdir_run/sv_plot.pdf" --width 16 --height 16 --types exon,intron,intergenic,TE-COPY,TE-CUT --link-types exon,TE-COPY,TE-CUT --gap 50000
    Rscript "$allele_freq_script" "$outdir_run/tracking.csv" "$outdir_run/allele_freq_plot" 3
done

# TE activity parameter sweep
# for config in "$config_dir"/TE-*.conf; do
#     run_name=$(basename "$config" .conf)
#     outdir_run=$outdir/$run_name
#     mkdir -p "$outdir_run"
#     echo "Running $run_name..."
#     $executable --config "$config" > "$outdir_run/$run_name.log"
#     #Rscript "$plotting_script" "$outdir_run/root_out.gff" "$outdir_run" --out "$outdir_run/sv_plot.pdf" --width 16 --height 16 --types exon,intron,intergenic,TE-COPY,TE-CUT --link-types exon,TE-COPY,TE-CUT --gap 50000
# done

# mutation rate/Recombination parameter sweep
for config in "$config_dir"/exon_mu_*_recombination_*.conf; do
    run_name=$(basename "$config" .conf)
    outdir_run=$outdir/$run_name
    mkdir -p "$outdir_run"
    echo "Running $run_name..."
    $executable --config "$config" > "$outdir_run/$run_name.log"
    #Rscript "$plotting_script" "$outdir_run/root_out.gff" "$outdir_run" --out "$outdir_run/sv_plot.pdf" --width 16 --height 16 --types exon,intron,intergenic,TE-COPY,TE-CUT --link-types exon,TE-COPY,TE-CUT --gap 50000
    Rscript "$allele_freq_script" "$outdir_run/tracking.csv" "$outdir_run/allele_freq_plot" 3
done