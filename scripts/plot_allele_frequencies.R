library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)

# Usage:
#   Rscript plot_allele_frequencies.R [tracking.csv] [output.pdf]
# Defaults to tracking.csv in the current directory and allele_freq_plot

args          <- commandArgs(trailingOnly = TRUE)
# Filter out R's own flags (e.g. --no-save, --no-restore) that leak through
args          <- args[!grepl("^--", args)]
tracking_file <- if (length(args) >= 1) args[1] else "tracking.csv"
outpref   <- if (length(args) >= 2) args[2] else "allele_freq_plot"
top_n_val <- if (length(args) >= 3) as.numeric(args[3]) else 3

#tracking_file <- "/Users/samhorsfield/Software/PansimNuc/parameter_sweep/baseline/tracking.csv"
#outpref <- "/Users/samhorsfield/Software/PansimNuc/parameter_sweep/baseline/plot"
#tracking_file <- "/Users/samhorsfield/Software/PansimNuc/parameter_sweep/exon_mu_high_recombination_none/tracking.csv"
#outpref <- "/Users/samhorsfield/Software/PansimNuc/parameter_sweep/exon_mu_high_recombination_none/plot"

if (!file.exists(tracking_file)) {
  stop("Cannot find tracking file: ", tracking_file)
}

message("Reading ", tracking_file)
df <- read.csv(tracking_file, stringsAsFactors = FALSE)

required_cols <- c("element_id", "feature_type", "generation",
                   "population_id", "genome_id", "sequence")
missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# For each group (element_id x generation x population_id), split all genome
# sequences into a position matrix and compute nucleotide frequencies and
# mean selection coefficients per nucleotide per position.
compute_position_freqs <- function(seqs, coefficients) {
  seqs         <- seqs[nchar(seqs) > 0]
  coefficients <- coefficients[nchar(coefficients) > 0]
  if (length(seqs) == 0) return(NULL)
  seq_chars   <- strsplit(seqs, "")
  coeff_chars <- strsplit(coefficients, ";")
  len         <- min(lengths(seq_chars), lengths(coeff_chars))
  seq_mat     <- do.call(rbind, lapply(seq_chars,   `[`, seq_len(len)))
  coeff_mat   <- do.call(rbind, lapply(coeff_chars, `[`, seq_len(len)))
  n           <- nrow(seq_mat)
  per_pos <- lapply(seq_len(len), function(p) {
    nucs   <- toupper(seq_mat[, p])
    coeffs <- suppressWarnings(as.numeric(coeff_mat[, p]))
    do.call(rbind, lapply(c("A", "C", "G", "T"), function(nt) {
      idx <- nucs == nt
      data.frame(
        position         = p,
        nucleotide       = nt,
        freq             = sum(idx) / n,
        mean_coefficient = if (any(idx)) mean(coeffs[idx], na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  })
  bind_rows(per_pos)
}

message("Computing per-position nucleotide frequencies...")
freq_data <- df %>%
  group_by(element_id, feature_type, generation, population_id) %>%
  summarise(
    pos_freqs = list(compute_position_freqs(sequence, log_selection_coefficients)),
    .groups   = "drop"
  ) %>%
  filter(!sapply(pos_freqs, is.null)) %>%
  unnest(pos_freqs)

# testing
subset.freq_data <- subset(freq_data, position == 1 & nucleotide == "T")

selection_coeff_data <- freq_data %>%
  group_by(element_id, feature_type, position, nucleotide, population_id) %>%
  summarise(mean_coeff = mean(mean_coefficient, na.rm = TRUE),
            std_dev_coeff = sd(mean_coefficient, na.rm = TRUE))

# determine positions with the greatest variation over time
# ── Select top N positions by the non-start allele that changes the most ─────
# Start allele = most frequent nucleotide at the earliest generation per position
start_alleles <- freq_data %>%
  group_by(element_id, population_id, position) %>%
  filter(generation == min(generation)) %>%
  slice_max(freq, n = 1, with_ties = FALSE) %>%
  select(element_id, population_id, position, start_nucleotide = nucleotide) %>%
  ungroup()

# For every non-start allele compute cumulative absolute frequency change
top_alleles <- freq_data %>%
  inner_join(start_alleles, by = c("element_id", "population_id", "position")) %>%
  filter(nucleotide != start_nucleotide) %>%
  arrange(element_id, population_id, position, nucleotide, generation) %>%
  group_by(element_id, population_id, position, nucleotide) %>%
  summarise(total_change = sum(abs(diff(freq))), .groups = "drop") %>%
  slice_max(total_change, n = top_n_val, with_ties = FALSE)

message(sprintf("Retaining top %d non-start allele(s) by cumulative frequency change.", top_n_val))

# Filter freq_data to only those positions
maf_data <- freq_data %>%
  semi_join(top_alleles, by = c("element_id", "population_id", "position"))

maf_data <- maf_data %>%
  inner_join(selection_coeff_data, by = c("element_id", "feature_type", "population_id", "position", "nucleotide"))

maf_data <- subset(maf_data, select=-c(mean_coefficient))

maf_data$Allele <- paste0("Pos: ", maf_data$position, ", Base: ", maf_data$nucleotide)

n_elements   <- length(unique(maf_data$element_id))
n_pops       <- length(unique(maf_data$population_id))
has_multi_pop <- n_pops > 1

# ── Plot 1: heatmap of minor allele frequency, position × generation ─────────
message("Plotting minor allele frequency heatmap...")

make_label <- function(x) paste0("element_id: ", x)

p_allele_freq <- ggplot(maf_data,
                    aes(x = generation, y = freq,
                        colour   = factor(position),
                        linetype = nucleotide,
                        group    = interaction(position, nucleotide))) +
  geom_line() +
  labs(
    x        = "Generation",
    y        = "Allele frequency",
    colour   = "Position",
    linetype = "Nucleotide"
  ) +
  scale_y_continuous(limits=c(0,1.0)) +
  scale_colour_npg() +
  theme_light(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.text.y  = element_text(size = 7),
    strip.text   = element_text(face = "bold")
  )

if (has_multi_pop) {
  p_allele_freq <- p_allele_freq +
    facet_grid(
      rows = vars(population_id),
      cols = vars(element_id),
      labeller = label_both
    )
} else {
  p_allele_freq <- p_allele_freq +
    facet_wrap(~ element_id, labeller = as_labeller(make_label), ncol = 1)
}

ggsave(paste0(outpref, "_allele_freq_top_", top_n_val, ".pdf"), plot=p_allele_freq, width=8, height=8)

## plot 1b selection coefficient vs. frequency
# plot only for the top alleles
maf_data_top_alleles <- maf_data %>%
  semi_join(top_alleles, by = c("element_id", "population_id", "position", "nucleotide"))

p_selection <- ggplot(maf_data_top_alleles,
                 aes(x = generation, y = freq,
                     colour   = mean_coeff,
                     group    = nucleotide)) +
  facet_grid(. ~ position) +
  geom_line() +
  labs(
    x        = "Generation",
    y        = "Allele frequency",
    colour   = "Selection coefficient",
  ) +
  scale_y_continuous(limits=c(0,1.0)) +
  theme_light(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.text.y  = element_text(size = 7),
    strip.text   = element_text(face = "bold")
  )

ggsave(paste0(outpref, "_allelic_selection_top_", top_n_val, ".pdf"), plot=p_selection, width=3*top_n_val, height=8)



