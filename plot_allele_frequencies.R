library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)

# Usage:
#   Rscript plot_allele_frequencies.R [tracking.csv] [output.pdf]
# Defaults to tracking.csv in the current directory and allele_freq_plot.pdf

args          <- commandArgs(trailingOnly = TRUE)
tracking_file <- if (length(args) >= 1) args[1] else "tracking.csv"
output_file   <- if (length(args) >= 2) args[2] else "allele_freq_plot.pdf"
top_n_val <- if (length(args) >= 3) as.numeric(args[3]) else 1

tracking_file <- "/Users/samhorsfield/Software/PansimNuc/parameter_sweep/baseline/tracking.csv"
output_file <- "/Users/samhorsfield/Software/PansimNuc/parameter_sweep/baseline/tracking_plot.pdf"

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

p_allele_freq

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

p_selection

# ── Plot 2: per-nucleotide frequency for polymorphic positions ────────────────
# Identify positions that are ever polymorphic (MAF > 0)
polymorphic_positions <- maf_data %>%
  group_by(element_id, population_id, position) %>%
  summarise(max_maf = max(minor_freq), .groups = "drop") %>%
  filter(max_maf > 0)

message(sprintf("Found %d polymorphic position(s).", nrow(polymorphic_positions)))

plots <- list(p_allele_freq)

if (nrow(polymorphic_positions) > 0) {
  poly_freq <- freq_data %>%
    semi_join(polymorphic_positions,
              by = c("element_id", "population_id", "position")) %>%
    mutate(pos_label = paste0("pos ", position))

  p_poly <- ggplot(poly_freq,
                   aes(x = generation, y = freq,
                       colour = nucleotide, group = nucleotide)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.2) +
    scale_colour_manual(
      values = c(A = "#2166AC", C = "#4DAC26", G = "#D6604D", T = "#762A83"),
      name   = "Nucleotide"
    ) +
    scale_y_continuous(limits = c(0, 1), name = "Frequency") +
    scale_x_continuous(name = "Generation") +
    labs(
      title    = "Nucleotide Frequencies at Polymorphic Positions",
      subtitle = "Only positions with MAF > 0 in at least one generation are shown"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text       = element_text(size = 7, face = "bold")
    )

  facet_vars <- c("element_id", "pos_label")
  if (has_multi_pop) facet_vars <- c("population_id", facet_vars)
  p_poly <- p_poly + facet_wrap(facet_vars, ncol = 5)

  plots <- c(plots, list(p_poly))
}

# plot 3 plot selection coefficient vs. frequency vs generation

# ── Save to PDF (one page per plot) ──────────────────────────────────────────
n_positions <- max(maf_data$position, na.rm = TRUE)
heatmap_width  <- max(10, n_positions / 50)
heatmap_height <- max(6,  length(unique(maf_data$generation)) / 5) * n_elements

message("Saving plots to ", output_file)
pdf(output_file, width = heatmap_width, height = heatmap_height)
for (pl in plots) print(pl)
dev.off()

message("Done.")

