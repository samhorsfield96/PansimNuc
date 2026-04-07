library(dplyr)
library(ggplot2)
library(tidyr)

# Usage:
#   Rscript plot_allele_frequencies.R [tracking.csv] [output.pdf]
# Defaults to tracking.csv in the current directory and allele_freq_plot.pdf

args          <- commandArgs(trailingOnly = TRUE)
tracking_file <- if (length(args) >= 1) args[1] else "tracking.csv"
output_file   <- if (length(args) >= 2) args[2] else "allele_freq_plot.pdf"

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
# sequences into a position matrix and compute nucleotide frequencies.
compute_position_freqs <- function(seqs) {
  seqs <- seqs[nchar(seqs) > 0]
  if (length(seqs) == 0) return(NULL)
  seq_chars <- strsplit(seqs, "")
  len       <- min(lengths(seq_chars))
  seq_mat   <- do.call(rbind, lapply(seq_chars, `[`, seq_len(len)))
  n         <- nrow(seq_mat)
  per_pos   <- lapply(seq_len(len), function(p) {
    tbl <- table(factor(toupper(seq_mat[, p]), levels = c("A", "C", "G", "T")))
    data.frame(
      position   = p,
      nucleotide = names(tbl),
      freq       = as.numeric(tbl) / n,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(per_pos)
}

message("Computing per-position nucleotide frequencies...")
freq_data <- df %>%
  group_by(element_id, feature_type, generation, population_id) %>%
  summarise(
    pos_freqs = list(compute_position_freqs(sequence)),
    .groups   = "drop"
  ) %>%
  filter(!sapply(pos_freqs, is.null)) %>%
  unnest(pos_freqs)

# Minor allele frequency: 1 - frequency of the most common nucleotide
maf_data <- freq_data %>%
  group_by(element_id, feature_type, generation, population_id, position) %>%
  summarise(
    minor_freq = 1 - max(freq),
    .groups    = "drop"
  )

n_elements   <- length(unique(maf_data$element_id))
n_pops       <- length(unique(maf_data$population_id))
has_multi_pop <- n_pops > 1

# ── Plot 1: heatmap of minor allele frequency, position × generation ─────────
message("Plotting minor allele frequency heatmap...")

make_label <- function(x) paste0("element_id: ", x)

p_heatmap <- ggplot(maf_data,
                    aes(x = position, y = factor(generation), fill = minor_freq)) +
  geom_tile() +
  scale_fill_gradient(
    low  = "white",
    high = "#D7191C",
    name = "Minor\nAllele\nFreq",
    limits = c(0, 1)
  ) +
  labs(
    title    = "Minor Allele Frequency per Position across Generations",
    subtitle = "Each tile is one genomic position; colour = frequency of the non-major allele",
    x        = "Position in Feature (bp)",
    y        = "Generation"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.text.y  = element_text(size = 7),
    strip.text   = element_text(face = "bold")
  )

if (has_multi_pop) {
  p_heatmap <- p_heatmap +
    facet_grid(
      rows = vars(population_id),
      cols = vars(element_id),
      labeller = label_both
    )
} else {
  p_heatmap <- p_heatmap +
    facet_wrap(~ element_id, labeller = as_labeller(make_label), ncol = 1)
}

# ── Plot 2: per-nucleotide frequency for polymorphic positions ────────────────
# Identify positions that are ever polymorphic (MAF > 0)
polymorphic_positions <- maf_data %>%
  group_by(element_id, population_id, position) %>%
  summarise(max_maf = max(minor_freq), .groups = "drop") %>%
  filter(max_maf > 0)

message(sprintf("Found %d polymorphic position(s).", nrow(polymorphic_positions)))

plots <- list(p_heatmap)

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

# ── Save to PDF (one page per plot) ──────────────────────────────────────────
n_positions <- max(maf_data$position, na.rm = TRUE)
heatmap_width  <- max(10, n_positions / 50)
heatmap_height <- max(6,  length(unique(maf_data$generation)) / 5) * n_elements

message("Saving plots to ", output_file)
pdf(output_file, width = heatmap_width, height = heatmap_height)
for (pl in plots) print(pl)
dev.off()

message("Done.")
