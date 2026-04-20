library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)

# Usage:
#   Rscript plot_haplotypes.R [tracking.csv] [output_prefix]
# Defaults to tracking.csv in the current directory and haplotype_plot

args          <- commandArgs(trailingOnly = TRUE)
# Filter out R's own flags (e.g. --no-save, --no-restore) that leak through
args          <- args[!grepl("^--", args)]
tracking_file <- if (length(args) >= 1) args[1] else "tracking.csv"
outpref       <- if (length(args) >= 2) args[2] else "haplotype_plot"

tracking_file <- "/Users/samhorsfield/Library/CloudStorage/OneDrive-Personal/Work/Postdoc_Unine/Analysis/PansimNuc/parameter_sweep/baseline/tracking.csv"

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

# ── Helper functions ──────────────────────────────────────────────────────────

# Build reference sequence as the majority nucleotide at each position,
# computed from the sequences in the founding generation.
build_reference <- function(seqs) {
  valid <- seqs[nchar(seqs) > 0]
  if (length(valid) == 0) return(character(0))
  chars <- strsplit(valid, "")
  len   <- min(lengths(chars))
  mat   <- do.call(rbind, lapply(chars, `[`, seq_len(len)))
  apply(mat, 2, function(col) {
    tb <- table(toupper(col))
    names(tb)[which.max(tb)]
  })
}

# Compute a mutation signature string "pos:base;pos:base;..." for a sequence
# relative to the reference. Returns "" for sequences identical to the reference.
get_mutation_sig <- function(seq, reference) {
  chars <- toupper(strsplit(seq, "")[[1]])
  len   <- min(length(chars), length(reference))
  idx   <- which(chars[seq_len(len)] != reference[seq_len(len)])
  if (length(idx) == 0) return("")
  paste(paste0(idx, ":", chars[idx]), collapse = ";")
}

# Parse a mutation signature string into a character vector of "pos:base" tokens.
parse_sig <- function(s) if (nchar(s) == 0) character(0) else strsplit(s, ";")[[1]]

# Determine whether sig_h is a recombinant of any two signatures in known_sigs.
# A haplotype H is a recombinant of A and B when its mutation set equals exactly
# the union of A's and B's mutation sets, and each parent contributes at least
# one exclusive mutation (neither is a subset of the other).
is_recombinant <- function(sig_h, known_sigs) {
  h_muts <- parse_sig(sig_h)
  # A reference-identical haplotype or fewer than two known haplotypes: skip.
  if (length(h_muts) == 0 || length(known_sigs) < 2) return(FALSE)
  known_muts <- lapply(known_sigs, parse_sig)
  n <- length(known_muts)
  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      a <- known_muts[[i]]
      b <- known_muts[[j]]
      # Both parents must each contribute at least one mutation the other lacks.
      if (length(setdiff(a, b)) == 0 || length(setdiff(b, a)) == 0) next
      if (setequal(h_muts, union(a, b))) return(TRUE)
    }
  }
  FALSE
}

# ── Per-group haplotype classification ───────────────────────────────────────
#
# For a single (element_id, feature_type, population_id) group:
#   - Initialise unique haplotypes from generation 1 (type = "founder").
#   - For each subsequent generation, any new sequence is classified as:
#       "recombinant" if its mutation set = union of two prior haplotypes, or
#       "mutant"      otherwise.
#   - Returns a tidy data frame: generation, haplotype_id, mutation_sig, freq, type.
classify_haplotypes <- function(group_df) {
  generations <- sort(unique(group_df$generation))
  first_gen   <- generations[1]

  reference <- build_reference(group_df$sequence[group_df$generation == first_gen])

  known_haps <- list()  # mutation_sig -> "founder" | "mutant" | "recombinant"
  hap_labels <- list()  # mutation_sig -> short label e.g. "F1", "M2", "R1"
  counter    <- 0L
  new_label  <- function(prefix) { counter <<- counter + 1L; paste0(prefix, counter) }

  rows <- list()

  gen_counter <- 0
  gen <- 1
  pb <- txtProgressBar(min = min(generations), max = max(generations), style = 3)
  message(paste0("Parsing generations: ", max(generations)))
  for (gen in generations) {
    gen_counter <- gen_counter + 1
    seqs    <- group_df$sequence[group_df$generation == gen]
    seqs    <- seqs[nchar(seqs) > 0]
    n_total <- length(seqs)
    if (n_total == 0) next

    seq_counter <- 0
    for (seq in unique(seqs)) {
      seq_counter <- seq_counter + 1
      sig  <- get_mutation_sig(seq, reference)
      freq <- sum(seqs == seq) / n_total

      if (!sig %in% names(known_haps)) {
        if (gen == first_gen) {
          htype <- "founder";     prefix <- "F"
        } else if (is_recombinant(sig, names(known_haps))) {
          htype <- "recombinant"; prefix <- "R"
        } else {
          htype <- "mutant";      prefix <- "M"
        }
        known_haps[[sig]] <- htype
        hap_labels[[sig]] <- new_label(prefix)
      }
      
      if (sig == "") {
        haplotype_id = "R"
        htype = "reference"
      } else {
        haplotype_id = hap_labels[[sig]]
        htype         = known_haps[[sig]]
      }

      df_tmp <- data.frame(
        generation   = gen,
        haplotype_id = haplotype_id,
        mutation_sig = sig,
        freq         = freq,
        type         = htype,
        stringsAsFactors = FALSE
      )
      
      rows[[length(rows) + 1]] <- df_tmp
    }
    setTxtProgressBar(pb, gen)
  }
  close(pb)

  bind_rows(rows)
}

# ── Run classification across all groups ─────────────────────────────────────
message("Tracking haplotype frequencies...")

group_df <- subset(df, element_id == "402" & population_id == 0)

hap_data <- df %>%
  group_by(element_id, feature_type, population_id) %>%
  group_modify(~ classify_haplotypes(.x)) %>%
  ungroup()

# Fill in zero-frequency rows so that every haplotype appears in every
# generation (needed for correct stacked areas and unbroken lines).
hap_data <- hap_data %>%
  group_by(element_id, feature_type, population_id) %>%
  complete(
    generation   = unique(generation),
    haplotype_id = unique(haplotype_id),
    fill         = list(freq = 0)
  ) %>%
  group_by(element_id, feature_type, population_id, haplotype_id) %>%
  fill(type, mutation_sig, .direction = "downup") %>%
  ungroup()

# ── Plotting helpers ──────────────────────────────────────────────────────────
n_pops        <- length(unique(hap_data$population_id))
has_multi_pop <- n_pops > 1

type_colour_values <- c(
  reference = "#3C5488FF",
  founder     = "#4DBBD5",
  mutant      = "#E64B35",
  recombinant = "#00A087"
)

type_colour_scale <- scale_colour_manual(values = type_colour_values, name = "Haplotype")
type_fill_scale   <- scale_fill_manual(values = type_colour_values,   name = "Haplotype")

add_facets <- function(p) {
  if (has_multi_pop) {
    p + facet_grid(
      rows     = vars(population_id),
      cols     = vars(element_id),
      labeller = label_both
    )
  } else {
    p + facet_wrap(
      ~ element_id,
      labeller = as_labeller(function(x) paste0("element_id: ", x)),
      ncol     = 1
    )
  }
}

base_theme <- theme_light(base_size = 11) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"))

# ── Plot 1: haplotype frequency lines, coloured by type ──────────────────────
message("Plotting haplotype frequency lines...")

p_lines <- ggplot(
  hap_data,
  aes(
    x      = generation,
    y      = freq,
    colour = type,
    group  = interaction(haplotype_id, population_id)
  )
) +
  geom_line(alpha = 0.8) +
  labs(x = "Generation", y = "Haplotype frequency") +
  scale_y_continuous(limits = c(0, 1)) +
  type_colour_scale +
  base_theme

p_lines <- add_facets(p_lines)
p_lines
ggsave(paste0(outpref, "_haplotype_freq.pdf"), plot = p_lines, width = 8, height = 6)

# ── Plot 2: stacked area chart of haplotype composition ──────────────────────
message("Plotting haplotype composition stacked areas...")

p_area <- ggplot(
  hap_data,
  aes(
    x    = generation,
    y    = freq,
    fill = type,
    group = haplotype_id
  )
) +
  geom_area(position = "stack", colour = NA, alpha = 0.8) +
  labs(x = "Generation", y = "Cumulative haplotype frequency") +
  scale_y_continuous(limits = c(0, 1)) +
  type_fill_scale +
  base_theme

p_area <- add_facets(p_area)
p_area
ggsave(paste0(outpref, "_haplotype_composition.pdf"), plot = p_area, width = 8, height = 6)

# ── Summary table ─────────────────────────────────────────────────────────────
hap_summary <- hap_data %>%
  group_by(element_id, feature_type, population_id, haplotype_id, type, mutation_sig) %>%
  summarise(
    first_generation = min(generation[freq > 0]),
    peak_freq        = max(freq),
    .groups          = "drop"
  ) %>%
  arrange(element_id, population_id, first_generation)

write.csv(hap_summary,
          file      = paste0(outpref, "_haplotype_summary.csv"),
          row.names = FALSE)

message(sprintf(
  "Done. %d haplotypes tracked (%d founder, %d mutant, %d recombinant).",
  nrow(hap_summary),
  sum(hap_summary$type == "founder"),
  sum(hap_summary$type == "mutant"),
  sum(hap_summary$type == "recombinant")
))
message("Output files:")
message("  ", outpref, "_haplotype_freq.pdf")
message("  ", outpref, "_haplotype_composition.pdf")
message("  ", outpref, "_haplotype_summary.csv")
