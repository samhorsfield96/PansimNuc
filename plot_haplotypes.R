library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(ggpattern)

# Usage:
#   Rscript plot_haplotypes.R [tracking.csv] [output_prefix]
# Defaults to tracking.csv in the current directory and haplotype_plot

args          <- commandArgs(trailingOnly = TRUE)
# Filter out R's own flags (e.g. --no-save, --no-restore) that leak through
args          <- args[!grepl("^--", args)]
tracking_file <- if (length(args) >= 1) args[1] else "tracking.csv"
outpref       <- if (length(args) >= 2) args[2] else "haplotype_analysis"
top_n         <- if (length(args) >= 3) as.integer(args[3]) else 5L  # 0 = keep all

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

# Determine whether h_muts (a parsed mutation vector) is a recombinant of any
# two entries in known_muts_list (a list of pre-parsed mutation vectors).
# A haplotype H is a recombinant of A and B when its mutation set equals exactly
# the union of A's and B's mutation sets, and each parent contributes at least
# one exclusive mutation (neither is a subset of the other).
is_recombinant <- function(h_muts, known_muts_list) {
  if (length(h_muts) == 0 || length(known_muts_list) < 2) return(FALSE)
  h_len <- length(h_muts)
  # Pre-filter: a valid parent must be a non-empty strict subset of h_muts
  candidates <- Filter(
    function(a) length(a) > 0 && length(a) < h_len && all(a %in% h_muts),
    known_muts_list
  )
  if (length(candidates) < 2) return(FALSE)
  n <- length(candidates)
  for (i in seq_len(n - 1)) {
    a <- candidates[[i]]
    for (j in seq(i + 1, n)) {
      b <- candidates[[j]]
      # Both parents must each contribute at least one mutation the other lacks.
      if (!any(!a %in% b) || !any(!b %in% a)) next
      ab <- union(a, b)
      if (length(ab) == h_len && setequal(ab, h_muts)) return(TRUE)
    }
  }
  FALSE
}

# Compute the mean selection coefficient for a set of semicolon-separated
# log-coefficient strings (one string per genome).
# Per genome: sum the per-site log coefficients, get total log selection coefficient.
# Returns the mean of those per-genome values, or NA if unavailable.
hap_sel_coeff <- function(coeff_strs) {
  coeff_strs <- coeff_strs[!is.na(coeff_strs) & nchar(coeff_strs) > 0]
  if (length(coeff_strs) == 0) return(NA_real_)
  per_genome <- sapply(coeff_strs, function(s) {
    vals <- suppressWarnings(as.numeric(strsplit(s, ";")[[1]]))
    sum(vals, na.rm = TRUE)
  })
  mean(per_genome, na.rm = TRUE)
}

# ── Per-group haplotype classification ───────────────────────────────────────
#
# For a single (element_id, feature_type, population_id) group:
#   - Initialise unique haplotypes from generation 1 (type = "founder").
#   - For each subsequent generation, any new sequence is classified as:
#       "recombinant" if its mutation set = union of two prior haplotypes, or
#       "mutant"      otherwise.
#   - Returns a tidy data frame: generation, haplotype_id, mutation_sig, freq,
#     type, sel_coeff
classify_haplotypes <- function(group_df) {
  generations   <- sort(unique(group_df$generation))
  first_gen     <- generations[1]
  has_sel_coeff <- "log_selection_coefficients" %in% colnames(group_df)

  reference <- build_reference(group_df$sequence[group_df$generation == first_gen])

  known_haps  <- list()  # mutation_sig -> "founder" | "mutant" | "recombinant"
  hap_labels  <- list()  # mutation_sig -> short label e.g. "F1", "M2", "R1"
  known_muts  <- list()  # mutation_sig -> pre-parsed mutation vector (cache)
  counter    <- 0L
  new_label  <- function(prefix) { counter <<- counter + 1L; paste0(prefix, counter) }

  rows <- list()

  gen_counter <- 0
  gen <- 1
  pb <- txtProgressBar(min = min(generations), max = max(generations), style = 3)
  message(paste0("Parsing generations: ", max(generations)))
  for (gen in generations) {
    gen_counter <- gen_counter + 1

    gen_mask  <- group_df$generation == gen & nchar(group_df$sequence) > 0
    sel_all   <- if (has_sel_coeff) group_df$log_selection_coefficients[gen_mask] else NULL
    seqs_all  <- group_df$sequence[gen_mask]
    n_total   <- length(seqs_all)
    if (n_total == 0) next

    # Pre-build frequency table and per-sequence coefficient lookup once per gen
    seq_tbl    <- table(seqs_all)
    sel_by_seq <- if (has_sel_coeff) split(sel_all, seqs_all) else NULL

    for (seq in names(seq_tbl)) {
      sig       <- get_mutation_sig(seq, reference)
      sig_muts  <- parse_sig(sig)
      freq      <- seq_tbl[[seq]] / n_total

      # Mean selection coefficient for genomes carrying this haplotype
      if (has_sel_coeff) {
        sel_coeff <- hap_sel_coeff(sel_by_seq[[seq]])
      } else {
        sel_coeff <- NA_real_
      }
      
      # start.time <- Sys.time()
      if (!sig %in% names(known_haps)) {
        if (gen == first_gen) {
          htype <- "founder";     prefix <- "F"
        } else if (is_recombinant(sig_muts, known_muts)) {
          htype <- "recombinant"; prefix <- "R"
        } else {
          htype <- "mutant";      prefix <- "M"
        }
        known_haps[[sig]] <- htype
        hap_labels[[sig]] <- new_label(prefix)
        known_muts[[sig]] <- sig_muts  # cache parsed form
      }
      # end.time <- Sys.time()
      # time.taken <- end.time - start.time
      # time.taken
      # message(paste0("Time taken: ", time.taken))
      
      if (sig == "") {
        haplotype_id = "R"
        htype = "reference"
      } else {
        haplotype_id = hap_labels[[sig]]
        htype         = known_haps[[sig]]
      }

      df_tmp <- data.frame(
        generation     = gen,
        haplotype_id   = haplotype_id,
        mutation_sig   = sig,
        freq           = freq,
        type           = htype,
        sel_coeff = sel_coeff,
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
  fill(type, mutation_sig, sel_coeff, .direction = "downup") %>%
  ungroup()

# ── Summary table ─────────────────────────────────────────────────────────────
hap_summary <- hap_data %>%
  group_by(element_id, feature_type, population_id, haplotype_id, type, mutation_sig) %>%
  summarise(
    first_generation = min(generation[freq > 0]),
    peak_freq        = max(freq),
    mean_sel_coeff   = mean(sel_coeff, na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  arrange(element_id, population_id, first_generation)

write.csv(hap_summary,
          file      = paste0(outpref, "_haplotype_summary.csv"),
          row.names = FALSE)

# ── Filter to top N haplotypes per type (0 = keep all) ───────────────────────
if (top_n > 0L) {
  message(sprintf("Retaining top %d haplotype(s) per type by cumulative frequency change.", top_n))
  top_haps <- hap_data %>%
    arrange(element_id, feature_type, population_id, haplotype_id, generation) %>%
    group_by(element_id, feature_type, population_id, haplotype_id, type) %>%
    summarise(total_change = sum(abs(diff(freq))), .groups = "drop") %>%
    group_by(element_id, feature_type, population_id, type) %>%
    slice_max(total_change, n = top_n, with_ties = FALSE) %>%
    ungroup()

  hap_data <- hap_data %>%
    semi_join(top_haps, by = c("element_id", "feature_type", "population_id", "haplotype_id"))
}

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
  labs(x = "Generation", y = "Haplotype frequency", colour = "Haplotype") +
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
  labs(x = "Generation", y = "Cumulative haplotype frequency", fill = "Haplotype") +
  scale_y_continuous(limits = c(0, 1)) +
  type_fill_scale +
  base_theme

p_area <- add_facets(p_area)
p_area
ggsave(paste0(outpref, "_haplotype_composition.pdf"), plot = p_area, width = 8, height = 6)

# ── Plot 3: stacked area chart of top changing haplotype composition ──────────────────────
message("Plotting haplotype composition stacked areas...")

p_area <- ggplot(
  hap_data,
  aes(
    x    = generation,
    y    = freq,
    fill = haplotype_id,
    group = haplotype_id
  )
) +
  geom_area(position = "stack", colour = NA, alpha = 0.8) +
  labs(x = "Generation", y = "Cumulative haplotype frequency", fill = "Haplotype ID") +
  scale_y_continuous(limits = c(0, 1)) +
  base_theme

p_area <- add_facets(p_area)
p_area
ggsave(paste0(outpref, "_per_haplotype_composition.pdf"), plot = p_area, width = 8, height = 6)

# ── Plot 4: top hits with selection coefficients + haplotype-type hatching ───
type_pattern_values <- c(
  reference   = "none",
  founder     = "none",
  mutant      = "stripe",
  recombinant = "crosshatch"
)

p_sel <- ggplot(
  hap_data,
  aes(
    x              = generation,
    y              = freq,
    fill           = sel_coeff,
    colour           = sel_coeff,
    pattern        = type,
    pattern_colour = type,
    group          = haplotype_id
  )
) +
  geom_area_pattern(
    position        = "stack",
    colour          = NA,
    alpha           = 0.8,
    pattern_density = 0.35,
    pattern_spacing = 0.025,
    pattern_fill    = NA
  ) +
  scale_pattern_manual(
    values = type_pattern_values,
    name   = "Haplotype type"
  ) +
  scale_pattern_colour_manual(
    values = c(
      reference   = "grey30",
      founder     = "grey30",
      mutant      = "grey30",
      recombinant = "grey30"
    ),
    name = "Haplotype type"
  ) +
  labs(x = "Generation", y = "Cumulative haplotype frequency", fill = "Selection coefficient") +
  scale_y_continuous(limits = c(0, 1)) +
  base_theme

p_sel <- add_facets(p_sel)
p_sel
ggsave(paste0(outpref, "_sel_coeff_composition.pdf"), plot = p_sel, width = 8, height = 6)

message(sprintf(
  "Done. %d haplotypes tracked (%d founder, %d mutant, %d recombinant).",
  nrow(hap_summary),
  sum(hap_summary$type == "founder"),
  sum(hap_summary$type == "mutant"),
  sum(hap_summary$type == "recombinant")
))

