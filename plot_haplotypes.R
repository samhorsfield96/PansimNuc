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

# ── Whole-genome recombinant detection ────────────────────────────────────────
# profile_c      : named character vector  element_key -> mut_sig  for one genome
# known_profiles : list of such named vectors for all previously-seen genome haplotypes
# C is a recombinant of A and B when every element allele in C matches A or B,
# A contributes >=1 element that differs from B, and B contributes >=1 that differs from A.
is_recombinant_genome <- function(profile_c, known_profiles) {
  if (length(known_profiles) < 2) return(FALSE)
  elements <- names(profile_c)
  n <- length(known_profiles)
  for (i in seq_len(n - 1)) {
    a <- known_profiles[[i]]
    for (j in seq(i + 1, n)) {
      b <- known_profiles[[j]]
      common_els <- Reduce(intersect, list(elements, names(a), names(b)))
      if (length(common_els) < 2) next
      c_sub <- profile_c[common_els]
      a_sub <- a[common_els]
      b_sub <- b[common_els]
      matches_a <- c_sub == a_sub
      matches_b <- c_sub == b_sub
      if (!all(matches_a | matches_b)) next
      if (any(matches_a & !matches_b) && any(matches_b & !matches_a)) return(TRUE)
    }
  }
  FALSE
}

# ── Assign per-element mutation signatures ────────────────────────────────────
# For a single (element_id, feature_type, population_id) group, build the
# reference from the first generation and annotate every row with its
# mutation signature relative to that reference.
assign_element_sigs <- function(group_df) {
  first_gen <- min(group_df$generation)
  reference <- build_reference(group_df$sequence[group_df$generation == first_gen])
  group_df$mut_sig <- vapply(group_df$sequence, function(s) {
    if (is.na(s) || nchar(s) == 0) return(NA_character_)
    get_mutation_sig(s, reference)
  }, character(1L))
  group_df
}

# ── Classify whole-genome haplotypes for one population ──────────────────────
# pop_df has per-element mutation signatures (mut_sig) for all genomes in one
# population across all generations.
# Returns a tidy data frame: generation, haplotype_id, profile_str, freq, type, sel_coeff
classify_genome_haplotypes <- function(pop_df) {
  has_sel     <- "log_selection_coefficients" %in% colnames(pop_df)
  generations <- sort(unique(pop_df$generation))
  first_gen   <- generations[1]

  # Build per-genome profile: a named vector (element_key -> mut_sig) per (genome_id, gen)
  genome_profiles <- pop_df %>%
    filter(!is.na(mut_sig)) %>%
    mutate(element_key = paste(element_id, feature_type, sep = ":")) %>%
    group_by(genome_id, generation) %>%
    summarise(
      profile_str = paste(sort(paste(element_key, mut_sig, sep = "=")), collapse = "|"),
      profile_vec = list(setNames(mut_sig, element_key)),
      sel_coeff   = if (has_sel) {
        sum(vapply(log_selection_coefficients, function(s) {
          vals <- suppressWarnings(as.numeric(strsplit(s, ";")[[1]]))
          sum(vals, na.rm = TRUE)
        }, numeric(1L)), na.rm = TRUE)
      } else NA_real_,
      .groups = "drop"
    )

  known_profiles <- list()   # profile_str -> named vector
  known_types    <- list()   # profile_str -> type string
  hap_labels     <- list()   # profile_str -> short label
  counter        <- 0L
  new_label <- function(prefix) { counter <<- counter + 1L; paste0(prefix, counter) }

  rows <- list()
  pb <- txtProgressBar(min = min(generations), max = max(generations), style = 3)
  message(paste0("Classifying genome haplotypes across ", length(generations), " generations"))

  for (gen in generations) {
    gen_data <- genome_profiles[genome_profiles$generation == gen, ]
    n_total  <- nrow(gen_data)
    if (n_total == 0) next

    prof_tbl    <- table(gen_data$profile_str)
    sel_by_prof <- split(gen_data$sel_coeff, gen_data$profile_str)
    vec_by_prof <- lapply(
      split(seq_len(nrow(gen_data)), gen_data$profile_str),
      function(idx) gen_data$profile_vec[[idx[1]]]
    )

    for (prof_str in names(prof_tbl)) {
      freq      <- prof_tbl[[prof_str]] / n_total
      sel_coeff <- mean(unlist(sel_by_prof[[prof_str]]), na.rm = TRUE)
      prof_vec  <- vec_by_prof[[prof_str]]

      if (!prof_str %in% names(known_profiles)) {
        if (gen == first_gen) {
          htype <- "founder"; prefix <- "F"
        } else if (is_recombinant_genome(prof_vec, unname(known_profiles))) {
          htype <- "recombinant"; prefix <- "R"
        } else {
          htype <- "mutant"; prefix <- "M"
        }
        known_profiles[[prof_str]] <- prof_vec
        known_types[[prof_str]]    <- htype
        hap_labels[[prof_str]]     <- new_label(prefix)
      }

      rows[[length(rows) + 1]] <- data.frame(
        generation   = gen,
        haplotype_id = hap_labels[[prof_str]],
        profile_str  = prof_str,
        freq         = freq,
        type         = known_types[[prof_str]],
        sel_coeff    = sel_coeff,
        stringsAsFactors = FALSE
      )
    }
    setTxtProgressBar(pb, gen)
  }
  close(pb)
  bind_rows(rows)
}

# ── Run classification ────────────────────────────────────────────────────────
message("Assigning per-element mutation signatures...")

element_sig_df <- df %>%
  group_by(element_id, feature_type, population_id) %>%
  group_modify(~ assign_element_sigs(.x)) %>%
  ungroup()

message("Classifying whole-genome haplotypes per population...")

hap_data <- element_sig_df %>%
  group_by(population_id) %>%
  group_modify(~ classify_genome_haplotypes(.x)) %>%
  ungroup()

# Fill in zero-frequency rows so that every haplotype appears in every
# generation (needed for correct stacked areas and unbroken lines).
hap_data <- hap_data %>%
  group_by(population_id) %>%
  complete(
    generation   = unique(generation),
    haplotype_id = unique(haplotype_id),
    fill         = list(freq = 0)
  ) %>%
  group_by(population_id, haplotype_id) %>%
  fill(type, profile_str, sel_coeff, .direction = "downup") %>%
  ungroup()

# ── Summary table ─────────────────────────────────────────────────────────────
hap_summary <- hap_data %>%
  group_by(population_id, haplotype_id, type, profile_str) %>%
  summarise(
    first_generation = min(generation[freq > 0]),
    peak_freq        = max(freq),
    mean_sel_coeff   = mean(sel_coeff, na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  arrange(population_id, first_generation)

write.csv(hap_summary,
          file      = paste0(outpref, "_haplotype_summary.csv"),
          row.names = FALSE)

# ── Filter to top N haplotypes per type (0 = keep all) ───────────────────────
if (top_n > 0L) {
  message(sprintf("Retaining top %d haplotype(s) per type by cumulative frequency change.", top_n))
  top_haps <- hap_data %>%
    arrange(population_id, haplotype_id, generation) %>%
    group_by(population_id, haplotype_id, type) %>%
    summarise(total_change = sum(abs(diff(freq))), .groups = "drop") %>%
    group_by(population_id, type) %>%
    slice_max(total_change, n = top_n, with_ties = FALSE) %>%
    ungroup()

  hap_data <- hap_data %>%
    semi_join(top_haps, by = c("population_id", "haplotype_id"))
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
    p + facet_wrap(~ population_id, labeller = label_both, ncol = 1)
  } else {
    p
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
    group  = haplotype_id
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

