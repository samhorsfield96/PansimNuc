library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)
library(ggpattern)

# Usage:
#   Rscript plot_haplotype_network.R [gff_dir] [output_prefix] [generation]
#
# Arguments:
#   gff_dir                 – directory containing GFF + FASTA files (default: .)
#   output_prefix           – prefix for output files       (default: haplotypes)
#   top_n                   - Number of top-changing haplotypes to return. 0 keeps all.
#  recombinant_threshold    - proportion of mutations in a haplotype that must be present to 
#                             identify a new haplotype as being a recombinant. Default 0.9 (90%). 
#                             Set to 0 to disable recombinant detection and classification.
#
# GFF files must match: pop_<pop>_gen_<gen>_genome_<id>.gff
# FASTA files must share the same basename with a .fasta extension.
#
# Output: multiple PDFs per population_id group, describing the changes in haplotypes
#         across generations.

args    <- commandArgs(trailingOnly = TRUE)
args    <- args[!grepl("^--", args)]
gff_dir <- if (length(args) >= 1) args[1] else "."
outpref <- if (length(args) >= 2) args[2] else "haplotypes"
top_n   <- if (length(args) >= 3) as.integer(args[3]) else 5L  # 0 = keep all
recombination_threshold <- if (length(args) >= 4) as.numeric(args[4]) else 0.9  # 0 = keep all

# ── GFF / FASTA reading (adapted from ld_analysis.R) ─────────────────────────

parse_attrs <- function(attr_str) {
  pairs <- strsplit(attr_str, ";", fixed = TRUE)[[1L]]
  kv    <- strsplit(pairs, "=", fixed = TRUE)
  keys  <- vapply(kv, `[[`, character(1L), 1L)
  vals  <- vapply(kv, function(x) if (length(x) >= 2L) x[[2L]] else NA_character_,
                  character(1L))
  setNames(vals, keys)
}

# ── GFF reader ────────────────────────────────────────────────────────────────
# Returns a data.frame with one row per feature (including selection coeff).
# checked, all good
read_sim_gff <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nchar(lines) > 0L & !startsWith(lines, "#")]
  if (length(lines) == 0L) return(NULL)
  rows <- lapply(lines, function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(f) < 9L) return(NULL)
    a <- parse_attrs(f[9L])
    # contig_N in col-1 → 0-based index for FASTA lookup
    contig_name  <- f[1L]
    contig_index <- suppressWarnings(
      as.integer(sub("contig_", "", contig_name)) - 1L
    )
    data.frame(
      contig_index = contig_index,
      start        = as.integer(f[4L]),   # GFF is 1-based
      end          = as.integer(f[5L]),
      strand       = f[7L],
      element_id   = suppressWarnings(as.integer(a[["element_id"]])),
      feature_type = a[["feature_type"]],
      log_sel_coeff = suppressWarnings(as.numeric(a[["log_element_selection_coefficient"]])),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(Filter(Negate(is.null), rows))
}

# ── FASTA reader ──────────────────────────────────────────────────────────────
# checked, all good
read_fasta <- function(path) {
  lines   <- readLines(path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  seqs    <- vector("list", length(headers))
  for (i in seq_along(headers)) {
    h_line  <- lines[headers[i]]
    # Extract suffix _contig<N> from the header
    m <- regmatches(h_line, regexpr("_contig(\\d+)", h_line, perl = TRUE))
    if (length(m) == 0L) next
    idx <- as.integer(sub("_contig", "", m))
    body_start <- headers[i] + 1L
    body_end   <- if (i < length(headers)) headers[i + 1L] - 1L else length(lines)
    seq_str    <- paste(lines[body_start:body_end], collapse = "")
    seqs[[i]]  <- list(idx = idx, seq = toupper(seq_str))
  }
  result <- Filter(Negate(is.null), seqs)
  setNames(
    vapply(result, `[[`, character(1L), "seq"),
    vapply(result, function(x) as.character(x[["idx"]]), character(1L))
  )
}

# ── Reverse complement ────────────────────────────────────────────────────────
# checked, all good
rev_comp <- function(seq) {
  comp <- chartr("ACGTN", "TGCAN", seq)
  paste(rev(strsplit(comp, "")[[1L]]), collapse = "")
}

# checked, all good
extract_window <- function(fasta_seqs, contig_index, start, end, strand) {
  key <- as.character(contig_index)
  if (!key %in% names(fasta_seqs)) return(NA_character_)
  full <- fasta_seqs[[key]]
  start <- max(1L, start)
  end   <- min(nchar(full), end)
  if (start > end) return(NA_character_)
  sub_seq <- substr(full, start, end)
  if (strand == "-") sub_seq <- rev_comp(sub_seq)
  sub_seq
}

# ── Discover and parse GFF files ──────────────────────────────────────────────

gff_files <- list.files(gff_dir,
                        pattern    = "^pop_\\d+_gen_\\d+_genome_\\d+\\.gff$",
                        full.names = TRUE)
if (length(gff_files) == 0L) {
  stop("No GFF files matching pop_<pop>_gen_<gen>_genome_<id>.gff found in: ", gff_dir)
}
message("Found ", length(gff_files), " GFF file(s) in: ", gff_dir)

file_meta <- lapply(gff_files, function(fp) {
  bn <- sub("\\.gff$", "", basename(fp))
  m  <- regmatches(bn, regexpr("^pop_(\\d+)_gen_(\\d+)_genome_(\\d+)$", bn, perl = TRUE))
  if (length(m) == 0L) return(NULL)
  parts <- as.integer(strsplit(sub("^pop_", "", m), "_gen_|_genome_")[[1L]])
  data.frame(
    gff_path   = fp,
    fasta_path = sub("\\.gff$", ".fasta", fp),
    pop_id     = parts[1L],
    gen_id     = parts[2L],
    genome_id  = parts[3L],
    stringsAsFactors = FALSE
  )
})
file_meta <- bind_rows(Filter(Negate(is.null), file_meta))

# ── Extract element sequences from each genome ────────────────────────────────
message("Extracting element sequences from GFF + FASTA files...")

# read existing dataset if present
gff_df_rds_file <- paste0(outpref, "_gff_df.rds")
if (!file.exists(gff_df_rds_file)) {
  df_rows <- lapply(seq_len(nrow(file_meta)), function(i) {
    row <- file_meta[i, ]
    if (!file.exists(row$fasta_path)) {
      message("  FASTA not found, skipping: ", row$fasta_path)
      return(NULL)
    }
    gff   <- read_sim_gff(row$gff_path)
    if (is.null(gff) || nrow(gff) == 0L) return(NULL)
    fasta <- read_fasta(row$fasta_path)
    
    gff$sequence <- mapply(
      extract_window,
      contig_index = gff$contig_index,
      start        = gff$start,
      end          = gff$end,
      strand       = gff$strand,
      MoreArgs     = list(fasta_seqs = fasta)
    )
    
    gff$generation   <- row$gen_id
    gff$population_id <- row$pop_id
    gff$genome_id    <- row$genome_id
    gff[!is.na(gff$sequence) & !is.na(gff$element_id), ]
  })
  
  df <- bind_rows(Filter(Negate(is.null), df_rows))
  saveRDS(df, gff_df_rds_file)
} else {
  df <- readRDS(gff_df_rds_file)
}

if (nrow(df) == 0L) stop("No element sequences could be extracted.")
message(sprintf("Extracted sequences for %d element × genome records.", nrow(df)))

# ── Helper functions ──────────────────────────────────────────────────────────
# checked, all good
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

# determine mutated sites between reference and sequences
# checked, all good
get_mutation_sig <- function(seq, reference, element_id) {
  chars <- toupper(strsplit(seq, "")[[1]])
  len   <- min(length(chars), length(reference))
  idx   <- which(chars[seq_len(len)] != reference[seq_len(len)])
  if (length(idx) == 0) return(NA_character_)
  paste(paste0(element_id, ":", idx, ":", chars[idx]), collapse = ";")
}

# checked, all good
parse_sig <- function(s) if (nchar(s) == 0) character(0) else strsplit(s, ";")[[1]]

# ── Whole-genome haplotype functions ─────────────────────────────────────────

# Assign per-element mutation signatures relative to each element's founding
# generation consensus. Called per (element_id, feature_type, population_id).
# checked, all good
assign_element_sigs <- function(group_df) {
  first_gen <- min(group_df$generation)
  reference <- build_reference(group_df$sequence[group_df$generation == first_gen])
  element_id <- unique(group_df$element_id_tmp)
  group_df$mut_sig <- vapply(group_df$sequence, function(s) {
    get_mutation_sig(s, reference, element_id)
  }, character(1L))
  group_df
}

# Genome-level recombinant detection.
# profile_c           : parsed character vector of mutation tokens for one genome
# all_profiles_named  : named list (profile_str -> parsed vec) of all known profiles
# threshold           : proportion of mutations in a haplotype that must be present to 
#                       identify a new haplotype as being a recombinant. Default 0.9 (90%). 
#                       Set to 0 to disable recombinant detection and classification.
# Returns a length-2 character vector of the two parent profile strings when C
# is a recombinant (union of their mutation sets equals C's, each contributes
# at least one exclusive mutation), or NULL otherwise.
# checked, all good
find_recombinant_parents <- function(profile_c, all_profiles_named, threshold = 0.9) {
  if (threshold <= 0 || length(all_profiles_named) < 2) return(NULL)

  profile_c_len <- length(profile_c)
  candidate_names <- Filter(
    function(nm) {
      a <- all_profiles_named[[nm]]
      length(a) > 0 &&
        (sum(a %in% profile_c) / length(a)) >= threshold
    },
    names(all_profiles_named)
  )
  if (length(candidate_names) < 2) return(NULL)

  n <- length(candidate_names)
  for (i in seq_len(n - 1)) {
    a_name <- candidate_names[[i]]
    a      <- all_profiles_named[[a_name]]
    for (j in seq(i + 1, n)) {
      b_name <- candidate_names[[j]]
      b      <- all_profiles_named[[b_name]]
      if (!any(!a %in% b) || !any(!b %in% a)) next
      # A and B each contribute >= threshold proportion of their mutations to C
      return(c(a_name, b_name))
    }
  }
  NULL
}

# Classify whole-genome haplotypes for one population.
# pop_df must have mut_sig (from assign_element_sigs) and log_sel_coeff columns.
# Returns: generation, haplotype_id, profile_str, sequence (concatenated), freq, type, sel_coeff
# checked, all good
classify_genome_haplotypes <- function(pop_df) {
  generations <- sort(unique(pop_df$generation))

  genome_profiles <- pop_df %>%
    group_by(genome_id, generation) %>%
    summarise(
      profile_str = {
        sigs <- mut_sig[!is.na(mut_sig)]
        if (length(sigs) == 0L) "NA" else paste(sigs, collapse = ";")
      },
      sel_coeff   = sum(log_sel_coeff, na.rm = TRUE),
      .groups     = "drop"
    )
  
  # Build a named list of all profiles (across all generations) for recombinant
  # detection – not reliant on the order in which generations are processed.

  known_profiles <- list()   # profile_str -> parsed vec  (profiles seen so far)
  known_types    <- list()   # profile_str -> haplotype type string
  known_parents  <- list()   # profile_str -> "P1,P2" or NA
  hap_labels     <- list()   # profile_str -> short label
  counter        <- 0L
  new_label <- function(prefix) { counter <<- counter + 1L; paste0(prefix, counter) }
  
  # determine how many generations present, adjust which generation to look for recombinants
  if (length(generations) > 1)
  {
    adjustment = 1
  } else {
    adjustment = 0
  }

  rows <- list()
  for (gen in generations) {
    gen_data <- genome_profiles[genome_profiles$generation == gen, ]
    n_total  <- nrow(gen_data)
    if (n_total == 0) next

    # look for recombinants in the context of profiles only in prior generation
    all_profile_names  <- names(table(genome_profiles[genome_profiles$generation == (gen - adjustment), ]$profile_str))
    all_profiles_named <- setNames(
      lapply(all_profile_names, parse_sig),
      all_profile_names
    )

    prof_tbl    <- table(gen_data$profile_str)
    sel_by_prof <- split(gen_data$sel_coeff, gen_data$profile_str)

    for (prof_str in names(prof_tbl)) {
      freq      <- prof_tbl[[prof_str]] / n_total
      sel_coeff <- mean(unlist(sel_by_prof[[prof_str]]), na.rm = TRUE)

      if (!prof_str %in% names(known_profiles)) {
        prof_vec       <- parse_sig(prof_str)
        parent_str     <- NA_character_
        if (prof_str == "NA") {
          htype <- "reference"; prefix <- "REF"
        } else if (gen == 0) {
          htype <- "founder"; prefix <- "F"
        } else {
          parent_profiles <- find_recombinant_parents(prof_vec, all_profiles_named, threshold = recombination_threshold)
          if (!is.null(parent_profiles)) {
            # Store parent profile strings now; labels are resolved after the loop.
            parent_str <- paste(parent_profiles, collapse = "||")
            htype  <- "recombinant"; prefix <- "R"
          } else {
            htype  <- "mutant"; prefix <- "M"
          }
        }
        known_profiles[[prof_str]] <- prof_vec
        known_types[[prof_str]]    <- htype
        known_parents[[prof_str]]  <- parent_str
        if (prefix != "REF") {
          hap_labels[[prof_str]]     <- new_label(prefix)
        } else {
          hap_labels[[prof_str]] <- "REF"
        }
      }

      rows[[length(rows) + 1]] <- data.frame(
        generation   = gen,
        haplotype_id = hap_labels[[prof_str]],
        profile_str  = prof_str,
        freq         = freq,
        type         = known_types[[prof_str]],
        parents      = known_parents[[prof_str]],
        sel_coeff    = sel_coeff,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- bind_rows(rows)

  # Resolve parent profile strings to haplotype labels now that all labels are assigned.
  result$parents <- vapply(result$parents, function(ps) {
    if (is.na(ps) || nchar(ps) == 0) return(NA_character_)

    profs  <- strsplit(ps, "\\|\\|")[[1]]
    labels <- vapply(profs, function(p) {
      lbl <- hap_labels[[p]]
      if (is.null(lbl)) NA_character_ else lbl
    }, character(1L))
    paste(labels[!is.na(labels)], collapse = ",")
  }, character(1L))

  result
}

# ── Classify whole-genome haplotypes ─────────────────────────────────────────
message("Assigning per-element mutation signatures...")

#group_df <- subset(df, element_id == 8 & population_id == 0) # TESTING
element_sig_df_rds_file <- paste0(outpref, "_element_sig_df.rds")
if (!file.exists(element_sig_df_rds_file)) { 
  element_sig_df <- df %>%
    mutate(element_id_tmp = element_id) %>%
    group_by(element_id, feature_type, population_id) %>%
    group_modify(~ assign_element_sigs(.x)) %>%
    mutate(element_id_tmp = NULL) %>%
    ungroup()
  
  saveRDS(element_sig_df, element_sig_df_rds_file)
} else {
  element_sig_df <- readRDS(element_sig_df_rds_file)
}

message("Classifying whole-genome haplotypes per population...")

#pop_df <- subset(element_sig_df, population_id == 0) # TESTING
hap_data_rds_file <- paste0(outpref, "_hap_data.rds")
if (!file.exists(hap_data_rds_file)) { 
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
  
  saveRDS(hap_data, hap_data_rds_file)
} else {
  hap_data <- readRDS(hap_data_rds_file)
}

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
    p + facet_grid(
      ~ population_id,
      labeller = as_labeller(function(x) paste0("population_id: ", x)),
      ncol     = 1
    )
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

