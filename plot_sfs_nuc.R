library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggsci)

# Usage:
#   Rscript plot_sfs_nuc.R [gff_dir] [output_prefix]
#
# For each (population, generation) group, reads all matching GFF + FASTA files
# (pop_<pop>_gen_<gen>_genome_<id>.gff / .fasta), extracts the nucleotide
# sequence at each annotated feature's coordinates, then computes a per-site
# allele frequency for every variable position.  Sites are grouped by the
# feature_type of the element they belong to, and the site frequency spectrum
# (SFS) is plotted as a histogram of derived allele frequencies.
#
# GFF col-1 naming:  contig_N  (1-based)
# FASTA header suffix: _contig{N-1}  (0-based)
# Both files are expected to reside in the same directory.

args       <- commandArgs(trailingOnly = TRUE)
args       <- args[!grepl("^--", args)]
gff_dir    <- if (length(args) >= 1) args[1] else "."
out_prefix <- if (length(args) >= 2) args[2] else "sfs_nuc"

gff_dir <- "/Users/samhorsfield/Library/CloudStorage/OneDrive-Personal/Work/Postdoc_Unine/Analysis/PansimNuc/parameter_sweep/exon_mu_high_recombination_high_selection_pos_high_neg_high_proppos_equal"
out_prefix <- "/Users/samhorsfield/Library/CloudStorage/OneDrive-Personal/Work/Postdoc_Unine/Analysis/PansimNuc/parameter_sweep/exon_mu_high_recombination_high_selection_pos_high_neg_high_proppos_equal/sfs_nuc"

# ── Attribute parser ──────────────────────────────────────────────────────────
parse_attrs <- function(attr_str) {
  pairs <- strsplit(attr_str, ";", fixed = TRUE)[[1L]]
  kv    <- strsplit(pairs, "=", fixed = TRUE)
  keys  <- vapply(kv, `[[`, character(1L), 1L)
  vals  <- vapply(kv, function(x) if (length(x) >= 2L) x[[2L]] else NA_character_,
                  character(1L))
  setNames(vals, keys)
}

# ── GFF reader ────────────────────────────────────────────────────────────────
# Returns a data.frame with one row per feature.
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
      stringsAsFactors = FALSE
    )
  })
  bind_rows(Filter(Negate(is.null), rows))
}

# ── FASTA reader ──────────────────────────────────────────────────────────────
# Returns a named character vector: contig_index (as string) → sequence.
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
rev_comp <- function(seq) {
  comp <- chartr("ACGT", "TGCA", seq)
  paste(rev(strsplit(comp, "")[[1L]]), collapse = "")
}

# ── Extract feature sequence ──────────────────────────────────────────────────
# Coordinates are 1-based inclusive (GFF convention).
extract_seq <- function(fasta_seqs, contig_index, start, end, strand) {
  key <- as.character(contig_index)
  if (!key %in% names(fasta_seqs)) return(NA_character_)
  full <- fasta_seqs[[key]]
  if (end > nchar(full)) return(NA_character_)
  sub_seq <- substr(full, start, end)
  if (strand == "-") sub_seq <- rev_comp(sub_seq)
  sub_seq
}

# ── Discover GFF files ────────────────────────────────────────────────────────
gff_files <- list.files(gff_dir,
                        pattern  = "^pop_\\d+_gen_\\d+_genome_\\d+\\.gff$",
                        full.names = TRUE)

if (length(gff_files) == 0L) {
  stop("No GFF files matching pop_<pop>_gen_<gen>_genome_<id>.gff found in: ",
       gff_dir)
}
message("Found ", length(gff_files), " GFF file(s) in: ", gff_dir)

# ── Parse filenames and build metadata table ──────────────────────────────────
file_meta <- lapply(gff_files, function(fp) {
  bn <- sub("\\.gff$", "", basename(fp))
  m  <- regmatches(bn, regexpr("^pop_(\\d+)_gen_(\\d+)_genome_(\\d+)$", bn,
                               perl = TRUE))
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

# ── Extract sequences per genome ──────────────────────────────────────────────
message("Extracting nucleotide sequences from ", nrow(file_meta), " genome(s)...")

genome_seqs <- lapply(seq_len(nrow(file_meta)), function(i) {
  row <- file_meta[i, ]
  if (!file.exists(row$fasta_path)) {
    warning("FASTA not found: ", row$fasta_path)
    return(NULL)
  }
  gff   <- read_sim_gff(row$gff_path)
  if (is.null(gff) || nrow(gff) == 0L) return(NULL)
  fasta <- read_fasta(row$fasta_path)

  gff$seq <- mapply(
    extract_seq,
    contig_index = gff$contig_index,
    start        = gff$start,
    end          = gff$end,
    strand       = gff$strand,
    MoreArgs     = list(fasta_seqs = fasta)
  )

  gff$pop_id    <- row$pop_id
  gff$gen_id    <- row$gen_id
  gff$genome_id <- row$genome_id
  gff
})

all_seqs <- bind_rows(Filter(Negate(is.null), genome_seqs))

if (nrow(all_seqs) == 0L) {
  stop("No sequences could be extracted.")
}

# ── Compute per-site allele frequencies ──────────────────────────────────────
# For each (pop, gen, element_id), collect all genome sequences.
# Sequences of the same element_id must be the same length (same feature
# coordinates across genomes).  For each position, count the most common
# allele as "ancestral" and report frequency of all others as derived.

message("Computing per-site allele frequencies...")

compute_site_freqs <- function(seqs_vec) {
  # Returns a data.frame with columns major_freq and minor_freq, one row per site.
  empty <- data.frame(major_freq = numeric(0), minor_freq = numeric(0))
  valid <- seqs_vec[!is.na(seqs_vec)]
  if (length(valid) < 2L) return(empty)

  # Fast path: skip expensive matrix build if all sequences are identical
  if (length(unique(valid)) == 1L) return(empty)

  # Split all sequences into a matrix (n_genomes x seq_len)
  mat <- do.call(rbind, strsplit(valid, ""))

  # Only keep sites without N/ambiguous bases
  is_valid_col <- apply(mat, 2L, function(col) all(col %in% c("A","C","G","T")))
  mat <- mat[, is_valid_col, drop = FALSE]
  if (ncol(mat) == 0L) return(empty)
  n_valid <- nrow(mat)

  # For each site return major (most common) and minor allele frequencies
  result <- apply(mat, 2L, function(col) {
    freq <- sort(as.numeric(table(col)) / n_valid, decreasing = TRUE)
    c(major_freq = freq[1L],
      minor_freq = if (length(freq) > 1L) freq[2L] else 0.0)
  })
  # apply returns a 2 x n_sites matrix; transpose to data.frame
  data.frame(t(result), row.names = NULL)
}

# testing
# 5 genomes, 8 bp each
seqs_vec_mut <- c(
  "AGCTACGT",   # genome 1
  "AGCTACGT",   # genome 2  (monomorphic at all sites)
  "AACTATGT",   # genome 3  (site 5: A→T)
  "TAGTACGT",   # genome 4  (site 1: A→T)
  "ACGTACGT"    # genome 5
)
compute_site_freqs(seqs_vec_mut)

seqs_vec_no_mut <- c(
  "AGCTACGT",   # genome 1
  "AGCTACGT",   # genome 2
  "AGCTACGT",   # genome 3
  "AGCTACGT",   # genome 4
  "AGCTACGT"   # genome 5
)
compute_site_freqs(seqs_vec_no_mut)

sfs_data <- all_seqs |>
  filter(!is.na(seq), !is.na(element_id)) |>
  group_by(pop_id, gen_id, element_id, feature_type) |>
  summarise(
    site_freqs = list(compute_site_freqs(seq)),
    .groups    = "drop"
  ) |>
  unnest(cols = site_freqs) |>
  filter(minor_freq > 0)  # drop monomorphic sites

if (nrow(sfs_data) == 0L) {
  stop("No variable sites found.")
}

# ── Plot ──────────────────────────────────────────────────────────────────────
n_genomes_per_group <- all_seqs |>
  group_by(pop_id, gen_id) |>
  summarise(n = n_distinct(genome_id), .groups = "drop")

n_bins <- max(n_genomes_per_group$n)

p_minor <- ggplot(sfs_data, aes(x = minor_freq, fill = feature_type)) +
  geom_histogram(
    bins     = n_bins,
    boundary = 0,
    colour   = NA,
    position = "identity",
    alpha    = 0.8
  ) +
  facet_grid(
    rows = vars(interaction(pop_id, gen_id,
                             sep = " / gen=", lex.order = TRUE)),
    cols = vars(feature_type),
    labeller = labeller(
      .rows = function(x) paste0("pop=", sub(" / gen=", "  gen=", x))
    ),
    scales = "free_y"
  ) +
  scale_x_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x     = "Minor allele frequency",
    y     = "Number of sites"
  ) +
  scale_fill_npg() +
  theme_light(base_size = 11) +
  theme(legend.position = "none",
        strip.text      = element_text(size = 9))

p_minor

stacked_sfs_data <- melt(sfs_data, measure.vars = c("major_freq", "minor_freq"))
stacked_sfs_data$variable <- as.character(stacked_sfs_data$variable)
stacked_sfs_data$variable[stacked_sfs_data$variable == "major_freq"] <- "Major allele"
stacked_sfs_data$variable[stacked_sfs_data$variable == "minor_freq"] <- "Minor allele"
p_both <- ggplot(stacked_sfs_data, aes(x = value, fill = variable)) +
  geom_histogram(
    bins     = n_bins,
    boundary = 0,
    colour   = NA,
    position = "identity",
    alpha    = 0.8
  ) +
  facet_grid(
    rows = vars(interaction(pop_id, gen_id,
                            sep = " / gen=", lex.order = TRUE)),
    cols = vars(feature_type),
    labeller = labeller(
      .rows = function(x) paste0("pop=", sub(" / gen=", "  gen=", x))
    ),
    scales = "free_y"
  ) +
  scale_x_continuous(
    limits = c(0, 1.0),
    breaks = seq(0, 1.0, by = 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x     = "Allele frequency",
    y     = "Number of sites",
    fill  = "Allele type"
  ) +
  scale_fill_npg() +
  theme_light(base_size = 11) +
  theme(strip.text = element_text(size = 9))

p_both

# ── Save ──────────────────────────────────────────────────────────────────────
n_pops  <- n_distinct(sfs_data$pop_id)
n_gens  <- n_distinct(sfs_data$gen_id)
n_types <- n_distinct(sfs_data$feature_type)

pdf_w <- max(9, 4.5 * n_types)
pdf_h <- max(6, 3 * n_pops * n_gens)

out_pdf <- paste0(out_prefix, "_minor_alleles.pdf")
ggsave(out_pdf, plot = p_minor, width = pdf_w, height = pdf_h)
message("Saved: ", out_pdf)
out_pdf <- paste0(out_prefix, "_all_alleles.pdf")
ggsave(out_pdf, plot = p_both, width = pdf_w, height = pdf_h)
message("Saved: ", out_pdf)

out_csv <- paste0(out_prefix, "_sfs.csv")
write.csv(sfs_data, out_csv, row.names = FALSE)
message("Saved: ", out_csv)
