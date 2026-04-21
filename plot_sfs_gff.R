library(dplyr)
library(ggplot2)

# Usage:
#   Rscript plot_sfs_gff.R [gff_dir] [output_prefix]
#
# Reads all GFF files matching pop_<pop>_gen_<gen>_genome_<id>.gff in gff_dir.
# For each (population_id, generation) group, computes the minor allele
# frequency (MAF) of each element_id across all sampled genomes, then plots
# the site frequency spectrum (SFS) as a histogram of MAF values.

args       <- commandArgs(trailingOnly = TRUE)
args       <- args[!grepl("^--", args)]
gff_dir    <- if (length(args) >= 1) args[1] else "."
out_prefix <- if (length(args) >= 2) args[2] else "sfs_plot"

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
read_sim_gff <- function(path, pop_id, gen_id, genome_id) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nchar(lines) > 0L & !startsWith(lines, "#")]
  if (length(lines) == 0L) {
    warning("No records in: ", path)
    return(NULL)
  }
  rows <- lapply(lines, function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(f) < 9L) return(NULL)
    a <- parse_attrs(f[9L])
    data.frame(
      population_id = pop_id,
      generation    = gen_id,
      genome_id     = genome_id,
      element_id    = suppressWarnings(as.integer(a[["element_id"]])),
      feature_type  = a[["feature_type"]],
      stringsAsFactors = FALSE
    )
  })
  bind_rows(Filter(Negate(is.null), rows))
}

# ── Discover files ────────────────────────────────────────────────────────────
gff_files <- list.files(gff_dir, pattern = "^pop_\\d+_gen_\\d+_genome_\\d+\\.gff$",
                        full.names = TRUE)

if (length(gff_files) == 0L) {
  stop("No files matching pop_<pop>_gen_<gen>_genome_<id>.gff found in: ", gff_dir)
}

message("Found ", length(gff_files), " GFF file(s) in: ", gff_dir)

# ── Parse filename metadata and read ─────────────────────────────────────────
all_feats <- lapply(gff_files, function(fp) {
  bn  <- sub("\\.gff$", "", basename(fp))
  m   <- regmatches(bn, regexpr("^pop_(\\d+)_gen_(\\d+)_genome_(\\d+)$", bn,
                                perl = TRUE))
  if (length(m) == 0L) {
    warning("Skipping unrecognised filename: ", basename(fp))
    return(NULL)
  }
  parts  <- as.integer(strsplit(sub("^pop_", "", m), "_gen_|_genome_")[[1L]])
  pop_id <- parts[1L]
  gen_id <- parts[2L]
  gid    <- parts[3L]
  read_sim_gff(fp, pop_id, gen_id, gid)
}) |> bind_rows()

if (nrow(all_feats) == 0L) {
  stop("No features were read from any GFF file.")
}

# ── Compute SFS ───────────────────────────────────────────────────────────────
# For each group (population x generation), determine the number of sampled
# genomes, then count how many carry each element_id.  Minor allele frequency
# is min(count / n, 1 - count / n).

genome_counts <- all_feats |>
  distinct(population_id, generation, genome_id) |>
  group_by(population_id, generation) |>
  summarise(n_genomes = n_distinct(genome_id), .groups = "drop")

# One row per (population, generation, element_id) – deduplicate within a genome
element_presence <- all_feats |>
  select(population_id, generation, genome_id, element_id, feature_type) |>
  distinct()

sfs <- element_presence |>
  group_by(population_id, generation, element_id, feature_type) |>
  summarise(count = n_distinct(genome_id), .groups = "drop") |>
  left_join(genome_counts, by = c("population_id", "generation")) |>
  mutate(
    frequency = count / n_genomes,
    maf       = pmin(frequency, 1.0 - frequency)
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
n_bins   <- max(genome_counts$n_genomes) / 2L  # bin per MAF step of 1/(2n)
pop_gens <- sort(unique(paste0("pop=", sfs$population_id,
                               "  gen=", sfs$generation)))

p <- ggplot(sfs, aes(x = maf, fill = feature_type)) +
  geom_histogram(bins     = n_bins,
                 boundary = 0,
                 colour   = NA,
                 position = "identity",
                 alpha    = 0.8) +
  facet_grid(
    rows = vars(interaction(population_id, generation,
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
    y     = "Number of loci",
    title = "Site frequency spectrum by feature type"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        strip.text      = element_text(size = 9))

# ── Save ──────────────────────────────────────────────────────────────────────
n_pops  <- n_distinct(sfs$population_id)
n_gens  <- n_distinct(sfs$generation)
n_types <- n_distinct(sfs$feature_type)

pdf_w <- max(4, 2.5 * n_types)
pdf_h <- max(4, 2.5 * n_pops * n_gens)

out_pdf <- paste0(out_prefix, ".pdf")
ggsave(out_pdf, plot = p, width = pdf_w, height = pdf_h)
message("Saved: ", out_pdf)

out_csv <- paste0(out_prefix, "_sfs.csv")
write.csv(sfs, out_csv, row.names = FALSE)
message("Saved: ", out_csv)
