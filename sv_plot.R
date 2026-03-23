#!/usr/bin/env Rscript
# sv_plot.R — Visualise PansimNuc structural variation with gggenomes
#
# Usage:
#   Rscript sv_plot.R <root.gff> <sim0.gff> [sim1.gff ...] \
#                     [--out plot.pdf] [--width 16] [--height auto] \
#                     [--types exon,intron,intergenic,TE-CUT,TE-COPY] \
#                     [--link-types exon,intron] \
#                     [--no-links]
#
# Arguments:
#   root.gff      Root / reference GFF (first positional argument)
#   sim*.gff      One or more simulation GFF files
#   --out         Output file path (default: sv_plot.pdf)
#   --width       Plot width in inches (default: 16)
#   --height      Plot height in inches (default: auto based on genome count)
#   --types       Comma-separated feature types to display (default: all)
#   --link-types  Comma-separated feature types to draw links for (default: all)
#   --no-links    Suppress synteny links entirely

suppressPackageStartupMessages({
  library(gggenomes)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# ── helpers ───────────────────────────────────────────────────────────────────

parse_attrs <- function(attr_string) {
  pairs <- strsplit(attr_string, ";", fixed = TRUE)[[1L]]
  keys  <- sub("=.*$",    "", pairs)
  vals  <- sub("^[^=]+=", "", pairs)
  setNames(vals, keys)
}

read_pansimnuc_gff <- function(path, bin_label) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nchar(lines) > 0L & !startsWith(lines, "#")]
  if (length(lines) == 0L) {
    warning("No records found in: ", path)
    return(NULL)
  }

  line <- lines[1]
  
  rows <- lapply(lines, function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(f) < 9L) return(NULL)
    a <- parse_attrs(f[9L])
    data.frame(
      bin_id        = bin_label,
      contig_id     = f[1L],
      start         = as.integer(f[4L]),
      end           = as.integer(f[5L]),
      strand        = f[7L],
      feature_type  = a[["feature_type"]],
      element_id    = suppressWarnings(as.integer(a[["element_id"]])),
      feature_id    = suppressWarnings(as.integer(a[["feature_id"]])),
      multiplier    = suppressWarnings(as.numeric(a[["multiplier"]])),
      log_sel_coeff = suppressWarnings(as.numeric(a[["log_element_selection_coefficient"]])),
      log_genome_sel_coeff = suppressWarnings(as.numeric(a[["log_genome_selection_coefficient"]])),
      stringsAsFactors = FALSE
    )
  })
  

  bind_rows(Filter(Negate(is.null), rows))
}

# ── CLI argument parsing ──────────────────────────────────────────────────────

# args <- commandArgs(trailingOnly = TRUE)
# 
# if (length(args) < 2L) {
#   cat(
#     "Usage: Rscript sv_plot.R <root.gff> <sim0.gff> [sim1.gff ...]\n",
#     "  [--out plot.pdf] [--width 16] [--height auto]\n",
#     "  [--types exon,intron,intergenic,TE-CUT,TE-COPY]\n",
#     "  [--link-types exon,intron] [--no-links]\n"
#   )
#   quit(status = 1L)
# }

# # helper to consume a named flag and its value from args vector
# take_flag <- function(flag, args, default = NULL) {
#   i <- match(flag, args)
#   if (is.na(i)) return(list(val = default, args = args))
#   val  <- args[i + 1L]
#   args <- args[-c(i, i + 1L)]
#   list(val = val, args = args)
# }
# take_switch <- function(flag, args) {
#   i <- match(flag, args)
#   if (is.na(i)) return(list(val = FALSE, args = args))
#   list(val = TRUE, args = args[-i])
# }
# 
# r <- take_flag("--out",        args, "sv_plot.pdf"); out_file   <- r$val; args <- r$args
# r <- take_flag("--width",      args, "16");          p_width    <- as.numeric(r$val); args <- r$args
# r <- take_flag("--height",     args, NULL);          p_height   <- if (is.null(r$val)) NULL else as.numeric(r$val); args <- r$args
# r <- take_flag("--types",      args, NULL);          keep_types <- if (is.null(r$val)) NULL else strsplit(r$val, ",")[[1L]]; args <- r$args
# r <- take_flag("--link-types", args, NULL);          link_types <- if (is.null(r$val)) NULL else strsplit(r$val, ",")[[1L]]; args <- r$args
# r <- take_switch("--no-links", args);                no_links   <- r$val; args <- r$args

# root_path <- args[1L]
# sim_paths <- args[-1L]

# ── read data ─────────────────────────────────────────────────────────────────

root_path <- "/Users/samhorsfield/Software/PansimNuc/output/root_test_output.gff"

message("Reading root GFF: ", root_path)
all_feats <- read_pansimnuc_gff(root_path, "root")

sim_directory <- "/Users/samhorsfield/Software/PansimNuc/output"
sim_paths <- list.files(sim.directory, pattern = "*.gff", full.names = TRUE)
sim_paths <- sim.files[sim.files != root_path]

i <- 1
for (i in seq_along(sim.paths)) {
  filename <- basename(sim_paths[i])
  label <- as.character(as.numeric(gsub("([0-9]+).*$", "\\1", filename)))
  message("Reading ", label, ": ", sim_paths[i])
  block <- read_pansimnuc_gff(sim_paths[i], label)
  if (!is.null(block)) all_feats <- bind_rows(all_feats, block)
}

if (is.null(all_feats) || nrow(all_feats) == 0L) {
  stop("No features loaded. Check that the GFF files are valid PansimNuc output.")
}

# Filter to requested feature types
# if (!is.null(keep_types)) {
#   all_feats <- filter(all_feats, feature_type %in% keep_types)
# }

# ── seqs table ────────────────────────────────────────────────────────────────
# One row per contig per genome: seq_id, bin_id, length

seqs <- all_feats |>
  group_by(bin_id, seq_id) |>
  summarise(length = max(end), .groups = "drop")

# ── links table ───────────────────────────────────────────────────────────────
# For each element_id, connect every root occurrence to every simulated
# occurrence. This exposes duplications (multiple copies in a derived genome)
# as well as translocations (changed position) and inversions (changed strand).

if (!no_links) {
  root_anchors <- all_feats |>
    filter(bin_id == "root") |>
    select(seq_id1 = seq_id, start1 = start, end1 = end, element_id, feature_type)

  sim_anchors <- all_feats |>
    filter(bin_id != "root") |>
    select(seq_id2 = seq_id, start2 = start, end2 = end, element_id)

  # Restrict links to requested feature types (useful to reduce clutter)
  if (!is.null(link_types)) {
    root_anchors <- filter(root_anchors, feature_type %in% link_types)
  }

  links <- inner_join(root_anchors, sim_anchors,
                      by = "element_id",
                      relationship = "many-to-many")
} else {
  links <- NULL
}

# ── feature colour palette ────────────────────────────────────────────────────

feature_colors <- c(
  exon       = "#4DAF4A",
  intron     = "#984EA3",
  intergenic = "#999999",
  "TE-CUT"   = "#E41A1C",
  "TE-COPY"  = "#FF7F00"
)

# Expand palette for any feature types not listed above
extra_types <- setdiff(unique(all_feats$feature_type), names(feature_colors))
if (length(extra_types) > 0L) {
  extra_cols <- setNames(
    scales::hue_pal()(length(extra_types)),
    extra_types
  )
  feature_colors <- c(feature_colors, extra_cols)
}

# ── plot dimensions ───────────────────────────────────────────────────────────

n_bins <- length(unique(all_feats$bin_id))
if (is.null(p_height)) p_height <- max(4, n_bins * 2.5)

# ── build gggenomes plot ──────────────────────────────────────────────────────

if (!no_links && !is.null(links) && nrow(links) > 0L) {
  p <- gggenomes(seqs = seqs, genes = all_feats, links = links)
} else {
  p <- gggenomes(seqs = seqs, genes = all_feats)
}

p <- p +
  geom_seq() +
  geom_seq_label()

if (!no_links && !is.null(links) && nrow(links) > 0L) {
  p <- p + geom_link(aes(fill = feature_type), alpha = 0.35)
}

p <- p +
  geom_gene(aes(fill = feature_type), size = 3) +
  scale_fill_manual(values = feature_colors, na.value = "grey60") +
  labs(
    title    = "PansimNuc structural variation",
    subtitle = sprintf(
      "%d simulated genome(s) vs root | linked by element_id",
      length(sim_paths)
    ),
    fill = "Feature type"
  ) +
  theme_gggenomes_clean()

# ── save ──────────────────────────────────────────────────────────────────────

message(sprintf("Writing %s (%.0f x %.0f in)", out_file, p_width, p_height))
ggsave(out_file, p, width = p_width, height = p_height)
message("Done.")
