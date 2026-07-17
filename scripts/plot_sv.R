#!/usr/bin/env Rscript
# sv_plot.R — Progressive Mauve-style SV plot for PansimNuc
#
# Each genome is drawn as a single horizontal track. Contigs within each
# genome are concatenated in numeric order with a small gap between them,
# mirroring the Progressive Mauve display style. Homologous elements
# (shared element_id) are connected by synteny ribbons across genomes.
#
# Usage:
#   Rscript sv_plot.R <root.gff> <sim.directory> [options]
#
# Options:
#   --out FILE       output file (default: sv_plot.pdf)
#   --width N        plot width in inches (default: 16)
#   --height N       plot height in inches (default: n_genomes x 2.5)
#   --types T,...    comma-separated feature types to display (default: all)
#   --link-types T,. comma-separated feature types to draw links for (default: all)
#   --no-links       suppress synteny ribbons
#   --gap N          bp gap inserted between contigs (default: 500)

suppressPackageStartupMessages({
  library(gggenomes)
  library(ggnewscale)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

# ── helpers ───────────────────────────────────────────────────────────────────

parse_attrs <- function(attr_string) {
  pairs <- strsplit(attr_string, ";", fixed = TRUE)[[1L]]
  keys  <- sub("=.*$",    "", pairs)
  vals  <- sub("^[^=]+=", "", pairs)
  setNames(vals, keys)
}

read_pansimnuc_gff <- function(path, bin_label) {
  lines <- read_lines(path)
  lines <- lines[nchar(lines) > 0L & !startsWith(lines, "#")]
  if (length(lines) == 0L) {
    warning("No records found in: ", path)
    return(NULL)
  }
  rows <- lapply(lines, function(line) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(f) < 9L) return(NULL)
    a <- parse_attrs(f[9L])
    data.frame(
      bin_id        = bin_label,
      contig_name   = f[1L],
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

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  cat(
    "Usage: Rscript sv_plot.R <root.gff> <sim0.gff> [sim1.gff ...]\n",
    "  [--out sv_plot.pdf] [--width 16] [--height auto]\n",
    "  [--types exon,intron,intergenic,TE-CUT,TE-COPY]\n",
    "  [--link-types exon,intron] [--no-links] [--gap 500]\n"
  )
  quit(status = 1L)
}

take_flag <- function(flag, args, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(list(val = default, args = args))
  val  <- args[i + 1L]
  args <- args[-c(i, i + 1L)]
  list(val = val, args = args)
}
take_switch <- function(flag, args) {
  i <- match(flag, args)
  if (is.na(i)) return(list(val = FALSE, args = args))
  list(val = TRUE, args = args[-i])
}

r <- take_flag("--out",        args, "sv_plot.pdf"); out_file   <- r$val; args <- r$args
r <- take_flag("--width",      args, "16");          p_width    <- as.numeric(r$val); args <- r$args
r <- take_flag("--height",     args, NULL);          p_height   <- if (is.null(r$val)) NULL else as.numeric(r$val); args <- r$args
r <- take_flag("--types",      args, NULL);          keep_types <- if (is.null(r$val)) NULL else strsplit(r$val, ",")[[1L]]; args <- r$args
r <- take_flag("--link-types", args, NULL);          link_types <- if (is.null(r$val)) NULL else strsplit(r$val, ",")[[1L]]; args <- r$args
r <- take_flag("--gap",        args, "500");         contig_gap <- as.integer(r$val); args <- r$args
r <- take_switch("--no-links", args);                no_links   <- r$val; args <- r$args

root_path <- args[1L]
sim_directory <- args[-1L]

# ── read data ─────────────────────────────────────────────────────────────────

# root_path <- "/Users/samhorsfield/Software/PansimNuc/parameter_sweep/baseline/root_out.gff"
# sim_directory <- "/Users/samhorsfield/Software/PansimNuc/parameter_sweep/baseline"
# contig_gap <- 500
# no_links <- FALSE
# link_types <- c("TE-COPY", "TE-CUT")
# keep_types <- c("exon", "intron", "intergenic", "TE-COPY", "TE-CUT")
# p_width <- 16
# p_height <- 16
# out_file <- "/Users/samhorsfield/Software/PansimNuc/parameter_sweep/baseline/sv_plot.pdf"

message("Reading root GFF: ", root_path)
all_feats <- read_pansimnuc_gff(root_path, "root")

sim_paths <- list.files(sim_directory, pattern = "*.gff", full.names = TRUE)
sim_paths <- sim_paths[sim_paths != root_path]

for (i in seq_along(sim_paths)) {
  filename <- basename(sim_paths[i])
  label <- as.character(as.numeric(gsub("([0-9]+).*$", "\\1", filename)))
  message("Reading ", label, ": ", sim_paths[i])
  block <- read_pansimnuc_gff(sim_paths[i], label)
  if (!is.null(block)) all_feats <- bind_rows(all_feats, block)
}

if (is.null(all_feats) || nrow(all_feats) == 0L) {
  stop("No features loaded. Check that the GFF files are valid PansimNuc output.")
}

if (!is.null(keep_types)) {
  all_feats <- filter(all_feats, feature_type %in% keep_types)
}

# ── linearize contigs (Progressive Mauve style) ──────────────────────────────
# Within each genome, sort contigs by their numeric suffix (contig_1 < contig_2
# < ...), concatenate them end-to-end with `contig_gap` bp between each pair,
# and shift all feature coordinates into that linearized space. Each genome
# becomes a single seq_id so gggenomes draws it as one track.

# Per-contig length and numeric sort key
contig_info <- all_feats |>
  mutate(contig_num = as.integer(str_extract(contig_name, "[0-9]+"))) |>
  group_by(bin_id, contig_name, contig_num) |>
  summarise(contig_len = max(end), .groups = "drop") |>
  arrange(bin_id, contig_num)

# Cumulative left-edge offset within each genome
contig_info <- contig_info |>
  group_by(bin_id) |>
  mutate(
    offset      = cumsum(lag(contig_len + contig_gap, default = 0L)),
    contig_rank = row_number()          # alternating shade index
  ) |>
  ungroup()

# Contig block rectangles for background shading
contig_blocks <- contig_info |>
  transmute(
    bin_id,
    seq_id       = bin_id,
    start        = offset + 1L,
    end          = offset + contig_len,
    contig_name,
    contig_shade = factor((contig_rank - 1L) %% 2L)  # "0" / "1" alternating
  )

# Shift all feature coordinates into linearized space
all_feats <- all_feats |>
  mutate(contig_num = as.integer(str_extract(contig_name, "[0-9]+"))) |>
  left_join(select(contig_info, bin_id, contig_name, offset),
            by = c("bin_id", "contig_name")) |>
  mutate(
    seq_id = bin_id,
    start  = start + offset,
    end    = end   + offset
  )

# ── seqs table ────────────────────────────────────────────────────────────────
# One row per genome, ordered root first then numerically by genome index.

seqs <- all_feats |>
  group_by(bin_id) |>
  summarise(seq_id = first(bin_id), length = max(end), .groups = "drop") |>
  mutate(
    sort_key = if_else(
      bin_id == "root", -1L,
      suppressWarnings(as.integer(str_extract(bin_id, "[0-9]+")))
    )
  ) |>
  arrange(sort_key) |>
  select(-sort_key)

# ── links table ───────────────────────────────────────────────────────────────
# Connect every root element to matching elements in simulated genomes via
# element_id. Crossed links = translocations; fan-out = duplications;
# absent link = deletion.

if (!no_links) {
  ordered_genomes <- seqs$seq_id  # root, 0, 1, 2, ... in figure order

  links_list <- vector("list", length(ordered_genomes) - 1L)
  for (i in seq_len(length(ordered_genomes) - 1L)) {
    upper_anchors <- all_feats |>
      filter(seq_id == ordered_genomes[i]) |>
      select(seq_id = seq_id, start = start, end = end,
             strand1 = strand, element_id, feature_type)

    if (!is.null(link_types)) {
      upper_anchors <- filter(upper_anchors, feature_type %in% link_types)
    }

    lower_anchors <- all_feats |>
      filter(seq_id == ordered_genomes[i + 1L]) |>
      select(seq_id2 = seq_id, start2 = start, end2 = end,
             strand2 = strand, element_id, feature_type2 = feature_type)

    links_list[[i]] <- inner_join(upper_anchors, lower_anchors,
                                  by = "element_id",
                                  relationship = "many-to-many") |>
      select(element_id, seq_id, start, end, strand1,
             seq_id2, start2, end2, strand2, feature_type, feature_type2)
  }

  links <- bind_rows(links_list)
  links$strand <- ifelse(links$strand1 == links$strand2, "+", "-")
  if (nrow(links) == 0L) links <- NULL
} else {
  links <- NULL
}

# ── colour palette ────────────────────────────────────────────────────────────

feature_colors <- c(
  exon       = "#4DAF4A",
  intron     = "#984EA3",
  intergenic = "#999999",
  "TE-CUT"   = "#E41A1C",
  "TE-COPY"  = "#FF7F00"
)
extra_types <- setdiff(unique(all_feats$feature_type), names(feature_colors))
if (length(extra_types) > 0L) {
  feature_colors <- c(
    feature_colors,
    setNames(hue_pal()(length(extra_types)), extra_types)
  )
}

# ── plot dimensions ───────────────────────────────────────────────────────────

n_bins <- nrow(seqs)
if (is.null(p_height)) p_height <- max(4.0, n_bins * 2.5)

# ── build gggenomes plot ──────────────────────────────────────────────────────

if (!no_links && !is.null(links) && nrow(links) > 0L) {
  p <- gggenomes(
    seqs  = seqs,
    genes = all_feats,
    links = links,
    feats = list(contigs = contig_blocks)
  )
} else {
  p <- gggenomes(
    seqs  = seqs,
    genes = all_feats,
    feats = list(contigs = contig_blocks)
  )
}

# Layer 1 — alternating contig background shading
p <- p +
  geom_feat(
    data        = feats("contigs"),
    aes(fill    = contig_shade),
    alpha       = 0.12,
    linewidth   = NA,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("0" = "#AAAAAA", "1" = "#555555"),
    guide  = "none"
  ) +
  new_scale_fill()

# Layer 2 — sequence backbone and genome labels
p <- p +
  geom_seq() +
  geom_seq_label()

# Layer 3 — synteny ribbons coloured by feature type
if (!no_links && !is.null(links) && nrow(links) > 0L) {
  p <- p +
    geom_link(
      aes(fill = feature_type),
      alpha  = 0.28,
      colour = NA
    ) + scale_fill_manual(
      values = feature_colors, na.value = "grey60",
      name   = "Feature type"
    ) +
    new_scale_fill()
}

# Layer 4 — gene features
p <- p +
  geom_gene(aes(fill = feature_type)) +
  scale_fill_manual(
    values = feature_colors, na.value = "grey60",
    name   = "Feature type"
  ) +
  theme_gggenomes_clean()

# ── save ──────────────────────────────────────────────────────────────────────

message(sprintf("Writing %s  (%.0f x %.0f in)", out_file, p_width, p_height))
ggsave(out_file, p, width = p_width, height = p_height)
message("Done.")
