library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)

# Usage:
#   Rscript ld_analysis.R [gff_dir] [output_prefix] [top_n] [flank_bp]
#
# For each (population, generation) group, reads all matching GFF + FASTA files
# (pop_<pop>_gen_<gen>_genome_<id>.gff / .fasta), identifies the N elements with
# the highest (least negative) log_element_selection_coefficient, extracts all
# variable sites within those elements plus a flanking region, and computes
# pairwise linkage disequilibrium (r²) between the focal element variants and
# all variants in the surrounding locus.
#
# GFF attribute fields parsed:
#   element_id, feature_type, log_element_selection_coefficient
# GFF col-1: contig_N  (1-based)
# FASTA header suffix: _contig{N-1}  (0-based)

args       <- commandArgs(trailingOnly = TRUE)
args       <- args[!grepl("^--", args)]
gff_dir    <- if (length(args) >= 1) args[1] else "."
out_prefix <- if (length(args) >= 2) args[2] else "ld_analysis"
top_n      <- if (length(args) >= 3) as.integer(args[3]) else 5L
flank_bp   <- if (length(args) >= 4) as.integer(args[4]) else 10000L

gff_dir <- "/Users/samhorsfield/OneDrive/Work/Postdoc_Unine/Analysis/PansimNuc_results/testing_high_mu_low_selection"
out_prefix <- "/Users/samhorsfield/OneDrive/Work/Postdoc_Unine/Analysis/PansimNuc_results/ld_analysis"
top_n <- 2
flank_bp <- 1000L

message(sprintf("Parameters: top_n=%d  flank_bp=%d", top_n, flank_bp))

# ── Attribute parser ──────────────────────────────────────────────────────────
# checked, all good
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

# ── Extract a genomic window from a FASTA contig ──────────────────────────────
# Coordinates are 1-based inclusive (GFF convention); returns NA if out of range.
# checked, all good
extract_window <- function(fasta_seqs, contig_index, start, end, strand) {
  key  <- as.character(contig_index)
  if (!key %in% names(fasta_seqs)) return(NA_character_)
  full <- fasta_seqs[[key]]
  start <- max(1L, start)
  end   <- min(nchar(full), end)
  if (start > end) return(NA_character_)
  sub_seq <- substr(full, start, end)
  if (strand == "-") sub_seq <- rev_comp(sub_seq)
  sub_seq
}

# ── Discover GFF files ────────────────────────────────────────────────────────
gff_files <- list.files(gff_dir,
                        pattern    = "^pop_\\d+_gen_\\d+_genome_\\d+\\.gff$",
                        full.names = TRUE)

if (length(gff_files) == 0L) {
  stop("No GFF files matching pop_<pop>_gen_<gen>_genome_<id>.gff found in: ",
       gff_dir)
}
message("Found ", length(gff_files), " GFF file(s) in: ", gff_dir)

# ── Parse filenames → metadata ────────────────────────────────────────────────
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

# ── Read all GFF features (with selection coefficients) ───────────────────────
message("Reading GFF files to extract selection coefficients...")

all_gff <- lapply(seq_len(nrow(file_meta)), function(i) {
  row <- file_meta[i, ]
  gff <- read_sim_gff(row$gff_path)
  if (is.null(gff) || nrow(gff) == 0L) return(NULL)
  gff$pop_id    <- row$pop_id
  gff$gen_id    <- row$gen_id
  gff$genome_id <- row$genome_id
  gff
})
all_gff <- bind_rows(Filter(Negate(is.null), all_gff))

if (nrow(all_gff) == 0L) stop("No GFF records could be parsed.")

# ── Identify top-N elements per (pop, gen) by mean log selection coefficient ──
# Higher (less negative) log_sel_coeff → stronger positive selection.
# We rank by mean across genomes to get population-level signal.
message(sprintf("Identifying top %d elements per (pop, gen) group by selection coefficient...", top_n))

element_sel <- all_gff |>
  filter(!is.na(log_sel_coeff), !is.na(element_id)) |>
  group_by(pop_id, gen_id, element_id, feature_type,
           contig_index, start, end, strand) |>
  summarise(mean_log_sel = mean(log_sel_coeff, na.rm = TRUE), .groups = "drop")

top_elements <- element_sel |>
  group_by(pop_id, gen_id) |>
  slice_max(order_by = mean_log_sel, n = top_n, with_ties = FALSE) |>
  ungroup()

# testing
top_elements <- element_sel[element_sel$element_id == 8,]

message(sprintf("Selected %d element × (pop,gen) combinations total.", nrow(top_elements)))
print(top_elements)

# ── For each top element, extract per-genome sequences for element + flanks ───
message("Extracting sequences for focal elements and flanking regions...")

# Build site matrix for a locus (element ± flank_bp) across all genomes.
# Returns a list with:
#   $site_matrix   : integer matrix  (n_genomes × n_sites),  0/1 biallelic
#   $site_positions: integer vector  (n_sites)  genomic positions (1-based)
#   $in_element    : logical vector  (n_sites)  TRUE if inside element boundaries
# checked, all good
build_site_matrix <- function(pop_id_val, gen_id_val, el_row, file_meta, flank_bp) {
  genomes <- file_meta |>
    filter(pop_id == pop_id_val, gen_id == gen_id_val)

  contig_idx  <- el_row$contig_index
  el_start    <- el_row$start        # GFF 1-based
  el_end      <- el_row$end
  locus_start <- max(1L, el_start - flank_bp)
  locus_end   <- el_end + flank_bp

  i <- 1
  seqs_per_genome <- lapply(seq_len(nrow(genomes)), function(i) {
    row <- genomes[i, ]
    if (!file.exists(row$fasta_path)) return(NULL)
    fasta <- read_fasta(row$fasta_path)
    seq   <- extract_window(fasta, contig_idx,
                            locus_start, locus_end, "+")  # always + strand for LD window
    if (is.na(seq)) return(NULL)
    seq
  })

  valid <- Filter(Negate(is.null), seqs_per_genome)
  if (length(valid) < 2L) return(NULL)

  # Pad / trim to common length (should all be the same unless near contig edge)
  lens   <- nchar(valid)
  min_l  <- min(lens)
  valid  <- lapply(valid, substr, 1L, min_l)
  actual_end <- locus_start + min_l - 1L

  # Build nucleotide matrix  (n_genomes × min_l)
  nuc_mat <- do.call(rbind, strsplit(unlist(valid), ""))

  # Filter to ACGT-only columns
  valid_col <- apply(nuc_mat, 2L, function(col) all(col %in% c("A","C","G","T")))
  nuc_mat   <- nuc_mat[, valid_col, drop = FALSE]
  positions <- (locus_start:actual_end)[valid_col]

  if (ncol(nuc_mat) == 0L) return(NULL)

  # Keep only variable sites
  is_var <- apply(nuc_mat, 2L, function(col) length(unique(col)) > 1L)
  nuc_mat   <- nuc_mat[, is_var, drop = FALSE]
  positions <- positions[is_var]

  if (ncol(nuc_mat) == 0L) return(NULL)

  # Convert to biallelic 0/1 (major allele = 0)
  site_mat <- apply(nuc_mat, 2L, function(col) {
    tbl    <- table(col)
    major  <- names(which.max(tbl))
    as.integer(col != major)
  })

  in_element <- positions >= el_start & positions <= el_end

  list(
    site_matrix    = site_mat,
    site_positions = positions,
    in_element     = in_element,
    locus_start    = locus_start,
    locus_end      = actual_end,
    el_start       = el_start,
    el_end         = el_end,
    n_genomes      = nrow(site_mat)
  )
}


# function for R-squared
rsq <- function (x, y) cor(x, y) ^ 2


# ── Compute r² between each focal-element site and every locus site ───────────
#checked, all good
compute_r2 <- function(site_mat, in_element, positions) {
  focal_idx  <- which(in_element)
  if (length(focal_idx) == 0L || ncol(site_mat) < 2L) return(NULL)

  n <- nrow(site_mat)
  fi <- 1 
  results <- lapply(focal_idx, function(fi) {
    x <- site_mat[, fi]
    # mean of focal site data in terms of minor allele frequency
    px <- mean(x)
    if (px == 0 || px == 1) return(NULL)   # monomorphic (shouldn't happen after filter)

    r2_vec <- sapply(seq_len(ncol(site_mat)), function(j) {
      y  <- site_mat[, j]
      # mean of observed data
      py <- mean(y)
      if (py == 0 || py == 1) return(NA_real_)
      rsq(x, y)
    })

    data.frame(
      focal_pos    = positions[fi],
      partner_pos  = positions,
      r2           = r2_vec,
      in_focal_el  = in_element[fi],
      in_partner_el = in_element,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(Filter(Negate(is.null), results))
}

# ── Main loop ─────────────────────────────────────────────────────────────────
ld_results <- list()

# checked, all good
for (row_i in seq_len(nrow(top_elements))) {
  el_row <- top_elements[row_i, ]
  pg_label <- sprintf("pop=%s gen=%s element_id=%s (%s)",
                      el_row$pop_id, el_row$gen_id,
                      el_row$element_id, el_row$feature_type)
  message("Processing: ", pg_label)

  pop_id_val <- el_row$pop_id
  gen_id_val <- el_row$gen_id
  sm <- build_site_matrix(pop_id_val, gen_id_val, el_row,
                           file_meta, flank_bp)
  if (is.null(sm)) {
    message("  → skipped (insufficient data)")
    next
  }

  message(sprintf("  Locus %d–%d  |  %d variable sites  |  %d genomes",
                  sm$locus_start, sm$locus_end,
                  ncol(sm$site_matrix), sm$n_genomes))

  ld <- compute_r2(sm$site_matrix, sm$in_element, sm$site_positions)
  if (is.null(ld) || nrow(ld) == 0L) {
    message("  → no variable focal sites")
    next
  }

  ld$pop_id        <- el_row$pop_id
  ld$gen_id        <- el_row$gen_id
  ld$element_id    <- el_row$element_id
  ld$feature_type  <- el_row$feature_type
  ld$mean_log_sel  <- el_row$mean_log_sel
  ld$contig_index  <- el_row$contig_index
  ld$el_start      <- el_row$start
  ld$el_end        <- el_row$end
  ld$locus_start   <- sm$locus_start
  ld$locus_end     <- sm$locus_end
  ld$n_genomes     <- sm$n_genomes

  ld_results[[row_i]] <- ld
}

all_ld <- bind_rows(Filter(Negate(is.null), ld_results))

if (nrow(all_ld) == 0L) stop("No LD data computed.")

# ── Summary: mean r² per (pop, gen, element_id) × partner position ───────────
ld_summary <- all_ld |>
  group_by(pop_id, gen_id, element_id, feature_type,
           mean_log_sel, partner_pos, in_partner_el,
           el_start, el_end, locus_start, locus_end) |>
  summarise(mean_r2 = mean(r2, na.rm = TRUE), .groups = "drop")

# ── Plot: r² decay across locus for each top element ─────────────────────────
make_ld_label <- function(pop, gen, eid, ftype, sel) {
  sprintf("pop=%s  gen=%s  |  element %s (%s)  |  mean log-s=%.3f",
          pop, gen, eid, ftype, sel)
}

ld_summary$facet_label <- mapply(
  make_ld_label,
  ld_summary$pop_id, ld_summary$gen_id,
  ld_summary$element_id, ld_summary$feature_type,
  ld_summary$mean_log_sel
)

p_ld <- ggplot(ld_summary,
               aes(x = partner_pos, y = mean_r2,
                   colour = in_partner_el)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_smooth(data    = ld_summary |> filter(!in_partner_el),
              aes(group = 1),
              method  = "loess", span = 0.3, se = FALSE,
              colour  = "grey40", linewidth = 0.6, linetype = "dashed") +
  geom_rect(aes(xmin = el_start, xmax = el_end,
                ymin = -Inf,     ymax = Inf),
            fill = "yellow", alpha = 0.006, colour = NA,
            inherit.aes = FALSE,
            data = ld_summary |>
              distinct(facet_label, el_start, el_end)) +
  facet_wrap(~ facet_label, ncol = 1, scales = "free_x") +
  scale_colour_manual(
    values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
    labels = c("TRUE" = "Within focal element", "FALSE" = "Flanking region"),
    name   = "Variant Type"
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Genomic position (bp, 1-based GFF coords)",
    y = expression(Mean~r^2)
  ) +
  theme_light(base_size = 11) +
  theme(strip.text = element_text(size = 7.5),
        legend.position = "right")
p_ld

# ── Heatmap: pairwise r² within each locus (focal × all sites) ───────────────
p_heat <- ggplot(all_ld |> filter(!is.na(r2)),
                 aes(x = focal_pos, y = partner_pos, fill = r2)) +
  geom_tile() +
  facet_wrap(~interaction(pop_id, gen_id, element_id, sep = " / "),
             scales = "free", ncol = 2) +
  scale_fill_gradientn(
    colours = c("white", "#4DBBD5", "#E64B35"),
    limits  = c(0, 1),
    name    = expression(r^2)
  ) +
  labs(
    x     = "Focal-element site position",
    y     = "Partner site position",
  ) +
  theme_light(base_size = 10) +
  theme(strip.text = element_text(size = 7))
p_heat

# ── Save ──────────────────────────────────────────────────────────────────────
n_elements <- n_distinct(ld_summary$element_id)
pdf_h_decay <- max(4, 3.5 * nrow(top_elements))
pdf_h_heat  <- max(6, 4  * ceiling(n_elements / 2))
pdf_w       <- 10

out_decay <- paste0(out_prefix, "_ld_decay.pdf")
ggsave(out_decay, plot = p_ld, width = pdf_w, height = pdf_h_decay, limitsize = FALSE)
message("Saved: ", out_decay)

out_heat <- paste0(out_prefix, "_ld_heatmap.pdf")
ggsave(out_heat, plot = p_heat, width = pdf_w, height = pdf_h_heat, limitsize = FALSE)
message("Saved: ", out_heat)

out_csv_ld  <- paste0(out_prefix, "_ld_pairwise.csv")
write.csv(all_ld, out_csv_ld, row.names = FALSE)
message("Saved: ", out_csv_ld)

out_csv_top <- paste0(out_prefix, "_top_elements.csv")
write.csv(top_elements, out_csv_top, row.names = FALSE)
message("Saved: ", out_csv_top)
