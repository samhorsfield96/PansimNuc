library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)

# Usage:
#   Rscript ld_analysis.R [gff_dir] [output_prefix] [flank_bp] [generation]
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
outpref    <- if (length(args) >= 2) args[2] else "ld_analysis"
flank_bp   <- if (length(args) >= 3) as.integer(args[3]) else 100000L
gen_arg    <- if (length(args) >= 4) args[4] else "last"


gff_dir <- "/Users/samhorsfield/Library/CloudStorage/OneDrive-Personal/Work/Postdoc_Unine/Analysis/PansimNuc_results/baseline_uniform_selection_no_demography_recomb_all_gens"
outpref <- "/Users/samhorsfield/Library/CloudStorage/OneDrive-Personal/Work/Postdoc_Unine/Analysis/PansimNuc_results/baseline_uniform_selection_no_demography_recomb_all_gens_ld_analysis"

message(sprintf("Parameters:  flank_bp=%d  generation=%s", flank_bp, gen_arg))

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
extract_window <- function(fasta_seqs, contig_index, seq_start, seq_end, strand) {
  key  <- as.character(contig_index)
  if (!key %in% names(fasta_seqs)) return(NA_character_)
  
  full <- fasta_seqs[[key]]
  
  # take full sequences if NA
  if (is.na(seq_start)) {
    start <- 1L
  } else {
    start <- max(1L, seq_start)
  }
  if (is.na(seq_end)) {
    end <- nchar(full)
  } else {
    end   <- min(nchar(full), seq_end)
  }
  
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

# ── Filter file_meta to requested generation(s) ───────────────────────────────
all_gens <- sort(unique(file_meta$gen_id))
if (tolower(gen_arg) == "all") {
  # keep all generations (no filtering)
} else if (tolower(gen_arg) == "last") {
  file_meta <- file_meta[file_meta$gen_id == max(all_gens), ]
} else {
  g <- suppressWarnings(as.integer(gen_arg))
  if (is.na(g)) stop("generation argument must be an integer, 'all', or 'last'.")
  if (!g %in% all_gens) stop("Generation ", g, " not found in data.")
  file_meta <- file_meta[file_meta$gen_id == g, ]
}
message(sprintf("Retaining %d file(s) after generation filter.", nrow(file_meta)))

# ── Read all GFF features (with selection coefficients) ───────────────────────
message("Reading GFF files to extract selection coefficients...")
gff_df_rds_file <- paste0(outpref, "_gff_df.rds")
if (!file.exists(gff_df_rds_file)) {
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
  saveRDS(all_gff, gff_df_rds_file)
} else {
  all_gff <- readRDS(gff_df_rds_file)
}

# ── For each top element, extract per-genome sequences for element + flanks ───
message("Extracting sequences for focal elements and flanking regions...")

# Build site matrix from pre-extracted locus sequences across genomes.
# seqs_per_element: list of window sequences, each N-padded so that the element
#                   always starts at offset flank_bp (1-based); NAs for missing.
# el_rows:          all GFF rows for this element across genomes (for el length).
# flank_bp:         flanking region size used during extraction.
# Positions are element-relative (0 = first base of element).
# checked, all good
build_site_matrix <- function(seqs_per_element, el_rows, flank_bp) {
  valid <- Filter(function(s) !is.null(s) && length(s) == 1L && !is.na(s),
                  seqs_per_element)
  if (length(valid) < 2L) return(NULL)

  # Trim all sequences to minimum length
  lens  <- nchar(valid)
  min_l <- min(lens)
  valid <- lapply(valid, substr, 1L, min_l)

  # Build nucleotide matrix (n_genomes × min_l)
  nuc_mat <- do.call(rbind, strsplit(unlist(valid), ""))

  # Filter to ACGT-only columns
  valid_col <- apply(nuc_mat, 2L, function(col) all(col %in% c("A","C","G","T")))
  nuc_mat   <- nuc_mat[, valid_col, drop = FALSE]
  if (ncol(nuc_mat) == 0L) return(NULL)

  # Keep only variable sites
  is_var  <- apply(nuc_mat, 2L, function(col) length(unique(col)) > 1L)
  nuc_mat <- nuc_mat[, is_var, drop = FALSE]
  if (ncol(nuc_mat) == 0L) return(NULL)

  # Element-relative positions: window index j (1-based) → j - 1 - flank_bp
  all_pos   <- seq_len(min_l) - 1L - flank_bp
  positions <- all_pos[valid_col][is_var]

  el_lens    <- el_rows$end - el_rows$start + 1L
  median_len <- as.integer(median(el_lens))
  in_element <- positions >= 0L & positions < median_len

  # Convert to biallelic 0/1 (major allele = 0)
  site_mat <- apply(nuc_mat, 2L, function(col) {
    tbl   <- table(col)
    major <- names(which.max(tbl))
    as.integer(col != major)
  })
  if (!is.matrix(site_mat)) site_mat <- matrix(site_mat, ncol = 1L)

  list(
    site_matrix    = site_mat,
    site_positions = positions,
    in_element     = in_element,
    locus_start    = -flank_bp,
    locus_end      = median_len + flank_bp - 1L,
    el_start       = 0L,
    el_end         = median_len - 1L,
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

# ── Pre-load all FASTA sequences into memory ─────────────────────────────────
message("Pre-loading FASTA sequences into memory...")
fasta_cache <- list()
for (.i in seq_len(nrow(file_meta))) {
  .row <- file_meta[.i, ]
  if (!file.exists(.row$fasta_path)) next
  .key <- paste(.row$pop_id, .row$gen_id, .row$genome_id, sep = "_")
  fasta_cache[[.key]] <- read_fasta(.row$fasta_path)
}
message(sprintf("Loaded %d FASTA file(s) into memory.", length(fasta_cache)))

# ── Main loop ─────────────────────────────────────────────────────────────────
ld_results <- list()

all_ld_rds_file <- paste0(outpref, "_all_ld.rds")
if (!file.exists(all_ld_rds_file)) {
  generations <- sort(unique(all_gff$gen_id))
  populations <- sort(unique(all_gff$pop_id))

  for (pop_id_val in populations) {
    message(sprintf("Processing Pop. %d", pop_id_val))
    for (gen_id_val in generations) {
      gff_pg <- all_gff |> filter(pop_id == pop_id_val, gen_id == gen_id_val)
      if (nrow(gff_pg) == 0L) next
      message(sprintf("  Processing Gen. %d", gen_id_val))

      genomes     <- file_meta |> filter(pop_id == pop_id_val, gen_id == gen_id_val)
      element_ids <- sort(unique(gff_pg$element_id))

      pb <- txtProgressBar(min = 1, max = length(element_ids), style = 3)
      for (ei in seq_along(element_ids)) {
        eid <- element_ids[ei]
        setTxtProgressBar(pb, ei)

        # All GFF rows for this element across all genomes in this pop/gen
        el_rows <- gff_pg |> filter(element_id == eid)

        # Extract locus window from each genome using that genome's coordinates.
        # Upstream truncation at contig start is compensated by N-padding so
        # the element always starts at offset flank_bp within the window.
        seqs_per_element <- lapply(seq_len(nrow(genomes)), function(i) {
          genome_row    <- genomes[i, ]
          cache_key     <- paste(genome_row$pop_id, genome_row$gen_id,
                                 genome_row$genome_id, sep = "_")
          fasta         <- fasta_cache[[cache_key]]
          if (is.null(fasta)) return(NA)
          el_g <- el_rows |> filter(genome_id == genome_row$genome_id)
          if (nrow(el_g) == 0L) return(NA)
          el_g          <- el_g[1L, ]
          desired_start <- el_g$start - flank_bp
          actual_start  <- max(1L, desired_start)
          seq           <- extract_window(fasta, el_g$contig_index,
                                          actual_start, el_g$end + flank_bp, "+")
          if (is.na(seq)) return(NA)
          n_pad <- actual_start - desired_start   # > 0 only when near contig start
          if (n_pad > 0L) seq <- paste0(strrep("N", n_pad), seq)
          seq
        })

        el_row <- el_rows[1L, ]   # representative row for metadata

        sm <- build_site_matrix(seqs_per_element, el_rows, flank_bp)
        if (is.null(sm)) next

        ld <- compute_r2(sm$site_matrix, sm$in_element, sm$site_positions)
        if (is.null(ld) || nrow(ld) == 0L) next

        ld$pop_id       <- pop_id_val
        ld$gen_id       <- gen_id_val
        ld$element_id   <- eid
        ld$feature_type <- el_row$feature_type
        ld$mean_log_sel <- mean(el_rows$log_sel_coeff, na.rm = TRUE)
        ld$contig_index <- el_row$contig_index
        ld$el_start     <- sm$el_start
        ld$el_end       <- sm$el_end
        ld$locus_start  <- sm$locus_start
        ld$locus_end    <- sm$locus_end
        ld$n_genomes    <- sm$n_genomes

        ld_results[[length(ld_results) + 1L]] <- ld
      }
      close(pb)
    }
  }

  all_ld <- bind_rows(Filter(Negate(is.null), ld_results))
  saveRDS(all_ld, all_ld_rds_file)
} else {
  all_ld <- readRDS(all_ld_rds_file)
}

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

locus_x_limits <- ld_summary |>
  distinct(facet_label, locus_start, locus_end) |>
  tidyr::pivot_longer(c(locus_start, locus_end),
                      names_to = "bound", values_to = "partner_pos") |>
  mutate(mean_r2 = 0)

p_ld <- ggplot(ld_summary,
               aes(x = partner_pos, y = mean_r2,
                   colour = in_partner_el)) +
  geom_blank(data = locus_x_limits,
             aes(x = partner_pos, y = mean_r2),
             inherit.aes = FALSE) +
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
heat_xy_limits <- all_ld |>
  distinct(pop_id, gen_id, element_id, locus_start, locus_end) |>
  tidyr::pivot_longer(c(locus_start, locus_end),
                      names_to = "bound", values_to = "pos") |>
  mutate(focal_pos = pos, partner_pos = pos, r2 = NA_real_)

p_heat <- ggplot(all_ld |> filter(!is.na(r2)),
                 aes(x = focal_pos, y = partner_pos, fill = r2)) +
  geom_blank(data = heat_xy_limits,
             aes(x = focal_pos, y = partner_pos),
             inherit.aes = FALSE) +
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
pdf_h_decay <- max(4, 3.5 * n_elements)
pdf_h_heat  <- max(6, 4  * ceiling(n_elements / 2))
pdf_w       <- 10

out_decay <- paste0(outpref, "_ld_decay.pdf")
ggsave(out_decay, plot = p_ld, width = pdf_w, height = pdf_h_decay, limitsize = FALSE)
message("Saved: ", out_decay)

out_heat <- paste0(outpref, "_ld_heatmap.pdf")
ggsave(out_heat, plot = p_heat, width = pdf_w, height = pdf_h_heat, limitsize = FALSE)
message("Saved: ", out_heat)

out_csv_ld  <- paste0(outpref, "_ld_pairwise.csv")
write.csv(all_ld, out_csv_ld, row.names = FALSE)
message("Saved: ", out_csv_ld)



