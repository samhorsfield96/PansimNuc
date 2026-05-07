library(dplyr)
library(ggplot2)
library(tidyr)
library(igraph)
library(ape)
library(pegas)

# Usage:
#   Rscript plot_haplotype_network.R [tracking.csv] [output_prefix] [generation]
#
# Arguments:
#   tracking.csv   – path to the tracking CSV (default: tracking.csv)
#   output_prefix  – prefix for output files       (default: haplotype_network)
#   generation     – generation to plot; "all" plots one PDF per generation,
#                    "last" plots the final generation (default: "last")
#
# Output: one PDF per (element_id, feature_type, population_id) group, per
#         requested generation. Nodes are sized by frequency, coloured by
#         haplotype type (reference / founder / mutant / recombinant), and
#         edge weights reflect the number of mutational steps between
#         haplotypes.

args          <- commandArgs(trailingOnly = TRUE)
args          <- args[!grepl("^--", args)]
tracking_file <- if (length(args) >= 1) args[1] else "tracking.csv"
outpref       <- if (length(args) >= 2) args[2] else "haplotype_network"
gen_arg       <- if (length(args) >= 3) args[3] else "last"

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

get_mutation_sig <- function(seq, reference) {
  chars <- toupper(strsplit(seq, "")[[1]])
  len   <- min(length(chars), length(reference))
  idx   <- which(chars[seq_len(len)] != reference[seq_len(len)])
  if (length(idx) == 0) return("")
  paste(paste0(idx, ":", chars[idx]), collapse = ";")
}

parse_sig <- function(s) if (nchar(s) == 0) character(0) else strsplit(s, ";")[[1]]

is_recombinant <- function(h_muts, known_muts_list) {
  if (length(h_muts) == 0 || length(known_muts_list) < 2) return(FALSE)
  h_len <- length(h_muts)
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
      if (!any(!a %in% b) || !any(!b %in% a)) next
      ab <- union(a, b)
      if (length(ab) == h_len && setequal(ab, h_muts)) return(TRUE)
    }
  }
  FALSE
}

hap_sel_coeff <- function(coeff_strs) {
  coeff_strs <- coeff_strs[!is.na(coeff_strs) & nchar(coeff_strs) > 0]
  if (length(coeff_strs) == 0) return(NA_real_)
  per_genome <- sapply(coeff_strs, function(s) {
    vals <- suppressWarnings(as.numeric(strsplit(s, ";")[[1]]))
    sum(vals, na.rm = TRUE)
  })
  mean(per_genome, na.rm = TRUE)
}

# Hamming distance between two equal-length sequences (character vectors or strings).
hamming <- function(a, b) {
  ca <- strsplit(toupper(a), "")[[1]]
  cb <- strsplit(toupper(b), "")[[1]]
  len <- min(length(ca), length(cb))
  sum(ca[seq_len(len)] != cb[seq_len(len)])
}

# ── Classify haplotypes for a single group ────────────────────────────────────
classify_haplotypes <- function(group_df) {
  generations   <- sort(unique(group_df$generation))
  first_gen     <- generations[1]
  has_sel_coeff <- "log_selection_coefficients" %in% colnames(group_df)

  reference <- build_reference(group_df$sequence[group_df$generation == first_gen])

  known_haps  <- list()
  hap_labels  <- list()
  known_muts  <- list()
  counter     <- 0L
  new_label   <- function(prefix) { counter <<- counter + 1L; paste0(prefix, counter) }

  rows <- list()
  for (gen in generations) {
    gen_mask <- group_df$generation == gen & nchar(group_df$sequence) > 0
    sel_all  <- if (has_sel_coeff) group_df$log_selection_coefficients[gen_mask] else NULL
    seqs_all <- group_df$sequence[gen_mask]
    n_total  <- length(seqs_all)
    if (n_total == 0) next

    seq_tbl    <- table(seqs_all)
    sel_by_seq <- if (has_sel_coeff) split(sel_all, seqs_all) else NULL

    for (seq in names(seq_tbl)) {
      sig      <- get_mutation_sig(seq, reference)
      sig_muts <- parse_sig(sig)
      freq     <- seq_tbl[[seq]] / n_total

      sel_coeff <- if (has_sel_coeff) hap_sel_coeff(sel_by_seq[[seq]]) else NA_real_

      if (!sig %in% names(known_haps)) {
        if (gen == first_gen) {
          htype  <- "founder";     prefix <- "F"
        } else if (is_recombinant(sig_muts, known_muts)) {
          htype  <- "recombinant"; prefix <- "R"
        } else {
          htype  <- "mutant";      prefix <- "M"
        }
        known_haps[[sig]] <- htype
        hap_labels[[sig]] <- new_label(prefix)
        known_muts[[sig]] <- sig_muts
      }

      if (sig == "") {
        haplotype_id <- "Ref"
        htype        <- "reference"
      } else {
        haplotype_id <- hap_labels[[sig]]
        htype        <- known_haps[[sig]]
      }

      rows[[length(rows) + 1]] <- data.frame(
        generation   = gen,
        haplotype_id = haplotype_id,
        mutation_sig = sig,
        sequence     = seq,
        freq         = freq,
        type         = htype,
        sel_coeff    = sel_coeff,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(rows)
}

# ── Build haplotype network for one generation snapshot ──────────────────────
# Returns a ggplot object (or NULL if < 2 haplotypes).
build_network_plot <- function(snap_df, title = "") {
  # snap_df: one row per haplotype, columns: haplotype_id, sequence, freq, type, sel_coeff
  snap_df <- snap_df[nchar(snap_df$sequence) > 0 & snap_df$freq > 0, , drop = FALSE]
  if (nrow(snap_df) < 2) {
    message("  Fewer than 2 haplotypes — skipping network for: ", title)
    return(NULL)
  }

  n <- nrow(snap_df)

  # ── Build pairwise distance matrix (Hamming steps) ───────────────────────
  dist_mat <- matrix(0L, n, n, dimnames = list(snap_df$haplotype_id, snap_df$haplotype_id))
  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      d <- hamming(snap_df$sequence[i], snap_df$sequence[j])
      dist_mat[i, j] <- d
      dist_mat[j, i] <- d
    }
  }

  # ── Minimum spanning tree via igraph ─────────────────────────────────────
  g_full <- graph_from_adjacency_matrix(
    dist_mat,
    mode    = "undirected",
    weighted = TRUE,
    diag    = FALSE
  )
  g_mst <- mst(g_full, weights = E(g_full)$weight)

  # ── Layout ────────────────────────────────────────────────────────────────
  set.seed(42)
  layout_mat <- layout_with_fr(g_mst)
  colnames(layout_mat) <- c("x", "y")

  node_df <- snap_df
  node_df$x <- layout_mat[, "x"]
  node_df$y <- layout_mat[, "y"]

  # ── Edge data frame ───────────────────────────────────────────────────────
  el       <- as_edgelist(g_mst)
  e_weight <- E(g_mst)$weight
  edge_df  <- data.frame(
    from   = el[, 1],
    to     = el[, 2],
    weight = e_weight,
    stringsAsFactors = FALSE
  )
  # Mid-point for edge label (mutation step count)
  from_xy <- node_df[match(edge_df$from, node_df$haplotype_id), c("x", "y")]
  to_xy   <- node_df[match(edge_df$to,   node_df$haplotype_id), c("x", "y")]
  edge_df$mx <- (from_xy$x + to_xy$x) / 2
  edge_df$my <- (from_xy$y + to_xy$y) / 2
  edge_df$x    <- from_xy$x
  edge_df$y    <- from_xy$y
  edge_df$xend <- to_xy$x
  edge_df$yend <- to_xy$y

  # ── Colour palette ────────────────────────────────────────────────────────
  type_colours <- c(
    reference   = "#3C5488FF",
    founder     = "#4DBBD5",
    mutant      = "#E64B35",
    recombinant = "#00A087"
  )
  # Ensure all types present in palette
  present_types <- unique(node_df$type)
  use_colours   <- type_colours[names(type_colours) %in% present_types]

  # ── Plot ──────────────────────────────────────────────────────────────────
  # Node radius proportional to frequency (sqrt for area)
  node_df$size <- 4 + 10 * sqrt(node_df$freq)

  p <- ggplot() +
    # Edges
    geom_segment(
      data = edge_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      colour    = "grey50",
      linewidth = 0.8
    ) +
    # Edge weight tick marks – small cross-ticks every mutational step along each edge
    {
      # Build tick-mark data: for each step along an edge, one small perpendicular segment
      ticks <- do.call(rbind, lapply(seq_len(nrow(edge_df)), function(k) {
        steps <- edge_df$weight[k]
        if (steps <= 1L) return(NULL)
        dx  <- edge_df$xend[k] - edge_df$x[k]
        dy  <- edge_df$yend[k] - edge_df$y[k]
        len <- sqrt(dx^2 + dy^2)
        if (len == 0) return(NULL)
        # perpendicular unit vector
        px  <- -dy / len * 0.15
        py  <-  dx / len * 0.15
        # positions at 1/(steps) intervals, excluding endpoints
        ts  <- seq_len(steps - 1) / steps
        data.frame(
          x    = edge_df$x[k] + ts * dx - px,
          y    = edge_df$y[k] + ts * dy - py,
          xend = edge_df$x[k] + ts * dx + px,
          yend = edge_df$y[k] + ts * dy + py
        )
      }))
      if (!is.null(ticks) && nrow(ticks) > 0) {
        geom_segment(data = ticks,
                     aes(x = x, y = y, xend = xend, yend = yend),
                     colour = "grey40", linewidth = 0.5)
      } else {
        NULL
      }
    } +
    # Edge step labels (only when > 1 step)
    geom_text(
      data = edge_df[edge_df$weight > 1, , drop = FALSE],
      aes(x = mx, y = my, label = weight),
      size   = 2.8,
      colour = "grey30",
      vjust  = -0.4
    ) +
    # Nodes
    geom_point(
      data = node_df,
      aes(x = x, y = y, fill = type, size = freq),
      shape  = 21,
      colour = "white",
      stroke = 0.6,
      alpha  = 0.9
    ) +
    # Node labels
    geom_text(
      data  = node_df,
      aes(x = x, y = y, label = haplotype_id),
      size  = 2.8,
      vjust = -1.2,
      fontface = "bold"
    ) +
    scale_fill_manual(values = use_colours, name = "Haplotype type") +
    scale_size_continuous(
      name   = "Frequency",
      range  = c(3, 14),
      limits = c(0, 1)
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 10)
    )

  p
}

# ── Classify all groups ───────────────────────────────────────────────────────
message("Classifying haplotypes across all groups...")

hap_data <- df %>%
  group_by(element_id, feature_type, population_id) %>%
  group_modify(~ classify_haplotypes(.x)) %>%
  ungroup()

# ── Resolve requested generations ────────────────────────────────────────────
all_gens <- sort(unique(hap_data$generation))

if (tolower(gen_arg) == "all") {
  target_gens <- all_gens
} else if (tolower(gen_arg) == "last") {
  target_gens <- max(all_gens)
} else {
  g <- suppressWarnings(as.integer(gen_arg))
  if (is.na(g)) stop("generation argument must be an integer, 'all', or 'last'.")
  if (!g %in% all_gens) stop("Generation ", g, " not found in data.")
  target_gens <- g
}

# ── Generate plots ────────────────────────────────────────────────────────────
groups <- hap_data %>%
  distinct(element_id, feature_type, population_id)

for (gen in target_gens) {
  message(sprintf("Plotting generation %d ...", gen))

  gen_plots <- list()
  gen_titles <- character(0)

  for (row_i in seq_len(nrow(groups))) {
    eid   <- groups$element_id[row_i]
    ftype <- groups$feature_type[row_i]
    pid   <- groups$population_id[row_i]

    snap <- hap_data %>%
      filter(
        element_id    == eid,
        feature_type  == ftype,
        population_id == pid,
        generation    == gen
      ) %>%
      # Keep one row per unique haplotype (take first occurrence, freq is per-gen)
      distinct(haplotype_id, .keep_all = TRUE)

    title_str <- sprintf("element %s | %s | pop %s | gen %d", eid, ftype, pid, gen)
    p <- build_network_plot(snap, title = title_str)
    if (!is.null(p)) {
      gen_plots[[length(gen_plots) + 1]] <- p
      gen_titles <- c(gen_titles, title_str)
    }
  }

  if (length(gen_plots) == 0) {
    message("  No plottable groups for generation ", gen)
    next
  }

  # Write one PDF per generation with one page per group
  out_pdf <- sprintf("%s_gen%04d.pdf", outpref, gen)
  pdf(out_pdf, width = 8, height = 7)
  for (p in gen_plots) print(p)
  dev.off()
  message("  Saved: ", out_pdf)
}

message("Done.")
