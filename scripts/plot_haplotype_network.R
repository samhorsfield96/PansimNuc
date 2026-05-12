library(dplyr)
library(ggplot2)
library(igraph)

# Usage:
#   Rscript plot_haplotype_network.R [gff_dir] [output_prefix] [generation]
#
# Arguments:
#   gff_dir        – directory containing GFF + FASTA files (default: .)
#   output_prefix  – prefix for output files       (default: haplotype_network)
#   generation     – generation to plot; "all" plots one PDF per generation,
#                    "last" plots the final generation (default: "last")
#   recombination_threshold – proportion of mutations in a haplotype that must be present to 
#                             identify a new haplotype as being a recombinant. Default 0.9 (90%). 
#
# GFF files must match: pop_<pop>_gen_<gen>_genome_<id>.gff
# FASTA files must share the same basename with a .fasta extension.
#
# Output: one PDF per (element_id, feature_type, population_id) group, per
#         requested generation. Nodes are sized by frequency, coloured by
#         haplotype type (reference / founder / mutant / recombinant), and
#         edge weights reflect the number of mutational steps between
#         haplotypes.

args    <- commandArgs(trailingOnly = TRUE)
args    <- args[!grepl("^--", args)]
gff_dir <- if (length(args) >= 1) args[1] else "."
outpref <- if (length(args) >= 2) args[2] else "haplotype_network"
gen_arg <- if (length(args) >= 3) args[3] else "all"
recombination_threshold <- if (length(args) >= 4) as.numeric(args[4]) else 0.9

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

# ── Build haplotype network for one generation snapshot ──────────────────────
# Returns a ggplot object (or NULL if < 2 haplotypes).
# checked, all good
build_network_plot <- function(snap_df, title = "") {
  # snap_df: one row per haplotype, columns: haplotype_id, sequence, freq, type, sel_coeff
  snap_df <- snap_df[snap_df$freq > 0, , drop = FALSE]

  # Collapse multiple generations in the window to one row per haplotype,
  # averaging frequency. profile_str / type / parents / sel_coeff are
  # constant per haplotype_id so take the first value.
  snap_df <- snap_df %>%
    group_by(haplotype_id) %>%
    summarise(
      freq        = mean(freq),
      profile_str = first(profile_str),
      type        = first(type),
      parents     = first(parents),
      sel_coeff   = mean(sel_coeff, na.rm = TRUE),
      .groups     = "drop"
    )

  if (nrow(snap_df) < 2) {
    message("  Fewer than 2 haplotypes — skipping network for: ", title)
    return(NULL)
  }
  
  # convert REF to empty string
  snap_df$profile_str[snap_df$profile_str == "NA"] <- ""
  n <- nrow(snap_df)

  # ── Build binary mutation-presence vectors ────────────────────────────────
  # Union of all mutation tokens across haplotypes defines the vector dimensions
  all_mutations <- unique(unlist(lapply(snap_df$profile_str, parse_sig)))
  mut_vectors   <- lapply(snap_df$profile_str, function(ps) {
    muts <- parse_sig(ps)
    as.integer(all_mutations %in% muts)
  })

  # ── Build pairwise distance matrix (Hamming steps) ───────────────────────
  dist_mat <- matrix(0L, n, n, dimnames = list(snap_df$haplotype_id, snap_df$haplotype_id))
  for (i in seq_len(n - 1)) {
    for (j in seq(i + 1, n)) {
      d <- sum(mut_vectors[[i]] != mut_vectors[[j]])
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

  has_sel <- any(!is.na(node_df$sel_coeff) & is.finite(node_df$sel_coeff))
  node_df$fill_var <- if (has_sel) node_df$sel_coeff else node_df$type

  # ── Recombinant parent edges (dashed) ─────────────────────────────────────
  recomb_nodes <- node_df[node_df$type == "recombinant" &
                           !is.na(node_df$parents) &
                           nchar(node_df$parents) > 0, , drop = FALSE]
  recomb_edge_df <- do.call(rbind, lapply(seq_len(nrow(recomb_nodes)), function(k) {
    parent_ids <- strsplit(recomb_nodes$parents[k], ",")[[1]]
    parent_ids <- parent_ids[parent_ids %in% node_df$haplotype_id]
    if (length(parent_ids) == 0) return(NULL)
    rx <- recomb_nodes$x[k]; ry <- recomb_nodes$y[k]
    do.call(rbind, lapply(parent_ids, function(pid) {
      pm <- node_df[node_df$haplotype_id == pid, ]
      data.frame(x = rx, y = ry, xend = pm$x, yend = pm$y,
                 stringsAsFactors = FALSE)
    }))
  }))

  p <- ggplot() +
    # Edges
    geom_segment(
      data = edge_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      colour    = "grey50",
      linewidth = 0.8
    ) +
    # Dashed recombinant–parent links
    {
      if (!is.null(recomb_edge_df) && nrow(recomb_edge_df) > 0) {
        geom_segment(
          data      = recomb_edge_df,
          aes(x = x, y = y, xend = xend, yend = yend),
          colour    = "#00A087",
          linewidth = 0.7,
          linetype  = "dashed"
        )
      } else NULL
    } +
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
      aes(x = x, y = y, fill = fill_var, colour = type, size = freq),
      shape  = 21,
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
    scale_colour_manual(values = use_colours, name = "Haplotype type") +
    { if (has_sel)
        scale_fill_continuous(name = "Log selection coefficient")
      else
        scale_fill_manual(values = use_colours, name = "Haplotype type", guide = "none")
    } +
    scale_size_continuous(
      name   = "Frequency",
      range  = c(3, 14)
    ) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 10)
    )

  p
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
  
  saveRDS(hap_data, hap_data_rds_file)
} else {
  hap_data <- readRDS(hap_data_rds_file)
}

# ── Resolve requested generations ────────────────────────────────────────────
all_gens <- sort(unique(hap_data$generation))

if (tolower(gen_arg) == "all") {
  # One window per consecutive pair of generations
  if (length(all_gens) < 2) {
    gen_windows <- list(all_gens)
  } else {
    gen_windows <- lapply(seq(2L, length(all_gens)), function(i) all_gens[c(i - 1L, i)])
  }
} else if (tolower(gen_arg) == "last") {
  gen_windows <- list(max(all_gens))
} else {
  g <- suppressWarnings(as.integer(gen_arg))
  if (is.na(g)) stop("generation argument must be an integer, 'all', or 'last'.")
  if (!g %in% all_gens) stop("Generation ", g, " not found in data.")
  gen_windows <- list(g)
}

# ── Generate and save plots ───────────────────────────────────────────────────
populations <- sort(unique(hap_data$population_id))

for (target_gens in gen_windows) {
  window_plots <- list()

  for (pid in populations) {
    snap <- hap_data %>%
      filter(population_id == pid, generation %in% target_gens)

    title_str <- sprintf("pop %s | gens %s", pid, paste(target_gens, collapse = ","))
    p <- build_network_plot(snap, title = title_str)
    if (!is.null(p)) {
      window_plots[[length(window_plots) + 1]] <- p
    }
  }

  if (length(window_plots) == 0) {
    message("  No plottable groups for gens: ", paste(target_gens, collapse = ","))
    next
  }

  gen_label <- paste(target_gens, collapse = "-")
  out_pdf   <- sprintf("%s_gen%s_hapnet.pdf", outpref, gen_label)
  pdf(out_pdf, width = 8, height = 7)
  for (p in window_plots) print(p)
  dev.off()
  message("  Saved: ", out_pdf)
}

message("Done.")
