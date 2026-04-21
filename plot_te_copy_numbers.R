library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)

# Usage:
#   Rscript plot_te_copy_numbers.R [output_dir] [output_prefix]
# output_dir defaults to current directory

args       <- commandArgs(trailingOnly = TRUE)
args       <- args[!grepl("^--", args)]
output_dir <- if (length(args) >= 1) args[1] else "."
outpref    <- if (length(args) >= 2) args[2] else "te_copy_numbers"

if (!dir.exists(output_dir)) {
  stop("Output directory does not exist: ", output_dir)
}

# ── Parse GFF files ───────────────────────────────────────────────────────────

gff_files <- list.files(output_dir, pattern = "^pop_\\d+_gen_\\d+_genome_\\d+\\.gff$",
                         full.names = TRUE)

if (length(gff_files) == 0) {
  stop("No GFF files matching pop_<N>_gen_<N>_genome_<N>.gff found in: ", output_dir)
}

message("Found ", length(gff_files), " GFF files in ", output_dir)

parse_attributes <- function(attr_str) {
  pairs <- strsplit(attr_str, ";")[[1]]
  kv <- strsplit(pairs, "=")
  vals <- sapply(kv, function(x) if (length(x) == 2) x[2] else NA_character_)
  names(vals) <- sapply(kv, `[[`, 1)
  vals
}

parse_gff <- function(path) {
  # Extract pop/gen/genome from filename
  bn  <- basename(path)
  m   <- regmatches(bn, regexpr("pop_(\\d+)_gen_(\\d+)_genome_(\\d+)", bn))
  parts <- strsplit(m, "_")[[1]]
  pop_id    <- as.integer(parts[2])
  gen_id    <- as.integer(parts[4])
  genome_id <- as.integer(parts[6])

  lines <- readLines(path, warn = FALSE)
  # Keep only data lines (not comments)
  data_lines <- lines[!startsWith(lines, "#") & nchar(trimws(lines)) > 0]

  if (length(data_lines) == 0) return(NULL)

  rows <- lapply(data_lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 9) return(NULL)

    feature_type <- fields[3]
    if (!feature_type %in% c("TE-CUT", "TE-COPY")) return(NULL)

    attrs <- parse_attributes(fields[9])
    element_id <- attrs["element_id"]
    log_genome_selection_coefficient = attrs["log_genome_selection_coefficient"]
    log_element_selection_coefficient = attrs["log_element_selection_coefficient"]
    multiplier = attrs["multiplier"]

    data.frame(
      pop_id     = pop_id,
      generation = gen_id,
      genome_id  = genome_id,
      feature_type = feature_type,
      element_id = element_id,
      log_genome_selection_coefficient = as.numeric(log_genome_selection_coefficient),
      log_element_selection_coefficient = as.numeric(log_element_selection_coefficient),
      multiplier = as.numeric(multiplier),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

message("Parsing GFF files...")
all_data <- do.call(rbind, Filter(Negate(is.null), lapply(gff_files, parse_gff)))

if (is.null(all_data) || nrow(all_data) == 0) {
  stop("No TE-CUT or TE-COPY features found across all GFF files.")
}

all_data$element_id <- as.integer(all_data$element_id)

# ── Copy number per genome: count distinct element_ids per genome ─────────────
# Each row is one element occurrence in one genome; element_id identifies the
# TE family/copy. Count occurrences (copy number) of each element_id per genome.

copy_counts <- all_data %>%
  group_by(pop_id, generation, genome_id, feature_type, element_id) %>%
  summarise(copies = n(), 
            log_genome_selection_coefficient = mean(log_genome_selection_coefficient),
            log_element_selection_coefficient = mean(log_element_selection_coefficient),
            multiplier = mean(multiplier),
            .groups = "drop")

# ── Distribution of copy numbers across genomes, by generation & population ───

copy_dist <- copy_counts %>%
  group_by(pop_id, generation, feature_type, copies) %>%
  summarise(frequency = n(), .groups = "drop")

# Also compute mean copy number per element across genomes
mean_copies <- copy_counts %>%
  group_by(pop_id, generation, feature_type, element_id) %>%
  summarise(mean_copies = mean(copies),
            sd_copies   = sd(copies),
            n_genomes   = n(),
            .groups = "drop")

# Total TE load per genome
total_load <- all_data %>%
  group_by(pop_id, generation, genome_id, feature_type) %>%
  summarise(total_copies = n(), .groups = "drop")

total_load_summary <- total_load %>%
  group_by(pop_id, generation, feature_type) %>%
  summarise(mean_load = mean(total_copies),
            sd_load   = sd(total_copies),
            median_load = median(total_copies),
            .groups = "drop")

# ── Write tables ──────────────────────────────────────────────────────────────

write.csv(copy_counts,        file.path(paste0(outpref, "_per_genome.csv")),   row.names = FALSE)
write.csv(copy_dist,          file.path(paste0(outpref, "_distribution.csv")), row.names = FALSE)
write.csv(mean_copies,        file.path(paste0(outpref, "_mean_per_element.csv")), row.names = FALSE)
write.csv(total_load_summary, file.path(paste0(outpref, "_total_load.csv")),   row.names = FALSE)

message("Tables written.")

# ── Plots ─────────────────────────────────────────────────────────────────────

generations <- sort(unique(all_data$generation))
populations <- sort(unique(all_data$pop_id))

# 1. Total TE load over generations, faceted by population

TE_types <- c("TE-COPY", "TE-CUT")
for (TE_type in TE_types)
{
  mean_copies_subset = subset(mean_copies, feature_type == TE_type)
  if (nrow(mean_copies_subset) == 0)
  {
    message(paste0("No ", TE_type, " present"))
    next
  }
  
  # mean copy number
  binwidth = (max(mean_copies_subset$mean_copies) - min(mean_copies_subset$mean_copies)) / 30
  if (binwidth == 0) {
    binwidth <- 1
  }
  p_dist <- ggplot(mean_copies_subset, 
                   aes(x = mean_copies, fill = feature_type, group = feature_type)) +
    #geom_histogram(aes(y = after_stat(density * nrow(mean_copies_subset))), bins = 30) +
    geom_histogram(binwidth = binwidth) +
    geom_density(aes(y = after_stat(density * (nrow(mean_copies_subset) * binwidth))), alpha = 0.25) +
    scale_fill_npg() +
    facet_grid(. ~ pop_id, labeller = label_both) +
    labs(
      x = "Mean copy number", y = "Count", fill = "TE type") +
    theme_light() +
    theme(legend.position = "none")
  
  p_dist
  ggsave(file.path(paste0(outpref, "_", TE_type, "_TE_mean_copy_dist.pdf")),
         p_dist, width = 10, height = 6)
  
  # SD of copy number
  binwidth = (max(mean_copies_subset$sd_copies) - min(mean_copies_subset$sd_copies)) / 30
  if (binwidth == 0) {
    binwidth <- 1
  }
  p_dist <- ggplot(mean_copies_subset, 
                   aes(x = sd_copies, fill = feature_type, group = feature_type)) +
    #geom_histogram(aes(y = after_stat(density * nrow(mean_copies_subset))), bins = 30) +
    geom_histogram(binwidth = binwidth) +
    geom_density(aes(y = after_stat(density * (nrow(mean_copies_subset) * binwidth))), alpha = 0.25) +
    scale_fill_npg() +
    facet_grid(. ~ pop_id, labeller = label_both) +
    labs(
      x = "SD copy number", y = "Count", fill = "TE type") +
    theme_light() +
    theme(legend.position = "none")
  
  p_dist
  ggsave(file.path(paste0(outpref, "_", TE_type, "_TE_SD_copy_dist.pdf")),
         p_dist, width = 10, height = 6)
}

