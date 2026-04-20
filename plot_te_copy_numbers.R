library(dplyr)
library(ggplot2)
library(tidyr)

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

    data.frame(
      pop_id     = pop_id,
      generation = gen_id,
      genome_id  = genome_id,
      feature_type = feature_type,
      element_id = element_id,
      stringsAsFactors = FALSE
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
  summarise(copies = n(), .groups = "drop")

# ── Distribution of copy numbers across genomes, by generation & population ───

copy_dist <- copy_counts %>%
  group_by(pop_id, generation, feature_type, copies) %>%
  summarise(n_genomes = n(), .groups = "drop")

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

write.csv(copy_counts,        file.path(output_dir, paste0(outpref, "_per_genome.csv")),   row.names = FALSE)
write.csv(copy_dist,          file.path(output_dir, paste0(outpref, "_distribution.csv")), row.names = FALSE)
write.csv(mean_copies,        file.path(output_dir, paste0(outpref, "_mean_per_element.csv")), row.names = FALSE)
write.csv(total_load_summary, file.path(output_dir, paste0(outpref, "_total_load.csv")),   row.names = FALSE)

message("Tables written.")

# ── Plots ─────────────────────────────────────────────────────────────────────

generations <- sort(unique(all_data$generation))
populations <- sort(unique(all_data$pop_id))

# 1. Total TE load over generations, faceted by population
p_load <- ggplot(total_load_summary,
                 aes(x = factor(generation), y = mean_load,
                     colour = feature_type, group = feature_type)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean_load - sd_load,
                    ymax = mean_load + sd_load), width = 0.2) +
  facet_wrap(~pop_id, labeller = label_both) +
  labs(title = "Mean total TE copy number per genome over generations",
       x = "Generation", y = "Mean copy number", colour = "TE type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, paste0(outpref, "_total_load.pdf")),
       p_load, width = 10, height = 6)

# 2. Copy number distribution (histogram) faceted by generation and feature type
# For readability, one plot per population
for (pop in populations) {
  df_pop <- copy_dist %>% filter(pop_id == pop)
  if (nrow(df_pop) == 0) next

  p_dist <- ggplot(df_pop,
                   aes(x = copies, y = n_genomes, fill = feature_type)) +
    geom_col(position = "dodge") +
    facet_wrap(~generation, labeller = label_both, scales = "free_y") +
    labs(title = paste0("TE copy-number distribution — population ", pop),
         x = "Copies per genome", y = "Number of genomes", fill = "TE type") +
    theme_bw()

  ggsave(file.path(output_dir,
                   paste0(outpref, "_dist_pop", pop, ".pdf")),
         p_dist, width = 12, height = 8)
}

# 3. Violin / boxplot of per-genome total TE load by generation, one per population
for (pop in populations) {
  df_pop <- total_load %>% filter(pop_id == pop)
  if (nrow(df_pop) == 0) next

  p_violin <- ggplot(df_pop,
                     aes(x = factor(generation), y = total_copies,
                         fill = feature_type)) +
    geom_violin(position = position_dodge(0.8), scale = "width", alpha = 0.6) +
    geom_boxplot(position = position_dodge(0.8), width = 0.15, outlier.size = 0.5) +
    labs(title = paste0("TE copy-number distribution — population ", pop),
         x = "Generation", y = "Total copies per genome", fill = "TE type") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(output_dir,
                   paste0(outpref, "_violin_pop", pop, ".pdf")),
         p_violin, width = 10, height = 6)
}

message("Done. Output written to: ", output_dir)
