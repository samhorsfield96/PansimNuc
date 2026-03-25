library(dplyr)
library(ggplot2)
library(stringr)

read_pansimnuc_gff <- function(path, bin_label) {
  lines <- readLines(path, warn = FALSE)
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

sim_dir   <- "/Users/samhorsfield/Software/PansimNuc/parameter_sweep/baseline"
root_path <- file.path(sim_dir, "root_out.gff")

# Read root to get ancestral element types
root_feats <- read_pansimnuc_gff(root_path, "root")

# Read all simulated genomes
sim_paths <- setdiff(
  list.files(sim_dir, pattern = "\\.gff$", full.names = TRUE),
  root_path
)
n_genomes <- length(sim_paths)

sim_feats <- lapply(seq_along(sim_paths), function(i) {
  label <- as.character(as.numeric(gsub("([0-9]+).*$", "\\1", basename(sim_paths[i]))))
  read_pansimnuc_gff(sim_paths[i], label)
}) |> bind_rows()

# Annotate element type from the root (many element_ids won't be in root
# if they are novel insertions — label those "novel")
root_types <- root_feats |>
  select(element_id, feature_type) |>
  distinct()

# Count how many genomes each element_id appears in
afs <- sim_feats |>
  select(bin_id, element_id, feature_type) |>
  distinct() |>                               # one row per genome per element
  group_by(element_id, feature_type) |>
  summarise(count = n_distinct(bin_id), .groups = "drop") |>
  mutate(frequency = count / n_genomes)

# Plot
ggplot(afs, aes(x = frequency, fill = feature_type)) +
  geom_histogram(bins = n_genomes, boundary = 0, colour = NA) +
  facet_wrap(~feature_type, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x     = "Allele frequency (fraction of genomes)",
    y     = "Number of elements",
    title = "Allele frequency spectrum by feature type"
  ) +
  theme_bw() +
  theme(legend.position = "none")
