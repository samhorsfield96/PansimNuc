library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)

# Usage:
#   Rscript plot_allele_frequencies.R [tracking.csv] [output.pdf]
# Defaults to tracking.csv in the current directory and allele_freq_plot

args          <- commandArgs(trailingOnly = TRUE)
# Filter out R's own flags (e.g. --no-save, --no-restore) that leak through
args          <- args[!grepl("^--", args)]
infile <- if (length(args) >= 1) args[1] else "selection_samples.csv"
outpref <- if (length(args) >= 2) args[2] else "DFE_plot"

if (!file.exists(infile)) {
  stop("Cannot find tracking file: ", infile)
}

message("Reading ", infile)
df <- read.csv(infile, stringsAsFactors = FALSE)
long_df <- stack(df)
colnames(long_df) <- c("Coefficient", "Annotation")
max_val <- max(df)
min_val <- min(df)

p_DFE <- ggplot(long_df, aes(x = Coefficient, y= after_stat(ndensity), fill = Annotation)) +
  facet_wrap(. ~ Annotation, scales = "free_x") +
  geom_density() + 
  labs(
    x        = "Selection Coefficient",
    y        = "Density",
  ) +
  geom_vline(xintercept = 0.0) +
  #scale_x_continuous(limits = c(min_val * 1.5, max_val * 1.5)) +
  scale_fill_npg() +
  theme_light(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    axis.text.y  = element_text(size = 7),
    strip.text   = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 0.5)
  ) 

p_DFE
ggsave(paste0(outpref, ".pdf"), plot=p_DFE, width=12, height=6)
