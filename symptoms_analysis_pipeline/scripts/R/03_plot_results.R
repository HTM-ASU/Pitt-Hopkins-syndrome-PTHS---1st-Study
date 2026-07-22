#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readxl)
  library(yaml)
  library(tidyr)
  library(ggplot2)
  library(dplyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b
root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(root, "config", "config.yaml"))) {
  stop("config/config.yaml not found. Run this script from the repository root.")
}

cfg <- yaml::read_yaml(file.path(root, "config", "config.yaml"))

input_file <- file.path(root, cfg$paths$input_excel)
change_file <- file.path(root, cfg$paths$processed_dir, cfg$files$ttest_results)
out_pdf <- file.path(root, cfg$paths$figures_dir, cfg$files$figure_pdf)

participant_col <- cfg$columns$participant
group_col <- cfg$columns$group
timepoint_col <- cfg$columns$timepoint
score_cols <- unlist(cfg$columns$scores)

x_order <- unlist(cfg$plot$x_order)
x_labels <- unlist(cfg$plot$x_labels)
names(x_labels) <- x_order

palette <- c("GROUP A" = "#27b85c", "GROUP B" = "purple")

df_symptoms_KN <- readxl::read_excel(input_file, sheet = cfg$paths$input_sheet)
p_values_df <- readxl::read_excel(change_file)

needed_cols <- c(participant_col, group_col, timepoint_col, score_cols)
missing <- setdiff(needed_cols, names(df_symptoms_KN))
if (length(missing) > 0) stop(paste("Missing columns in input file:", paste(missing, collapse = ", ")))

main_df <- df_symptoms_KN %>%
  select(all_of(needed_cols)) %>%
  filter(.data[[group_col]] %in% c(cfg$groups$treated, cfg$groups$control)) %>%
  mutate(TimeGroup = paste0(.data[[timepoint_col]], "_", .data[[group_col]]))

long_df <- main_df %>%
  pivot_longer(cols = all_of(score_cols), names_to = "Feature", values_to = "Value") %>%
  mutate(
    Value = suppressWarnings(as.numeric(Value)),
    TimePoint = factor(.data[[timepoint_col]], levels = x_order),
    GroupBlock = .data[[group_col]]
  )

p_values_df <- p_values_df %>%
  mutate(
    feature = as.character(feature),
    timepoint = as.character(timepoint)
  )

# Pretty axis labels
pretty_y <- c(
  "FLACC" = "FLACC",
  "DSR Total Events" = "DSR Total Events",
  "GSRS AVG" = "GSRS AVG"
)

figures_list <- list()

for (taxon in score_cols) {
  df.1 <- long_df %>%
    filter(Feature == taxon)

  if (nrow(df.1) == 0) next

  stat_test <- p_values_df %>%
    filter(feature == taxon, timepoint == cfg$plot$show_pvalue_only_for)

  p <- ggplot(df.1, aes(x = TimePoint, y = Value, fill = .data[[group_col]])) +
    geom_boxplot(position = position_dodge(width = 0.65), outlier.shape = NA, width = 0.6) +
    geom_point(
      aes(group = GroupBlock),
      position = position_jitterdodge(
        jitter.width = 0.5,
        jitter.height = 0,
        dodge.width = 0.8,
        seed = 42
      ),
      shape = 21, size = 3, stroke = 0.5, colour = "black"
    ) +
    scale_x_discrete(limits = x_order, labels = x_labels, drop = FALSE) +
    scale_fill_manual(values = palette) +
    labs(y = pretty_y[[taxon]] %||% taxon, x = NULL) +
    theme_classic() +
    theme(
      plot.margin = margin(5, 5, 10, 5),
      axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 18, face = "bold"),
      axis.text.y = element_text(color = "black", face = "bold", size = 18),
      axis.title.y = element_text(size = 18, face = "bold"),
      legend.position = "none"
    )

  if (nrow(stat_test) > 0) {
    max_y <- max(df.1$Value, na.rm = TRUE)
    if (!is.finite(max_y)) max_y <- 0
    y_pos <- max_y * 1.10 + ifelse(max_y == 0, 1, 0)
    label <- sprintf("p=%.2f", stat_test$p_value[1])

    p <- p + geom_text(
      data = data.frame(TimePoint = cfg$plot$show_pvalue_only_for, y = y_pos, label = label),
      aes(x = TimePoint, y = y, label = label),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    )
  }

  figures_list[[taxon]] <- p
}

pdf(out_pdf, width = 7, height = 4.5)
for (plot in figures_list) print(plot)
dev.off()
cat("Wrote:", out_pdf, "\n")
