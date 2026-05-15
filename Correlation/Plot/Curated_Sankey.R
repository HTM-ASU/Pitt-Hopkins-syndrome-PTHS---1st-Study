#########################
# Final Sankey Diagram #
########################

library(readxl)
library(dplyr)
library(ggalluvial)
library(ggplot2)
library(ggnewscale)
library(stringr)
library(tidyr)

# === Utility: Clean taxonomic labels ===
clean_feature_label <- function(feature_string) {
  if (grepl("^[dkpcofgs]__", feature_string)) {
    matches <- stringr::str_extract_all(feature_string, "[dkpcofgs]__[^.;|_]+")[[1]]
    if (length(matches) == 0) return(feature_string)
    last_level <- tail(matches, 1)
    return(last_level)
  } else {
    return(feature_string)
  }
}

# === Load and prepare data ===
df_t1 <- read_excel("Top_Omics_Symptom_Correlations_T1_CLR.xlsx", sheet = "Top20")
df_t2 <- read_excel("Top_Omics_Symptom_Correlations_T2_CLR.xlsx", sheet = "Top20")

df_t1 <- df_t1 %>%
  rename(mean_delta = `mean_delta (MTT-Baseline)`,
         mean_delta_score = `mean_delta_SymptomsScore(MTT-Baseline)`) %>%
  mutate(comparison = "10 wks vs. Baseline")

df_t2 <- df_t2 %>%
  mutate(comparison = "3 months vs. Baseline")

df_combined <- bind_rows(df_t1, df_t2)

# === Clean labels ===
df_plot <- df_combined %>%
  mutate(
    weight = abs(correlation),
    direction = ifelse(correlation > 0, "Positive", "Negative"),
    feature_display = case_when(
      grepl("^\\d+\\.\\d+\\.\\d+\\.\\d+:.*", feature) ~ {
        name_part <- sub("^\\d+\\.\\d+\\.\\d+\\.\\d+:(.*)", "\\1", feature)
        #name_part <- str_replace(name_part, "\\s+[Aa]ctivity$", "")
        name_part <- str_squish(name_part)
        name_part
      },
      grepl("^[dkpcofgs]__", feature) ~ sapply(feature, clean_feature_label),
      
      grepl("^GH\\d+", feature) ~ {
        number_part <- feature
        paste("Glycoside Hydrolase Family", number_part)
      },
  
      grepl("^GT\\d+", feature) ~ {
        number_part <- feature
        paste("GlycosylTransferase  Family", number_part)
      },
      
      TRUE ~ str_squish(feature)
    ),
    #Symptom_display = str_to_title(as.character(Symptom)),
    Symptom_display = case_when(
      Symptom == "GSRS Subscores- indigestion" ~ "Indigestion",
      Symptom == "Ataxia (lack of communication between brain and body)"    ~ "Ataxia",
      Symptom == "GSRS Subscores- abdominal pain" ~ "Abdominal Pain",
      Symptom == "# of days 1 or 2"          ~ "Number of Days (DSR)",
      TRUE ~ str_to_title(as.character(Symptom))  # fallback: keep existing
    ),
    
    flow_group = case_when(
      direction == "Positive" & comparison == "10 wks vs. Baseline" ~ "Pos. (10 wks vs. Baseline)",
      direction == "Positive" & comparison == "3 months vs. Baseline" ~ "Pos. (3 months vs. Baseline)",
      direction == "Negative" & comparison == "10 wks vs. Baseline" ~ "Neg. (10 wks vs. Baseline)",
      direction == "Negative" & comparison == "3 months vs. Baseline" ~ "Neg. (3 months vs. Baseline)"
    ),
    alluvium_id = paste(feature, Symptom_display, comparison, sep = "_") # <-- key fix
  )

# === Build lodes ===
cor_lodes <- bind_rows(
  df_plot %>%
    transmute(
      axis = "Feature",
      stratum = feature_display,
      alluvium = alluvium_id,
      weight,
      flow_group,
      featureType_plot = gsub("_", " ", featureType)
    ),
  df_plot %>%
    transmute(
      axis = "Symptom",
      stratum = Symptom_display,
      alluvium = alluvium_id,
      weight,
      flow_group,
      featureType_plot = "Symptom"
    )
)

# === Color palettes ===
flow_colors <- c(
  "Pos. (3 months vs. Baseline)" = "#33a02c",
  "Pos. (10 wks vs. Baseline)" = "#1f78b4",
  "Neg. (10 wks vs. Baseline)" = "#e31a1c",
  "Neg. (3 months vs. Baseline)" = "#ff7f00"
)

full_palette <- c(
  "16S Bacteria"      = "#1f77b4",
  "Shotgun Bacteria"  = "#aec7e8",
  "ITS Fungi"         = "#ff7f0e",
  "Shotgun Fungi"     = "#ffbb78",
  "AMR"               = "#2ca02c",
  "Virulence"         = "#98df8a",
  "Phage"             = "#d62728",
  "CAZy"              = "#ff9896",
  "EC"                = "#9467bd",
  "GO"                = "#c5b0d5",
  "PFAM"              = "#8c564b",
  "SCFA"              = "#e377c2",
  "Symptom"           = "grey80"
)

# Remove rows with NA modality to remove NA from legend and plot
cor_lodes <- cor_lodes %>% filter(!is.na(featureType_plot))

featureType_colors <- full_palette[intersect(names(full_palette), unique(cor_lodes$featureType_plot))]
# === Create combined plot with both legends shown ===
p <- ggplot(cor_lodes,
            aes(x = axis, stratum = stratum, alluvium = alluvium, y = weight)) +
  
  # Alluvium layer with flow_group colors
  geom_alluvium(aes(fill = flow_group), width = 0.1, alpha = 0.8) +
  scale_fill_manual(values = flow_colors, name = "Association") +
  
  ggnewscale::new_scale_fill() +
  
  # Stratum layer with featureType colors
  geom_stratum(aes(fill = featureType_plot), width = 0.1, color = "black") +
  scale_fill_manual(values = featureType_colors, name = "Modality", na.translate = FALSE) + 
  
  # Labels
  geom_text(
    data = cor_lodes %>% filter(axis == "Feature"),
    stat = "stratum",
    aes(label = paste0(stratum, " →")),
    size = 5,
    color = "black",
    hjust = 1,
    nudge_x = -0.06
  ) +
  geom_text(
    data = cor_lodes %>% filter(axis == "Symptom"),
    stat = "stratum",
    aes(label = paste0("← ", stratum)),
    size = 5,
    color = "black",
    hjust = 0,
    nudge_x = 0.06
  ) +
  
  theme_minimal(base_size = 20) +
  labs(
    title = "Top Omics–Symptom Correlations (Pre–Post Changes)",
    x = NULL, y = "Correlation between pre−post delta values"
  ) +
  theme(
    panel.grid = element_blank(), 
    plot.title = element_text(size = 24, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 24, color ="black"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "top",
    legend.box = "vertical" # stack the two legends cleanly
  )

p <- p +
  scale_x_discrete(expand = expansion(add = c(0.8, 0.4)))


# === Save high-resolution PDF ===
ggsave("Curated_Sankey.pdf", p, width = 20.6, height = 20, limitsize = FALSE)
