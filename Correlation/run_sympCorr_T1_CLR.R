###############################################
# Remove all the variables from the Workspace #
###############################################

rm(list = ls(all = TRUE))
ls()

#########################
# Set working directory #
#########################

workingDirectory <- '/Users/himelmallick/Library/CloudStorage/Dropbox/Insanity/Projects/Projects_Completed/ASU_PTHS/'
setwd(workingDirectory)

##################
# Load libraries #
##################

# Core data wrangling
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readxl)
library(tidyverse)
library(compositions)

# Visualization
library(ggplot2)
library(ggalluvial)
library(scales)
library(ggnewscale)

# Correlation + progress
library(pbapply)

###################
# Load All Layers #
###################

# Get all .RData files except 'Symptoms.RData'
file_paths <- list.files("./Data_Processed", pattern = "\\.RData$", full.names = TRUE)
file_paths <- file_paths[!grepl("Symptoms\\.RData$", file_paths)]

# Load each 'pcl' object and assign it to a named variable
for (path in file_paths) {
  var_name <- paste0("pcl_", tools::file_path_sans_ext(basename(path)))
  assign(var_name, get(load(path)))  # load returns the name of the loaded object
}

#################
# Load metadata #  
#################

metadata_full <- read_xlsx('Data/MASTER_OMICS_data_05302025.xlsx', sheet = 'Metadata')
metadata_full <- column_to_rownames(metadata_full, 'SampleID')

####################
# Cleanup metadata #
####################

metadata_full <- metadata_full %>% filter(Criteria_Longitudinal=='Include') %>% filter(Note=='NA') 

##################################
# Subset to subjects of interest #
##################################

metadata <- metadata_full %>% filter(ReGroup %in% c('Baseline', 'EndMTT'))

# Define Time_point based on ReGroup
metadata$Time_point <- ifelse(metadata$ReGroup == 'Baseline', 'T0', 'T1') 
metadata$Y <- ifelse(metadata$ReGroup == 'Baseline', 0, 1)

######################################
# PREPARE DATA FOR INTEGRATEDLEARNER #
######################################

###########################################
# Dynamically define all available layers #
###########################################

# Find all variables starting with 'pcl_'
pcl_vars <- ls(pattern = "^pcl_")

# Remove 'pcl_' prefix and create a mapping
layer_names <- gsub("^pcl_PTHS_", "", pcl_vars)
layer_map <- setNames(pcl_vars, layer_names)  # names(layer_map) are clean layer names


# Compostional Adjustment
# Set your prevalence threshold (e.g., present in at least 20% of samples)
prevalence_threshold <- 0.2
sample_ids <- rownames(metadata)
list_features <- list()

for (layer in names(layer_map)) {
  obj <- get(layer_map[[layer]])
  layer_data <- obj$features
  
  # Restrict to samples in metadata and reorder
  layer_data <- layer_data[sample_ids, , drop = FALSE]
  
  # Prevalence filtering (before CLR)
  prevalence <- colMeans(layer_data > 0)
  keep_features <- prevalence >= prevalence_threshold
  layer_data <- layer_data[, keep_features, drop = FALSE]
  
  message(sprintf("Layer %s: retained %d/%d features after %.0f%% prevalence filtering.",
                  layer, sum(keep_features), length(keep_features), prevalence_threshold * 100))
  
  # CLR transformation (add pseudocount to avoid log(0))
  layer_data <- layer_data + min(layer_data[layer_data > 0])
  layer_data <- compositions::clr(layer_data)
  
  # Transpose to feature × sample
  layer_data_t <- as.data.frame(t(layer_data))
  
  # Check sample consistency
  stopifnot(all(colnames(layer_data_t) == sample_ids))
  
  list_features[[layer]] <- layer_data_t
}

#####################################
# Merge all features into one table #
#####################################

feature_table <- as.data.frame(Reduce(rbind, list_features))

##########################################################
# Reformat metadata table according to IntegratedLearner #
##########################################################

sample_metadata <- metadata %>%
  dplyr::rename(subjectID = `Patient ID`)

################################################
# Make sure features and metadata are in order #
################################################

stopifnot(all(rownames(sample_metadata) == colnames(feature_table)))

######################################
# Create metadata table for features #
######################################

rowID <- rep(names(list_features), sapply(list_features, nrow))
feature_metadata <- cbind.data.frame(
  featureID = rownames(feature_table),
  featureType = rowID
)
rownames(feature_metadata) <- feature_metadata$featureID

######################
# Final sanity check #
######################

stopifnot(all(rownames(feature_metadata) == rownames(feature_table)))
stopifnot(all(rownames(sample_metadata) == colnames(feature_table)))

############## 
# Remove NAs #
##############

feature_table <- na.omit(feature_table)
newIDs <- rownames(feature_table)
feature_metadata <- feature_metadata[newIDs, ]

# Final consistency checks
stopifnot(all(rownames(feature_metadata) == rownames(feature_table)))
stopifnot(all(rownames(sample_metadata) == colnames(feature_table)))

######################
# Clean up workspace #
######################

rm(list = setdiff(ls(), c("feature_table", "feature_metadata", "sample_metadata")))

######################
# Load symptoms data #
######################

load("./Data_Processed/PTHS_Symptoms.RData")
sympdata<-pcl$features; rm(pcl)
sympdata<-sympdata[complete.cases(sympdata),]

######################################
# Merge to get consistent sample IDs #
######################################

sample_metadata<-rownames_to_column(sample_metadata, 'sampleID')
new_symp_metadata<-merge(sample_metadata, sympdata, 'sampleID')
new_symp_metadata <- new_symp_metadata %>% dplyr::select(-ends_with(".y"))
names(new_symp_metadata) <- sub("\\.x$", "", names(new_symp_metadata))

####################################################
# New Analysis: Delta Feature vs. Delta Symptom    #
# Correlation + Multi-omics Alluvial Visualization #
####################################################

# Compute delta symptom scores
symptom_wide <- new_symp_metadata %>%
  filter(!is.na(Score)) %>%
  pivot_wider(id_cols = c(subjectID, Symptom), 
              names_from = Time_point, values_from = Score) %>%
  filter(!is.na(T0) & !is.na(T1)) %>%
  mutate(delta_score = T1 - T0)

###########################################
# Compute delta features by subjectID     #
###########################################

sample_metadata<-column_to_rownames(sample_metadata, 'sampleID')

# Align sample IDs and transpose feature table
stopifnot(ncol(feature_table) == nrow(sample_metadata))
colnames(feature_table) <- rownames(sample_metadata)

# Transpose to sample × feature
omics_df <- as.data.frame(t(feature_table))
omics_df$sampleID <- rownames(omics_df)
sample_metadata <- rownames_to_column(sample_metadata, "sampleID")

# Join with metadata
omics_long <- sample_metadata %>%
  select(sampleID, subjectID, Time_point) %>%
  left_join(omics_df, by = "sampleID") %>%
  pivot_longer(-c(sampleID, subjectID, Time_point), names_to = "feature", values_to = "value")

# Ensure Time_point is factor with defined levels
omics_long$Time_point <- factor(omics_long$Time_point, levels = c("T0", "T1"))

# Compute delta only for subject-feature pairs with both timepoints
feature_deltas <- omics_long %>%
  group_by(subjectID, feature, Time_point) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(Time_point = trimws(toupper(as.character(Time_point)))) %>%
  filter(Time_point %in% c("T0", "T1")) %>%
  mutate(Time_point = factor(Time_point, levels = c("T0", "T1"))) %>%
  pivot_wider(names_from = Time_point, values_from = value) %>%
  filter(!is.na(T0) & !is.na(T1)) %>%
  mutate(delta_feature = T1 - T0)

############################################
# Correlation Testing and FDR Adjustment   #
# + Keep All Results, Filter Later         #
############################################

# Apply Spearman correlation
safe_cor_test <- function(df) {
  res <- tryCatch({
    test <- cor.test(df$delta_feature, df$delta_score, method = "spearman")
    tibble(correlation = as.numeric(test$estimate), p_value = test$p.value)
  }, error = function(e) {
    tibble(correlation = cor(df$delta_feature, df$delta_score, method = "spearman"), p_value = NA_real_)
  })
  return(res)
}

# Merge delta feature + symptom data
merged_df <- feature_deltas %>%
  inner_join(symptom_wide %>% dplyr::select(subjectID, Symptom, delta_score),
             by = "subjectID")

# Split by feature × symptom and apply correlation test
df_list <- merged_df %>%
  group_by(feature, Symptom) %>%
  filter(n() >= 3) %>%
  group_split()

# Apply correlation test with progress bar
results_list <- pbapply::pblapply(df_list, function(df) {
  stat <- safe_cor_test(df)
  tibble(
    feature = unique(df$feature),
    Symptom = unique(df$Symptom),
    correlation = stat$correlation,
    p_value = stat$p_value
  )
})

# Combine results and calculate FDR (no early filtering)
cor_results <- bind_rows(results_list) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  left_join(feature_metadata, by = c("feature" = "featureID"))

# Optional: label raw p-value significance
cor_results <- cor_results %>%
  mutate(p_label = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))


##################################
# Wilcoxon tests (nonparametric) #
##################################

# Wilcoxon test for delta features
feature_stats <- feature_deltas %>%
  group_by(feature) %>%
  summarise(p_value_feature = tryCatch(
    wilcox.test(delta_feature, mu = 0)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Wilcoxon test for delta symptoms
symptom_stats <- symptom_wide %>%
  group_by(Symptom) %>%
  summarise(p_value_symptom = tryCatch(
    wilcox.test(delta_score, mu = 0)$p.value, error = function(e) NA_real_),
    .groups = "drop")

# Merge p-values into results
cor_results <- cor_results %>%
  left_join(feature_stats, by = "feature") %>%
  left_join(symptom_stats, by = "Symptom") %>%
  mutate(
    sig_feature = p_value_feature < 0.05,
    sig_symptom = p_value_symptom < 0.05
  )


###############################################
# Select top features per modality + overall  #
###############################################

top_per_layer <- cor_results %>%
  group_by(featureType) %>%
  arrange(p_value, desc(abs(correlation))) %>%
  slice_head(n = 5) %>%
  ungroup()

# Take union of top 5 per layer, then pick top 30 by p-value
top30_all <- top_per_layer %>%
  arrange(p_value, desc(abs(correlation))) %>%
  slice_head(n = 30)

##################
# Label cleaning #
##################

extract_last_taxon <- function(taxon_string) {
  matches <- str_extract_all(taxon_string, "[a-z]__[^.;_]+")[[1]]
  if (length(matches) == 0) return(NA_character_)
  str_remove(tail(matches, 1), "^[a-z]__")
}

cor_results_filtered <- top30_all %>%
  mutate(
    direction = ifelse(correlation > 0, "Positive", "Negative"),
    feature_display = case_when(
      featureType == "EC" ~ str_remove(feature, "^\\d+\\.\\d+\\.\\d\\.\\d+:"),  # Removes EC number
      grepl("^d__.*__", feature) ~ map_chr(feature, extract_last_taxon),
      TRUE ~ str_trunc(feature, 30)
    ),
    Symptom_display = str_to_title(
      str_remove(as.character(Symptom), "\\s*\\(.*\\)")
    ) %>% str_trunc(30)
  )
# Prepare for plotting
cor_long <- cor_results_filtered %>%
  mutate(
    correlation_abs = abs(correlation),
    alluvium_id = paste(feature_display, Symptom_display, sep = "_"),
    featureType_plot = featureType
  )

cor_lodes <- bind_rows(
  cor_long %>%
    transmute(
      axis = "Feature",
      stratum = feature_display,
      alluvium = alluvium_id,
      correlation_abs,
      direction,
      featureType_plot
    ),
  cor_long %>%
    transmute(
      axis = "Symptom",
      stratum = Symptom_display,
      alluvium = alluvium_id,
      correlation_abs,
      direction,
      featureType_plot = "Symptom"
    )
)

########################
# Define color palettes
########################

full_palette <- c(
  "16S_Bacteria"      = "#1f77b4",
  "Shotgun_Bacteria"  = "#aec7e8",
  "ITS_Fungi"         = "#ff7f0e",
  "Shotgun_Fungi"     = "#ffbb78",
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

direction_colors <- c(
  "Positive" = "#E69F00",
  "Negative" = "#56B4E9"
)

featureType_colors <- full_palette[names(full_palette) %in% unique(cor_lodes$featureType_plot)]

#####################
# Plot alluvial chart
#####################

p <- ggplot(cor_lodes,
            aes(x = axis, stratum = stratum, alluvium = alluvium, y = correlation_abs)) +
  geom_alluvium(aes(fill = direction), width = 1/12, alpha = 0.8) +
  scale_fill_manual(values = direction_colors, guide = guide_legend(title = "Direction")) +
  new_scale_fill() +
  geom_stratum(aes(fill = featureType_plot), width = 1/12, color = "black") +
  scale_fill_manual(values = featureType_colors, guide = guide_legend(title = "Modality")) +
  geom_text(stat = "stratum", aes(label = paste0(stratum, " →")), size = 4.5,
            hjust = 1, data = cor_lodes %>% filter(axis == "Feature"), nudge_x = -0.03) +
  geom_text(stat = "stratum", aes(label = paste0("← ", stratum)), size = 4.5,
            hjust = 0, data = cor_lodes %>% filter(axis == "Symptom"), nudge_x = 0.03) +
  theme_minimal() +
  labs(
    title = "Top Omics–Symptom Correlations (10 wks vs. Baseline)",
    x = NULL, y = "Correlation between pre−post delta values"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # 🔼 Title size here
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 14),  # ← This increases y-axis title
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )

# Save output
# pdf("./Analysis/Analysis_Symptoms/Results/Top_Omics_Symptom_Correlations_T1_CLR.pdf", width = 15, height = 8)
# print(p)
# dev.off()


# Compute mean delta (T1 - T0) per feature
feature_mean_deltas <- feature_deltas %>%
  group_by(feature) %>%
  summarise(mean_delta = mean(delta_feature, na.rm = TRUE), .groups = "drop")

# Compute mean delta (T1 - T0) per symptom
symptom_mean_deltas <- symptom_wide %>%
  group_by(Symptom) %>%
  summarise(mean_delta_score = mean(delta_score, na.rm = TRUE), .groups = "drop")

# Join with correlation results
cor_results_export <- cor_results %>%
  filter(p_value < 0.001, abs(correlation) > 0.8) %>%
  left_join(feature_mean_deltas, by = "feature") %>%
  left_join(symptom_mean_deltas, by = "Symptom") %>%
  mutate(direction = ifelse(correlation > 0, "Positive", "Negative"))

# Write to CSV
write.csv(cor_results_export, "./Analysis/Analysis_Symptoms/Results/Top_Omics_Symptom_Correlations_T1_CLR.csv", row.names = FALSE)