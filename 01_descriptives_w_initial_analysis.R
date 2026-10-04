# ==============================================================================
# THESIS ANALYSIS SCRIPT
# ==============================================================================
#
# Author: Khaled Aboul Azm
#
# Purpose:
# This script performs the main statistical analyses for the thesis and
# generates the corresponding tables and figures.
#
# ==============================================================================
# ANALYSIS WORKFLOW
# ==============================================================================
#
# The complete preprocessing and instrument-selection workflow consists of:
#
# 1. 01_Datamanagement.R
#    Performs the initial data management and generates the cleaned datasets.
#
# 2. 02_Datacheck.R
#    Performs data-quality and consistency checks.
#
# 3. 03_SNP_MetaA.R
#    Performs SNP selection and associated filtering procedures.
#
# 4. 04_vitGeno_1000G.R
#    Processes and filters the vitamin D genetic data using the 1000 Genomes
#    reference data.
#
# 5. 05_Instrument_strength.R
#    evaluates the strength of the proposed genetic instruments, including
#    the calculation of the relevant F-statistics.
#
# Scripts 01_Datamanagement.R and 02_Datacheck.R have already been executed
# on the BISON server. Therefore, this script loads their resulting datasets
# directly from the project data directory.
#
# Scripts 03_SNP_MetaA.R, 04_vitGeno_1000G.R, and
# 05_Instrument_strength.R should normally be run before this script.
#
# The prerequisite scripts are executed separately and are not sourced
# automatically here. This avoids unintentionally repeating data management
# or computationally intensive genetic-data processing.
#
# When running the analysis on the GenEpi server, run this script from:
#
#   ~/ii2025MRVitD/analysis/Notebooks
#
#
#
# Notes:
# Ideally, when considering saving plots and outputs, one would save and then rm(). However, this script is a collection of all the different notebooks used
# in the analysis and it serves as a simplified way to parse through the analysis. Thus, there are no saving of any of the plots. Data might still be saved as well
# other objects that might be needed in separate scripts.
#
# ==============================================================================


# 1. INITIAL SETUP -------------------------------------------------------------

# Remove objects left over from previous R sessions
rm(list = ls())

# Set general R options
options(
  stringsAsFactors = FALSE,
  scipen = 999,
  digits = 4
)

# Set a seed for reproducible procedures involving random processes
set.seed(12345)


# 2. REQUIRED PACKAGES ---------------------------------------------------------

# Packages included in tidyverse, such as dplyr, tidyr, ggplot2, tibble, and
# purrr, do not need to be loaded separately.
required_packages <- c(
  "tidyverse",
  "lme4",
  "lmerTest",
  "hrbrthemes",
  "ExclusionTable",
  "viridis",
  "gt",
  "AER",
  "sandwich",
  "gtsummary",
  "ggdag",
  "broom",
  "nnet",
  "readxl",
  "data.table",
  "ggpubr",
  "ggExtra",
  "patchwork",
  "kableExtra"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following required packages are not installed: ",
      paste(missing_packages, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}

invisible(
  lapply(
    required_packages,
    library,
    character.only = TRUE
  )
)


# 3. PROJECT DIRECTORIES -------------------------------------------------------

# The project location can be changed through the IIVITD_PROJECT_DIR
# environment variable. The GenEpi/BISON project location is used by default.
project_dir <- Sys.getenv(
  "IIVITD_PROJECT_DIR",
  unset = path.expand("~/ii2025MRVitD")
)

analysis_dir <- file.path(project_dir, "analysis/")
notebooks_dir <- file.path(analysis_dir, "Notebooks/")
data_dir <- file.path(project_dir, "data/")
results_dir <- file.path(project_dir, "results/")

tables_dir <- file.path(results_dir, "tables/")
figures_dir <- file.path(results_dir, "figures/")
models_dir <- file.path(results_dir, "models/")


# 4. DIRECTORY CHECKS ----------------------------------------------------------

required_directories <- c(
  project_dir,
  analysis_dir,
  notebooks_dir,
  data_dir
)

missing_directories <- required_directories[
  !dir.exists(required_directories)
]

if (length(missing_directories) > 0) {
  stop(
    paste0(
      "The following required directories were not found:\n",
      paste0(" - ", missing_directories, collapse = "\n"),
      "\nCheck the project location or set IIVITD_PROJECT_DIR."
    ),
    call. = FALSE
  )
}

# Create output directories if necessary
output_directories <- c(
  results_dir,
  tables_dir,
  figures_dir,
  models_dir
)

invisible(
  lapply(
    output_directories,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

message("Project directory: ", project_dir)
message("Analysis directory: ", analysis_dir)
message("Notebook directory: ", notebooks_dir)
message("Data directory: ", data_dir)
message("Results directory: ", results_dir)




# 6. LOAD PREPROCESSED DATA ----------------------------------------------------

# Load the cleaned and imputed analysis data
vit_imp <- readRDS(paste0(data_dir, "analysisdata_imputed.rds")) %>%
  # Convert vitamin D to nmol/L and express the model variable per 12.5 nmol/L
  mutate(
    vitd25_nmol = vitd25 * 2.5,
    vitd25 = (vitd25 * 2.5) / 12.5
  )

# Combine the age variables and remove the old merge-generated columns
vit_imp <- vit_imp %>%
  mutate(
    age = coalesce(age.x, age.y)
  ) %>%
  select(
    -ends_with(".x"),
    -ends_with(".y")
  )

# Load the exclusion table and genetic data
exl_table <- read.csv(
  paste0(results_dir, "exclusion_table.txt")
)

vitGeno <- data.table::fread(
  file = paste0(data_dir, "vitGeno.tsv")
)

# Merge phenotype and genetic data
vitG <- merge(
  vit_imp,
  subset(vitGeno, select = -c(sex)),
  by.x = "id_cohort",
  by.y = "iid"
)

# Clean the full repeated-measurement dataset
vit_rep <- vitG %>%
  mutate(
    z_vit = as.numeric(scale(vitd25)),
    survey = as.factor(survey)
  ) %>%
  filter(
    abs(z_vit) < 4
  )

# Scale the principal components
vit_rep[16:25] <- scale(vit_rep[16:25])

# Store the number of observations after outlier exclusion
true_post_outliers_n <- nrow(vit_rep)

# Create a single-measurement dataset
# FU2 is preferred; baseline is used when FU2 is unavailable
vit <- vit_rep %>%
  filter(survey %in% c("B", "FU2")) %>%
  group_by(id_cohort) %>%
  arrange(
    match(survey, c("FU2", "B")),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

# Create separate baseline and FU2 datasets
vit_b <- vit_rep %>%
  filter(survey == "B")

vit_fu <- vit_rep %>%
  filter(survey == "FU2")

# Create BMI strata
vit <- vit %>%
  mutate(
    bmi_strata = cut(
      z_bmi,
      breaks = c(-Inf, -1, 1, Inf),
      labels = c(
        "Low (< -1 SD)",
        "Average (-1 to +1 SD)",
        "High (> +1 SD)"
      )
    )
  )

# Identify participants with both baseline and FU2 measurements
paired_ids <- vit_rep %>%
  filter(
    id_cohort %in% vit$id_cohort,
    survey %in% c("B", "FU2")
  ) %>%
  group_by(id_cohort) %>%
  filter(n_distinct(survey) == 2) %>%
  pull(id_cohort) %>%
  unique()

# Restrict the repeated-measurement dataset to paired participants
vit_rep <- vit_rep %>%
  filter(
    id_cohort %in% paired_ids,
    survey %in% c("B", "FU2")
  )

# Save the analytical datasets
save(
  vit,
  vit_rep,
  vit_b,
  vit_fu,
  file = paste0(data_dir, "vit.rda")
)


# Check if doesnt matter which is preffered FU2 or Baseline.
# 1. Baseline-preferred dataset (B first, FU2 fallback)
vit_b_preferred <- vit_rep %>%
  filter(survey %in% c("B", "FU2")) %>%
  group_by(id_cohort) %>%
  arrange(match(survey, c("B", "FU2")), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(
    bmi_strata = cut(
      z_bmi,
      breaks = c(-Inf, -1, 1, Inf),
      labels = c("Low (< -1 SD)", "Average (-1 to +1 SD)", "High (> +1 SD)")
    )
  )

# 2. Compare survey counts between both strategies
cat("--- FU2 Preferred Survey Counts ---\n")
print(table(vit$survey))

cat("\n--- Baseline Preferred Survey Counts ---\n")
print(table(vit_b_preferred$survey))

# 3. Cross-tabulate BMI strata changes among individuals present in both
compare_df <- inner_join(
  vit %>% select(id_cohort, strata_fu2 = bmi_strata, vitd_fu2 = vitd25),
  vit_b_preferred %>% select(id_cohort, strata_b = bmi_strata, vitd_b = vitd25),
  by = "id_cohort"
)

cat("\n--- BMI Strata Shift (FU2 vs Baseline) ---\n")
print(table(FU2_Preferred = compare_df$strata_fu2, Baseline_Preferred = compare_df$strata_b))

# 4. Calculate correlation for Vitamin D between timepoints
vitd_corr <- cor(compare_df$vitd_fu2, compare_df$vitd_b, use = "complete.obs")
cat("\nVitamin D correlation between FU2 and Baseline:", round(vitd_corr, 3), "\n")



###################################
#       Outliers
#################################


# 1. Calculate counts for the outlier exclusion (Full repeated dataset)
n_prior_outliers <- nrow(vitG)      # Raw merged dataset
n_post_outliers <- true_post_outliers_n
n_excluded_outliers <- n_prior_outliers - n_post_outliers

# 2. Calculate counts for selecting FU2 / Baseline (Reducing to 1 per subject)
n_prior_selection <- true_post_outliers_n
n_post_selection <- nrow(vit)       # Final single-survey dataset
n_excluded_selection <- n_prior_selection - n_post_selection

# 3. FIXED: Dynamic table parsing process
# Convert to character vector in case exl_table was read as a 1-column data.frame
exl_text <- if(is.data.frame(exl_table)) exl_table[[1]] else exl_table

# Dynamically identify data rows (lines starting with a number and a space)
data_lines <- exl_text[grepl("^\\s*[0-9]+\\s+", exl_text)]

# Remove the leading row numbers (1, 2, 3...) so we can parse the text cleanly
x_clean <- sub("^\\s*[0-9]+\\s+", "", data_lines)

df <- do.call(rbind, lapply(x_clean, function(line) {
  nums <- as.numeric(unlist(regmatches(line, gregexpr("[0-9]+", line))))
  text <- trimws(gsub("[0-9]+", "", line))
  data.frame(
    exclusion = text,
    n_prior = nums[1],
    n_post = nums[2],
    n_excluded = nums[3],
    stringsAsFactors = FALSE
  )
}))

df <- transform(df,
                n_prior = as.numeric(n_prior),
                n_post = as.numeric(n_post),
                n_excluded = as.numeric(n_excluded)
)

# A. Drop the old, outdated TOTAL row
# (Added trimws() to safely ensure spaces don't prevent the match)
df <- df[trimws(df$exclusion) != "TOTAL", ]

# B. Append the outlier exclusion row
outlier_row <- data.frame(
  exclusion = "Excluded as outliers",
  n_prior = n_prior_outliers,
  n_post = n_post_outliers,
  n_excluded = n_excluded_outliers,
  stringsAsFactors = FALSE
)

# C. Append the FU2/Baseline selection row
selection_row <- data.frame(
  exclusion = "Selected one survey per ID (FU2 preferred)",
  n_prior = n_prior_selection,
  n_post = n_post_selection,
  n_excluded = n_excluded_selection,
  stringsAsFactors = FALSE
)


# Bind them all together in order!
df <- rbind(df, outlier_row, selection_row)

# E. Generate a brand new, accurate TOTAL row at the very bottom
new_total_row <- data.frame(
  exclusion = "TOTAL",
  n_prior = df$n_prior[1],             # Your original starting sample size
  n_post = tail(df$n_post, 1),         # Automatically grabs the final number (nrow(vit) - 2)
  n_excluded = sum(df$n_excluded),     # Recalculates the sum including the 2 new ones
  stringsAsFactors = FALSE
)
df <- rbind(df, new_total_row)

# Clean up row names so the final output renders cleanly as 1, 2, 3...
rownames(df) <- NULL

# View final table
df




#####################
# Basic descriptives
######################

vit_small <- vit %>%
  dplyr::select(id_cohort, survey, vitd25_nmol, z_mets)

# --- 1. pairs within the same wave ---
within_wave <- vit_small %>%
  count(survey, name = "n_within") %>%
  mutate(pair_type = paste0("Exposure & Outcome in ", survey))

# --- 2. pairs across waves (baseline exposure, followup outcome) ---
cross_wave <- vit_small %>%
  pivot_wider(
    names_from = survey,
    values_from = c(vitd25_nmol, z_mets)
  ) %>%
  filter(!is.na(vitd25_nmol_B) & !is.na(z_mets_FU2)) %>%
  summarise(n_cross = n()) %>%
  mutate(pair_type = "Exposure at baseline & Outcome at follow-up")

# --- 3. combine ---
summary_table <- bind_rows(
  within_wave %>% select(pair_type, n = n_within),
  cross_wave %>% select(pair_type, n = n_cross)
)

gt(summary_table)

## Descriptive Analysis

# 1. Build the gtsummary table
tab1_ind <- vit %>%
  tbl_summary(
    include = c(age, sex, country, z_mets, vitd25_nmol, z_bmi,bmi_strata, isced, avm_1_week, pa, alc_life,
                smoke_occ, colldat_m, uvdvc_pre2, prs),
    by = survey,
    label = list(
      age        = "Age",
      sex        = "Sex",
      z_mets     = "Metabolic Syndrome (Z-Scores)",
      vitd25_nmol     = "25-hydroxyvitamin D (nmol/L)",
      z_bmi      = "Body Mass Index (Z-scores)",
      bmi_strata = "Body Mass Index (Categorical)",
      country    = "Country",
      isced      = "Parental Education",
      avm_1_week = "Audio-visual media consumption",
      pa         = "Physical activity",
      alc_life   = "Lifetime alcohol consumption",
      smoke_occ  = "Smoking status",
      colldat_m  = "Collection Date (Month)",
      uvdvc_pre2 = "UV exposure (2 months prior)",
      prs        = "PRS"
    ),
    type = list(prs ~ "continuous2"),
    statistic = list(
      prs ~ c("{median} ({p10}, {p90})",
              "{mean} ({sd})")
    )
  ) %>%
  bold_labels() %>%
  modify_header(all_stat_cols() ~ "**{level}**, N = {n}")

tab1_ind
# 2. Convert to raw LaTeX string
tex_raw <- tab1_ind %>%
  as_kable_extra(format = "latex", booktabs = TRUE) %>%
  as.character()

# 3. Robust Regex: Converts any bold <span> tag regardless of spacing
tex_clean <- gsub('<span[^>]*font-weight:\\s*bold[^>]*>(.*?)</span>', '\\\\textbf{\\1}', tex_raw)

# 4. Strip any remaining leftover HTML tags
tex_clean <- gsub('<[^>]+>', '', tex_clean)

# 5. Save clean LaTeX file directly
writeLines(
  tex_clean,
  "~/ii2025MRVitD/analysis/Notebooks/Supplementary Files/table1_characteristics.tex"
)

# Smoking and Alchol consumption are to be excluded from te analysis.


#################
# Initial Analyses
################


# 1. Define a consistent, publication-ready theme
theme_pub <- function() {
  theme_classic(base_size = 12, base_family = "sans") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.position = "top",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      # Add this to trim dead space so LaTeX handles the spacing:
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
}
# --- HISTOGRAMS ---

# Plot 1: Vitamin D (Imputed)
vitd_imputed <- ggplot(vit_imp, aes(x = vitd25_nmol)) +
  geom_histogram(bins = 30, color = "white", fill = "#3182bd", alpha = 0.85) +
  labs(
    x = "25-hydroxyvitamin D (nmol/L)",  # Updated unit
    y = "Count"
  ) +
  theme_pub()

vitd_imputed
# Plot 2 Enhanced: Histogram + Density Curve + Cutoff Annotation
vitd_wo_outliers <- ggplot(vit, aes(x = vitd25_nmol)) +
  geom_histogram(bins = 30, color = "white", fill = "#2ca25f", alpha = 0.85) +
  labs(
    x = "25-hydroxyvitamin D (nmol/L)",  # Updated unit
    y = "Count"
  ) +
  theme_pub()
vitd_wo_outliers
# Plot 3: z-Mets
z_mets_dist <- ggplot(vit, aes(x = z_mets)) +
  geom_histogram(bins = 30, color = "white", fill = "#756bb1", alpha = 0.85) +
  labs(
    x = "Metabolic Syndrome (z-score)",
    y = "Count"
  ) +
  theme_pub()
z_mets_dist


# --- SCATTERPLOT ---

# 1. Calculate the OLS beta and p-value for each survey group
beta_stats <- vit %>%
  group_by(survey) %>%
  group_modify(~ tidy(lm(z_mets ~ vitd25_nmol, data = .x))) %>%
  filter(term == "vitd25_nmol") %>%
  mutate(
    p_format = ifelse(p.value < 0.001, "< 0.001", sprintf("= %.3f", p.value)),
    plot_label = sprintf("\u03b2 = %.3f, p %s", estimate, p_format),
    
    # Set x to Inf to push it to the absolute right edge of the plot
    x_pos = Inf,
    
    # Adjust these Y values if needed so they sit cleanly in empty space
    y_pos = ifelse(survey == "B", 2.0, 1.7)
  )

# 2. Build the Scatterplot
p_scatter <- ggplot(vit, aes(x = vitd25_nmol, y = z_mets, color = survey)) +
  geom_point(size = 1.5, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = survey), alpha = 0.15, linewidth = 0.8) +
  
  # Add the labels
  geom_text(
    data = beta_stats,
    aes(x = x_pos, y = y_pos, label = plot_label, color = survey),
    # hjust = 1 right-aligns the text.
    # hjust = 1.05 adds a tiny bit of padding so it doesn't touch the axis line.
    hjust = 1.05,
    size = 4,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    name = "Survey Timepoint",
    values = c("B" = "#1f77b4", "FU2" = "#d62728"),
    labels = c("B" = "Baseline", "FU2" = "Follow-up 2")
  ) +
  scale_fill_manual(
    values = c("B" = "#1f77b4", "FU2" = "#d62728"),
    guide = "none"
  ) +
  labs(
    x = "25-hydroxyvitamin D (nmol/L)",
    y = "Metabolic Syndrome (z-score)"
  ) +
  theme_pub()

p_scatter

# Take the scatterplot object (p_scatter) created above and add marginal density plots
p_marginal <- ggMarginal(
  p_scatter,
  type = "density",     # Can also be "histogram" or "boxplot"
  groupColour = TRUE,   # Colors the margins by the survey group
  groupFill = TRUE,     # Fills the margins by the survey group
  alpha = 0.3,
  size = 4              # Controls the ratio of main plot to marginal plots
)

# To view the plot
print(p_marginal)



# 1. Classify the data and calculate percentages
vit_classified <- vit %>%
  filter(!is.na(vitd25_nmol)) %>% # Remove NAs to get accurate percentages
  mutate(
    status = case_when(
      vitd25_nmol < 25 ~ "Severe Deficiency (<25 nmol/L)",
      vitd25_nmol >= 25 & vitd25_nmol <= 50 ~ "Insufficiency (25-50 nmol/L)",
      vitd25_nmol > 50 ~ "Sufficiency (>50 nmol/L)"
    ),
    # Keep the logical order for the legend
    status = factor(status, levels = c(
      "Severe Deficiency (<25 nmol/L)",
      "Insufficiency (25-50 nmol/L)",
      "Sufficiency (>50 nmol/L)"
    ))
  )

# Calculate cohort percentages for annotations
stats <- vit_classified %>%
  count(status) %>%
  mutate(pct = round((n / sum(n)) * 100, 1))

# Extract individual percentages for easy plotting labels
pct_severe  <- stats$pct[stats$status == "Severe Deficiency (<25 nmol/L)"]
pct_insuff  <- stats$pct[stats$status == "Insufficiency (25-50 nmol/L)"]
pct_suff    <- stats$pct[stats$status == "Sufficiency (>50 nmol/L)"]

# 2. Build the Plot
vit_discret <- ggplot(vit_classified, aes(x = vitd25_nmol, fill = status)) +
  # Draw the histogram colored by clinical category
  geom_histogram(bins = 30, color = "white", alpha = 0.85) +
  
  # Add vertical boundary lines at 25 and 50 nmol/L
  geom_vline(xintercept = c(25, 50), linetype = "dashed", color = "grey40", linewidth = 0.6) +
  
  # Color palette: soft clinical red, orange, and green
  scale_fill_manual(
    name = "Clinical Status",
    values = c(
      "Severe Deficiency (<25 nmol/L)" = "#d9534f", # Muted Red
      "Insufficiency (25-50 nmol/L)" = "#f0ad4e",    # Muted Orange
      "Sufficiency (>50 nmol/L)" = "#5cb85c"         # Muted Green
    )
  ) +
  
  # Annotate percentages at the very top of each section (using y = Inf)
  annotate("text", x = 12.5, y = Inf, label = paste0(pct_severe, "%"),
           vjust = 2, size = 4.5, fontface = "bold", color = "#d9534f") +
  
  annotate("text", x = 37.5, y = Inf, label = paste0(pct_insuff, "%"),
           vjust = 2, size = 4.5, fontface = "bold", color = "#e67e22") +
  
  annotate("text", x = 75, y = Inf, label = paste0(pct_suff, "%"),
           vjust = 2, size = 4.5, fontface = "bold", color = "#4cae4c") +
  
  labs(
    x = "25-hydroxyvitamin D (nmol/L)",
    y = "Count"
  ) +
  theme_pub() +
  theme(
    legend.position = "none"
  )


vit_discret


# Based on ESPGHAN Committee on Nutrition thresholds for vitamin D status (Braegger et. al., 2013).



###############################
# Testing Associations
##############################

## Are Exposure and Outcome associated?

# We use OLS models identify all the possible confounders.

model0 <- lm(
  z_mets ~ survey + factor(country),
  data = vit
)

model1 <- lm(
  z_mets ~ vitd25 + survey + factor(country),
  data = vit
)

summary(model1)$coefficients
anova(model0, model1)



# Result: **Vitamin D** and **z_Mets** are associated.

## Confounder Identification

#Covariates: age, sex, weekly audio-visual media consumption, physical activity, lifetime alcohol consumption, smoking, colldat month, UV, ISCED , ancestry PCs

#Later we examine the association between those variables and SNPs.

# Reference model
model0 <- lm(
  z_mets ~ survey+ factor(country),
  data = vit
)

model2a <- lm(
  z_mets ~ age + survey+ factor(country),
  data = vit
)

model2s <- lm(
  z_mets ~ sex + survey+ factor(country),
  data = vit
)

model2avm <- lm(
  z_mets ~ avm_1_week + survey+ factor(country),
  data = vit
)

model2pa <- lm(
  z_mets ~ pa + survey+ factor(country),
  data = vit
)

model2uv <- lm(
  z_mets ~ uvdvc_pre2 + survey+ factor(country),
  data = vit
)

model2m <- lm(
  z_mets ~ colldat_m + survey+ factor(country),
  data = vit
)

model2i <- lm(
  z_mets ~ isced + survey+ factor(country),
  data = vit
)

model2_ancestry <- lm(
  z_mets ~ pc1 + pc2 + pc3 + pc4 + pc5 +
    pc6 + pc7 + pc8 + pc9 + pc10 + survey,
  data = vit
)

model2_bmi <- lm(
  z_mets ~ z_bmi + survey+ factor(country),
  data = vit
)


# Compare each candidate model with the reference model
covariates <- c(
  "age",
  "sex",
  "AVM",
  "PA",
  "UV",
  "colldat",
  "ISCED",
  "Ancestry (10 PCs)",
  "z_bmi"
)

res <- c(
  anova(model0, model2a)$`Pr(>F)`[2],
  anova(model0, model2s)$`Pr(>F)`[2],
  anova(model0, model2avm)$`Pr(>F)`[2],
  anova(model0, model2pa)$`Pr(>F)`[2],
  anova(model0, model2uv)$`Pr(>F)`[2],
  anova(model0, model2m)$`Pr(>F)`[2],
  anova(model0, model2i)$`Pr(>F)`[2],
  anova(model0, model2_ancestry)$`Pr(>F)`[2],
  anova(model0, model2_bmi)$`Pr(>F)`[2]
)

# FDR adjustment
fdr_adj <- p.adjust(res, method = "fdr")

# Build results table
dfm1 <- data.frame(
  covariates = covariates,
  p_val = round(res, 5),
  p_val_fdr = round(fdr_adj, 5),
  sig_fdr_5pct = fdr_adj < 0.05
)

# Order by unadjusted p-value
dfm1 <- dfm1[order(dfm1$p_val), ]

gt::gt(dfm1)


# The covariates **AVM**, **ISCED, and z_bmi** are associated with **z_mets**. Ancestry was shown to have an effect in the lmm model, while country (not shown here) has an effect in the OLS model.

## Vitamin D predictors


# Reference model
model30 <- lm(
  vitd25 ~ survey + factor(country),
  data = vit
)

# Candidate covariate models
model3a <- lm(
  vitd25 ~ age + survey + factor(country),
  data = vit
)

model3s <- lm(
  vitd25 ~ sex + survey + factor(country),
  data = vit
)

model3avm <- lm(
  vitd25 ~ avm_1_week + survey + factor(country),
  data = vit
)

model3pa <- lm(
  vitd25 ~ pa * survey + factor(country),
  data = vit
)

model3uv <- lm(
  vitd25 ~ uvdvc_pre2 + survey + factor(country),
  data = vit
)

model3m <- lm(
  vitd25 ~ colldat_m + survey + factor(country),
  data = vit
)

model3i <- lm(
  vitd25 ~ isced + survey + factor(country),
  data = vit
)

model3_ancestry <- lm(
  vitd25 ~ pc1 + pc2 + pc3 + pc4 + pc5 +
    pc6 + pc7 + pc8 + pc9 + pc10 +
    survey + factor(country),
  data = vit
)

model3_bmi <- lm(
  vitd25 ~ z_bmi + survey + factor(country),
  data = vit
)




res3 <- c(
  age = anova(model30, model3a)[["Pr(>F)"]][2],
  sex = anova(model30, model3s)[["Pr(>F)"]][2],
  avm_1_week = anova(model30, model3avm)[["Pr(>F)"]][2],
  pa = anova(model30, model3pa)[["Pr(>F)"]][2],
  uvdvc_pre2 = anova(model30, model3uv)[["Pr(>F)"]][2],
  colldat_m = anova(model30, model3m)[["Pr(>F)"]][2],
  isced = anova(model30, model3i)[["Pr(>F)"]][2],
  ancestry_pc1_pc10 = anova(model30, model3_ancestry)[["Pr(>F)"]][2],
  z_bmi = anova(model30, model3_bmi)[["Pr(>F)"]][2]
)

names(res3) <- c(
  "age",
  "sex",
  "avm_1_week",
  "pa",
  "uvdvc_pre2",
  "colldat_m",
  "isced",
  "ancestry_pc1_pc10",
  "z_bmi"
)
fdr_adj3 <- p.adjust(res3, method = "fdr")

dfm3 <- data.frame(
  covariates = names(res3),
  p_val = round(unname(res3), 5),
  p_val_fdr = round(unname(fdr_adj3), 5),
  sig_fdr_5pct = unname(fdr_adj3) < 0.05
)

dfm3 <- dfm3[order(dfm3$p_val), ]

gt::gt(dfm3)


# The covariates **sex** and **ISCED** are **NOT** associated with **vitamin D**.


##############################
# SNP selction
##############################

snp_summary <- readRDS("../../data/snp_summary.rds")
snp_summary


# Out of the full 218 SNPs provided:

# - Excluded 14 palindromic snps (AF\> 0.4)

# - 37 SNPs are missing from the Sunlight Consortium sumstats

# - 61 independent SNPs remain after LD clumping ($r^2 >0.001$)



snp_names <- c('rs11732044_t',"rs116970203_a", "prs")

covars <- c(
  "z_bmi",
  "vitd25",
  "age",
  "avm_1_week",
  "pa"
)

pcs <- paste0("pc", 1:10)

# All variables potentially needed (including ALL found SNPs)
needed_vars <- unique(c(
  "id_cohort",
  snp_names,
  covars,
  "survey",
  pcs
))

# Baseline dataset: one row per participant
vit_baseline <- vit %>%
  select(all_of(needed_vars)) %>%
  group_by(id_cohort) %>%
  slice_head(n = 1) %>%
  ungroup()


# 2. Updated Analysis Function (Accepts both covariate and specific SNP)
run_analysis_quad <- function(v, curr_snp) {
  
  # Skip non-numeric outcomes
  if (!is.numeric(vit[[v]])) {
    return(
      tibble(
        snp            = curr_snp,
        covariate      = v,
        p_base_total   = NA_real_,
        p_base_direct  = NA_real_,
        p_full_total   = NA_real_,
        p_full_direct  = NA_real_
      )
    )
  }
  
  adj_vars <- c("age", "survey", pcs)
  
  # Remove outcome if present
  adj_vars <- setdiff(adj_vars, v)
  
  # Direct model includes vitd25 unless outcome IS vitd25
  direct_vars <- adj_vars
  if (v != "vitd25") {
    direct_vars <- c("vitd25", direct_vars)
  }
  
  # Formulas dynamically using the current SNP
  form_total <- as.formula(
    paste(v, "~", paste(c(curr_snp, adj_vars), collapse = " + "))
  )
  
  form_direct <- as.formula(
    paste(v, "~", paste(c(curr_snp, direct_vars), collapse = " + "))
  )
  
  get_p <- function(formula, data) {
    tryCatch({
      fit <- lm(formula, data = data)
      broom::tidy(fit) %>%
        filter(term == curr_snp) %>%
        pull(p.value) %>%
        first()
    }, error = function(e) {
      NA_real_
    })
  }
  
  tibble(
    snp            = curr_snp,
    covariate      = v,
    p_base_total   = get_p(form_total, vit_baseline),
    p_base_direct  = get_p(form_direct, vit_baseline),
    p_full_total   = get_p(form_total, vit),
    p_full_direct  = get_p(form_direct, vit)
  )
}


# 3. Create a grid of all combinations and run the analysis
analysis_grid <- expand_grid(v = covars, curr_snp = snp_names)
final_results <- pmap_dfr(analysis_grid, run_analysis_quad)


# 4. Generate the GT Table
p_cols <- c(
  "p_base_total",
  "p_base_direct",
  "p_full_total",
  "p_full_direct"
)

tab <- final_results %>%
  group_by(snp) %>%  # Groups the table visually by each SNP
  gt() %>%
  tab_header(
    title = "SNP Association (Outcome ~ SNP)",
    subtitle = "Categorical outcomes skipped"
  ) %>%
  cols_label(
    covariate = "Outcome Variable"
  ) %>%
  fmt_number(
    columns = all_of(p_cols),
    decimals = 3
  ) %>%
  sub_missing(missing_text = "-")

# Bold significant p-values safely
for (col in p_cols) {
  tab <- tab %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = all_of(col),
        rows = .data[[col]] <= 0.05
      )
    )
}

tab


# We test the associations in the form of outcome \~ SNP + age + survey +pcs + X . The results are represented across 4 cases:

#  - p_base_total: Cross-sectional subset at baseline where X = null

# - p_base_direct: Cross-sectional subset at baseline where X = Vitamin D

# - p_full_total: Full longitudinal data where X = null

# - p_full_direct: Full longitudinal data where X = Vitamin D

analysis_data <- vit %>%
  dplyr::select(all_of(snp_names), survey, age, sex, country, isced) %>%
  mutate(
    across(all_of(snp_names), as.numeric),
    isced   = ordered(isced),
    sex     = factor(sex),
    country = factor(country),
    survey  = factor(survey)
  ) %>%
  droplevels()

# =====================================================
# VARIABLE TYPES
# =====================================================

ordered_vars <- c("isced")
binary_vars  <- c("sex")
nominal_vars <- c("country")

# =====================================================
# MODEL FUNCTION
# =====================================================

run_model <- function(outcome, curr_snp) {
  
  # Build covariates using the current SNP
  covars <- c(curr_snp, "survey", "age", "sex", "country")
  
  # Remove outcome from covariates if it's in the list
  covars <- setdiff(covars, outcome)
  
  f <- as.formula(
    paste(outcome, "~", paste(covars, collapse = " + "))
  )
  
  model_data <- model.frame(
    f,
    data = analysis_data,
    na.action = na.omit
  )
  
  if (outcome %in% ordered_vars) {
    
    fit <- MASS::polr(f, data = model_data, Hess = TRUE)
    
    coefs <- coef(summary(fit))
    pvals <- 2 * pnorm(abs(coefs[, "t value"]), lower.tail = FALSE)
    
    out <- tibble(
      outcome = outcome,
      comparison = NA_character_,
      term = rownames(coefs),
      estimate = coefs[, "Value"],
      std.error = coefs[, "Std. Error"],
      statistic = coefs[, "t value"],
      p.value = pvals,
      model = "ordinal_logistic"
    )
    
  } else if (outcome %in% binary_vars) {
    
    y <- droplevels(factor(model_data[[outcome]]))
    if (nlevels(y) < 2) return(tibble())
    
    fit <- glm(f, data = model_data, family = binomial)
    levs <- levels(y)
    
    out <- tidy(fit) %>%
      dplyr::mutate(
        outcome = .env$outcome,
        comparison = paste0(levs[2], " vs ", levs[1]),
        model = "binary_logistic"
      ) %>%
      dplyr::select(outcome, comparison, term, estimate, std.error, statistic, p.value, model)
    
  } else if (outcome %in% nominal_vars) {
    
    y <- droplevels(factor(model_data[[outcome]]))
    if (nlevels(y) < 2) return(tibble())
    
    fit <- nnet::multinom(f, data = model_data, trace = FALSE)
    s <- summary(fit)
    
    # Handle binary outcomes gracefully in multinomial logic (returns vector instead of matrix)
    if (is.null(dim(s$coefficients))) {
      coefs <- matrix(s$coefficients, nrow = 1, dimnames = list(rownames(s$standard.errors), names(s$coefficients)))
      ses   <- matrix(s$standard.errors, nrow = 1, dimnames = list(rownames(s$standard.errors), names(s$standard.errors)))
    } else {
      coefs <- as.matrix(s$coefficients)
      ses   <- as.matrix(s$standard.errors)
    }
    
    zvals <- coefs / ses
    pvals <- 2 * (1 - pnorm(abs(zvals)))
    
    out <- tibble(
      outcome = outcome,
      comparison = rep(rownames(coefs), each = ncol(coefs)),
      term = rep(colnames(coefs), times = nrow(coefs)),
      estimate = as.vector(t(coefs)),
      std.error = as.vector(t(ses)),
      statistic = as.vector(t(zvals)),
      p.value = as.vector(t(pvals)),
      model = "multinomial_logistic"
    )
    
  } else {
    out <- tibble()
  }
  
  # FIX 2: Filter specifically for the current SNP term
  out %>%
    filter(term == curr_snp) %>%
    dplyr::mutate(
      snp = curr_snp # Add a column to identify which SNP is being tested
    )
}

# =====================================================
# RUN MODELS
# =====================================================

all_outcomes <- c(ordered_vars, binary_vars, nominal_vars)

# Create grid of all outcomes x all SNPs
analysis_grid <- expand_grid(outcome = all_outcomes, curr_snp = snp_names)

# Map over the grid
results <- pmap_dfr(analysis_grid, run_model)

# =====================================================
# MULTIPLE TESTING + ODDS RATIOS
# =====================================================

snp_results <- results %>%
  group_by(snp) %>% # Calculate FDR within each SNP's family of tests
  dplyr::mutate(
    p.adjusted = p.adjust(p.value, method = "fdr"),
    odds_ratio = exp(estimate),
    OR_low = exp(estimate - 1.96 * std.error),
    OR_high = exp(estimate + 1.96 * std.error)
  ) %>%
  ungroup() %>%
  arrange(snp, outcome, model, term)

# =====================================================
# FINAL TABLE
# =====================================================

final_table <- snp_results %>%
  dplyr::select(
    snp,
    outcome,
    comparison,
    model,
    odds_ratio,
    OR_low,
    OR_high,
    p.value,
    p.adjusted
  )

# =====================================================
# DISPLAY
# =====================================================

final_table %>%
  group_by(snp) %>% # Group visually by SNP in the GT table
  gt() %>%
  tab_header(
    title = "Predictor Associations With Categorical Outcomes",
    subtitle = "Models adjusted for covariates; P-values FDR-corrected per predictor"
  ) %>%
  fmt_number(
    columns = c(odds_ratio, OR_low, OR_high, p.value, p.adjusted),
    decimals = 3
  ) %>%
  cols_label(
    outcome = "Outcome",
    comparison = "Comparison",
    model = "Model",
    odds_ratio = "OR",
    OR_low = "CI Low",
    OR_high = "CI High",
    p.value = "P",
    p.adjusted = "FDR P"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p.adjusted,
      rows = p.adjusted <= 0.05
    )
  ) %>%
  sub_missing(missing_text = "-") %>%
  opt_row_striping()


#In theory, Ancestry should capture all the effects that the variable country can have on the SNP.
# This might be both a blessing and a curse, as its far easier to falsify the SNP given a discrete country 
# variable than using the continuous ancestry PCs.

