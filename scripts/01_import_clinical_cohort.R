# ============================================================
# Clinical Cohort Import & QC Report
# Script: 01_import_clinical_cohort.R
# Purpose: Import the raw clinical cohort, perform initial
# quality control (QC), remove duplicate sample IDs,
# and generate a processed dataset for downstream analysis.
# ============================================================


# ------------------------------------------------------------
# Import dataset
# ------------------------------------------------------------

library(tidyverse)

Clinical_data <- read_csv(
  "data/raw/Clinical_cohort.csv"
)


# ------------------------------------------------------------
# Five-check reflex
# ------------------------------------------------------------

# Check structure
glimpse(Clinical_data)

# View first rows
head(Clinical_data)

# Dataset dimensions
dim(Clinical_data)

# Column names
names(Clinical_data)

# Basic summary
summary(Clinical_data)


#------------------------------------
# Missing Data Assessment
#------------------------------------

# Count missing values per column
missing_summary <- colSums(is.na(Clinical_data))

missing_summary


# Convert missingness results into a table
missing_table <- data.frame(
  variable = names(missing_summary),
  missing_count = as.numeric(missing_summary)
)

missing_table


write_csv(
  missing_table,
  "results/tables/missing_data_summary.csv"
)



#------------------------------------
# Disease Status QC
#------------------------------------

disease_status_counts <- table(
  Clinical_data$disease_status
)

disease_status_counts


# Convert disease status counts into a table
disease_status_table <- data.frame(
  disease_status = names(disease_status_counts),
  count = as.numeric(disease_status_counts)
)


write_csv(
  disease_status_table,
  "results/tables/disease_status_summary.csv"
)



#------------------------------------
# Sex QC
#------------------------------------

sex_counts <- table(
  Clinical_data$sex
)

sex_counts


sex_table <- data.frame(
  sex = names(sex_counts),
  count = as.numeric(sex_counts)
)


write_csv(
  sex_table,
  "results/tables/sex_summary.csv"
)


# Disease status by sex cross-tabulation

disease_sex_table <- table(
  Clinical_data$disease_status,
  Clinical_data$sex
)


disease_sex_table_df <- as.data.frame(
  disease_sex_table
)


write_csv(
  disease_sex_table_df,
  "results/tables/disease_status_by_sex.csv"
)



#------------------------------------
# Age QC
#------------------------------------

summary(Clinical_data$age)


age_summary <- summary(
  Clinical_data$age
)


age_summary_table <- data.frame(
  statistic = names(age_summary),
  value = as.numeric(age_summary)
)


write_csv(
  age_summary_table,
  "results/tables/age_summary.csv"
)



#------------------------------------
# Age by Disease Status
#------------------------------------

age_summary_by_status <- do.call(
  data.frame,
  aggregate(
    age ~ disease_status,
    data = Clinical_data,
    FUN = summary
  )
)


write_csv(
  age_summary_by_status,
  "results/tables/age_summary_by_disease_status.csv"
)



#------------------------------------
# Mean Age by Disease Status
#------------------------------------

mean_age_by_status <- aggregate(
  age ~ disease_status,
  data = Clinical_data,
  FUN = mean
)


write_csv(
  mean_age_by_status,
  "results/tables/mean_age_by_disease_status.csv"
)



#------------------------------------
# Sample ID QC
#------------------------------------


# Total number of samples

total_samples <- nrow(Clinical_data)

total_samples


# Number of unique sample IDs

unique_samples <- n_distinct(
  Clinical_data$sample_id
)

unique_samples


# Missing sample IDs

missing_sample_ids <- sum(
  is.na(Clinical_data$sample_id)
)

missing_sample_ids



# Detect duplicate sample IDs

duplicate_sample_ids <- Clinical_data %>%
  count(sample_id) %>%
  filter(n > 1)


duplicate_sample_ids


write_csv(
  duplicate_sample_ids,
  "results/tables/duplicate_sample_ids.csv"
)



# Create sample ID summary table

sample_id_summary <- data.frame(
  metric = c(
    "total rows",
    "unique sample IDs",
    "missing sample IDs",
    "duplicate sample IDs"
  ),
  value = c(
    total_samples,
    unique_samples,
    missing_sample_ids,
    nrow(duplicate_sample_ids)
  )
)


sample_id_summary


write_csv(
  sample_id_summary,
  "results/tables/sample_id_summary.csv"
)



# Inspect duplicate records

Clinical_data %>%
  filter(
    sample_id %in% duplicate_sample_ids$sample_id
  ) %>%
  arrange(sample_id)



#------------------------------------
# Create Processed Dataset
#------------------------------------

Clinical_data_clean <- Clinical_data %>%
  distinct(
    sample_id,
    .keep_all = TRUE
  )


write_csv(
  Clinical_data_clean,
  "data/processed/Clinical_cohort_clean.csv"
)



#------------------------------------
# Clinical Variable Validation
#------------------------------------


#------------------------------------
# Age QC
#------------------------------------

# Check age data type
class(Clinical_data_clean$age)

# Check age distribution and summary statistics
summary(Clinical_data_clean$age)

# Check missing age values
sum(is.na(Clinical_data_clean$age))

# Identify suspicious age values
Clinical_data_clean %>%
  filter(age > 100)


# QC note:
# CTRL059 has an age value of 141 years.
# This value is biologically suspicious and should be reviewed.
# The original value is retained and not modified.



#------------------------------------
# Sex QC
#------------------------------------

# Check sex data type
class(Clinical_data_clean$sex)

# Check missing sex values
sum(is.na(Clinical_data_clean$sex))

# Check sex categories and counts
table(Clinical_data_clean$sex)


# QC note:
# Sex variable contains two consistent categories:
# F = 43 samples
# M = 23 samples
# No missing values detected.



#------------------------------------
# Disease Status QC
#------------------------------------

# Check disease status data type
class(Clinical_data_clean$disease_status)

# Check missing values
sum(is.na(Clinical_data_clean$disease_status))

# Check disease status categories
table(Clinical_data_clean$disease_status)


# QC note:
# Disease status contains inconsistent capitalization.
# Case samples:
# case = 19
# Case = 7
# CASE = 8
#
# Control samples:
# control = 22
# Control = 5
# CONTROL = 5
#
# These categories should be standardized before
# downstream analysis.



#------------------------------------
# Variant Count QC
#------------------------------------

# Check variant count data type
class(Clinical_data_clean$variant_count)

# Check missing values
sum(is.na(Clinical_data_clean$variant_count))

# Check variant count distribution
table(Clinical_data_clean$variant_count)

# Summary statistics
summary(Clinical_data_clean$variant_count)


# QC note:
# Variant count is numeric with no missing values.
# Values range from 0 to 8 variants per sample.
# No negative or suspicious values detected.



#------------------------------------
# Onset Age QC
#------------------------------------

# Check onset age data type
class(Clinical_data_clean$onset_age)

# Check missing values
sum(is.na(Clinical_data_clean$onset_age))

# Check onset age distribution
table(Clinical_data_clean$onset_age)

# Summary statistics
summary(Clinical_data_clean$onset_age)


# Investigate missing onset ages

Clinical_data_clean %>%
  filter(is.na(onset_age)) %>%
  count(disease_status)


# QC note:
# Missing onset_age values are only present in control samples.
# This is expected because controls do not have disease onset.
# Case samples have recorded onset ages.



# ============================================================
# End of Script
#
# Clinical cohort successfully imported, quality controlled,
# and processed for downstream analyses.
#
# Next script:
# 02_data_cleaning.R
# ============================================================