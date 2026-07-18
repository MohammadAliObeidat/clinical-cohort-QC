# ============================================================
# Clinical Cohort Import & QC Report
# Script: 01_import_clinical_cohort.R
# Purpose: Import raw clinical cohort and perform initial QC
# ============================================================

# ------------------------------------------------------------
# Import dataset
# ------------------------------------------------------------

library(tidyverse)

library(readr)

Clinical_data <- read_csv(
  "data/raw/clinical_cohort.csv"
  )

# ------------------------------------------------------------
#  Five-check reflex
# ------------------------------------------------------------

#Check structure
glimpse(Clinical_data)


#View first rows
head(Clinical_data)

# Dataset dimensions
dim(Clinical_data)


# Column names
names(Clinical_data)


# Basic summary
summary(Clinical_data)

# -----------------------------
# Missing data assessment
# -----------------------------


# Count missing values per column
missing_summary <- colSums(is.na(Clinical_data))

missing_summary

# Convert missingness results to a table
missing_table <- data.frame(
  variable = names(missing_summary),
  missing_count = as.numeric(missing_summary)
)

missing_table

write_csv(
  missing_table,
  "results/tables/missing_data_summary.csv"
)

#-------------------------------------
#Disease Status QC
#-------------------------------------

disease_status_counts <- table(Clinical_data$disease_status)
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

#-------------------------------------
# Sex QC
#-------------------------------------

sex_counts  <- table(Clinical_data$sex)
sex_counts

sex_table <- data.frame(
  sex = names(sex_counts),
  count = as.numeric(sex_counts)
)

write.csv(
  sex_table,
  "results/tables/sex_summary.csv"
)

# Disease status by sex cross-tabulation

disease_sex_table <- table(
  Clinical_data$disease_status,
  Clinical_data$sex
)

disease_sex_table_df <- as.data.frame(disease_sex_table)

write_csv(
  disease_sex_table_df,
  "results/tables/disease_status_by_sex.csv"
)

#-------------------------------------
# Age QC
#-------------------------------------

summary(Clinical_data$age)

age_summary <- summary(Clinical_data$age)
age_summary_table <- data.frame(
  statistic = names(age_summary),
  value = as.numeric(age_summary)
)
write.csv(age_summary_table, "results/tables/age_summary.csv")

#-------------------------------------
#Age by disease status
#-------------------------------------

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

#-------------------------------------
# Save mean age by disease status
#-------------------------------------

mean_age_by_status <- aggregate(
  age ~ disease_status,
  data = Clinical_data,
  FUN = mean
)

write_csv(
  mean_age_by_status,
  "results/tables/mean_age_by_disease_status.csv"
)
