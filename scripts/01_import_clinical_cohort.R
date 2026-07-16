#import clinical cohort data


library(readr)

Clinical_data <- read_csv(
  "data/raw/clinical_cohort.csv"
  )

#Check structure
glimpse(Clinical_data)


#View first rows
head(Clinical_data)


