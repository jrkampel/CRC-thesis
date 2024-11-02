# CLINICAL DATA ANALYSIS 
# Load Libs
library(readxl)
library(dplyr)
library(readr)
install.packages("writexl")
library(writexl)


clin_data_TCGA <- read.delim("C:/Users/jrkam/OneDrive - University of Cape Town/Honours/thesis/Finding genes/coadread_tcga_pan_can_atlas_2018_clinical_data.tsv", header = TRUE, sep = "\t")


clin_data_SA <- readxl::read_excel("C:/Users/jrkam/OneDrive - University of Cape Town/Honours/thesis/Finding genes/clinical_data.xlsx")

# Calculate the average of Diagnosis.Age
average_age_SA <- mean(clin_data_SA$Age, na.rm = TRUE)

# Print the average age
average_age_SA

# Assuming your data is in a data frame named df
table(clin_data_SA$pT)

# Store the counts in a data frame
stage_counts <- as.data.frame(table(clin_data_SA$pT))

# View the data frame
print(stage_counts)

# Frequency table of stages
stage_counts <- table(clin_data_SA$pT)

# Total number of patients
total_patients <- sum(stage_counts)

# Subset for stages 3, 4a, and 4b
advanced_stages <- stage_counts[c("3", "4a", "4b")]

# Total number of patients in stages 3, 4a, and 4b
advanced_total <- sum(advanced_stages)

# Calculate the percentage
percentage_advanced <- (advanced_total / total_patients) * 100

# Print the result
percentage_advanced

# Create a matrix of gender data
gender_data <- matrix(c(51, 49, 56.3, 43.8), nrow = 2, byrow = TRUE)
colnames(gender_data) <- c("Male", "Female")
rownames(gender_data) <- c("Porto", "South_African")

# Perform chi-square test
chisq.test(gender_data)
# Create a matrix of CRC localization data
crc_localization <- matrix(c(65, 35, 65.6, 34.4), nrow = 2, byrow = TRUE)
colnames(crc_localization) <- c("Proximal", "Distal")
rownames(crc_localization) <- c("Porto", "South_African")

# Perform chi-square test
chisq.test(crc_localization)

# Create a matrix for TNM stages
tnm_stage_data <- matrix(c(68, 32, 6.25, 93.75), nrow = 2, byrow = TRUE)

# Add column and row names for better understanding
colnames(tnm_stage_data) <- c("Stage_I_II", "Stage_III_IV")
rownames(tnm_stage_data) <- c("Porto", "South_African")

# Perform Chi-square test
chisq.test(tnm_stage_data)

# If any expected frequencies are less than 5, you can use Fisher's Exact Test instead:
fisher.test(tnm_stage_data)

