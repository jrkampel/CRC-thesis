# Load Libs
library(readxl)
library(dplyr)
library(readr)
install.packages("writexl")
library(writexl)

# Load data 
likely_sporadic_data <- readxl::read_excel("C:/Users/jrkam/OneDrive - University of Cape Town/Honours/thesis/Finding genes/likley_sporadic_data.xlsx")
likely_ls_data <- readxl::read_excel("C:/Users/jrkam/OneDrive - University of Cape Town/Honours/thesis/Finding genes/likley_ls_data.xlsx")
allmuts2_CRC <- read_delim("C:/Users/jrkam/OneDrive - University of Cape Town/Honours/thesis/Finding genes/allmuts2_crc_new.csv", 
                           delim = ",", escape_double = FALSE, trim_ws = TRUE)
mutation_data <- read.delim("C:/Users/jrkam/OneDrive - University of Cape Town/Honours/thesis/Finding genes/data_mutations.txt", header = TRUE, sep = "\t")


# create new column combining all these rows 
mutation_data <- mutation_data %>%
  mutate(newID = paste(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = "_"))

mutation_data$newID
allmuts2_CRC$newID

# counts how many tumor samples in mutation_data are also listed in likely_ls_data and num_likely_sporadic
num_likely_LS <- length(intersect(unique(mutation_data$Tumor_Sample_Barcode), likely_ls_data$Patient_ID_LS))
num_likely_sporadic <- length(intersect(unique(mutation_data$Tumor_Sample_Barcode), likely_sporadic_data$Patient_ID_sporadic))

# Filter out pathogenic 
mutation_data <- mutation_data %>%
  filter(newID %in% allmuts2_CRC$newID)


# Filter out according to excels
mutation_data_likely_LS <- mutation_data %>% filter(Tumor_Sample_Barcode %in% likely_ls_data$Patient_ID_LS)
mutation_data_likely_sporadic <- mutation_data %>% filter(Tumor_Sample_Barcode %in% likely_sporadic_data$Patient_ID_sporadic)

length(unique(mutation_data_likely_sporadic$Tumor_Sample_Barcode))

# Extract sample IDs from the Excel spreadsheets
sporadic_samples <- likely_sporadic_data$Patient_ID_sporadic
lynch_samples <- likely_ls_data$Patient_ID_LS

# Count the genes in each sample
gene_count_ls <- mutation_data_likely_LS %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Unique_Genes = length(unique(Hugo_Symbol)))

# Count the genes in each sample
gene_count_sporadic <- mutation_data_likely_sporadic %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Unique_Genes = length(unique(Hugo_Symbol)))

# List of all genes present in the data
all_genes <- unique(c(mutation_data_likely_LS$Hugo_Symbol, mutation_data_likely_sporadic$Hugo_Symbol))

# Create a function to generate contingency tables for each gene
generate_contingency_table <- function(gene, lynch_data, sporadic_data) {
  # Subset the data for the specific gene
  lynch_gene_data <- lynch_data[lynch_data$Hugo_Symbol == gene, "Tumor_Sample_Barcode", drop = FALSE]
  sporadic_gene_data <- sporadic_data[sporadic_data$Hugo_Symbol == gene, "Tumor_Sample_Barcode", drop = FALSE]
  
  # Count the number of mutated and not mutated samples for each dataset
  mutated_lynch <- length(unique(lynch_gene_data$Tumor_Sample_Barcode))
  not_mutated_lynch <- num_likely_LS - mutated_lynch
  
  mutated_sporadic <- length(unique(sporadic_gene_data$Tumor_Sample_Barcode))
  not_mutated_sporadic <- num_likely_sporadic - mutated_sporadic
  
  # Create the contingency table
  contingency_table <- matrix(c(mutated_lynch, not_mutated_lynch, mutated_sporadic, not_mutated_sporadic), nrow = 2, byrow = TRUE)
  colnames(contingency_table) <- c("Mutated", "Not Mutated")
  rownames(contingency_table) <- c("Likely Lynch", "Likely Sporadic")
  
  return(contingency_table)
}

# Create contingency tables for each gene
contingency_tables <- lapply(all_genes, function(gene) {
  generate_contingency_table(gene, mutation_data_likely_LS, mutation_data_likely_sporadic)
})

# Show contingency tables for each gene
names(contingency_tables) <- all_genes


# Perform Fisher's exact test for each gene
fisher_results <- lapply(all_genes, function(gene) {
  # Access the specific contingency table for the gene
  gene_table <- contingency_tables[[gene]]
  
  # Perform Fisher's exact test on the specific contingency table
  fisher_test_result <- fisher.test(gene_table)
  
  # Return the result along with the gene name
  return(list(Gene = gene, Fisher_Test_Result = fisher_test_result))
})


# Display the results
fisher_results


# Extract gene names and p-values from fisher_results
gene_names <- sapply(fisher_results, function(x) x$Gene)
p_values <- sapply(fisher_results, function(x) x$Fisher_Test_Result$p.value)

# Create a data frame with gene names and p-values
results_table <- data.frame(Gene = gene_names, P_Value = p_values)

# Show the table
results_table

# GOOD UP UNTIL HERE

# To include number of LS and LL
# Get the number of samples with a mutation in each gene ####
gene_mutation_counts_likely_LS <- mutation_data_likely_LS %>%
  group_by(Hugo_Symbol) %>%                  # Group data by gene symbol
  summarise(Unique_Samples = n_distinct(Tumor_Sample_Barcode)) %>%  # Count unique sample IDs
  arrange(desc(Unique_Samples))              # Arrange the data in descending order of counts

gene_mutation_counts_likely_sporadic <- mutation_data_likely_sporadic %>%
  group_by(Hugo_Symbol) %>%                  # Group data by gene symbol
  summarise(Unique_Samples = n_distinct(Tumor_Sample_Barcode)) %>%  # Count unique sample IDs
  arrange(desc(Unique_Samples))              # Arrange the data in descending order of counts

# Merge the data frames on Hugo_Symbol
merged_mutation_counts <- merge(
  gene_mutation_counts_likely_LS, 
  gene_mutation_counts_likely_sporadic, 
  by = "Hugo_Symbol", 
  suffixes = c("_likely_LS", "_likely_Sporadic"), 
  all = TRUE  # Includes all records from both data frames
)

# Replace NA values with 0 in the merged data frame
merged_mutation_counts[is.na(merged_mutation_counts)] <- 0

names(merged_mutation_counts)[names(merged_mutation_counts) == "Hugo_Symbol"] <- "Gene"

# Merging data frames
results_table <- merge(results_table,
                       merged_mutation_counts[, c("Gene", "Unique_Samples_likely_LS", "Unique_Samples_likely_Sporadic")],
                       by = "Gene",
                       all.x = TRUE)

####### GOOD!!! DO NOT TOUCH 

# To add frequencies in tables

# Calculate the total number of unique samples for each category
total_samples_likely_LS <- n_distinct(mutation_data_likely_LS$Tumor_Sample_Barcode)
total_samples_likely_sporadic <- n_distinct(mutation_data_likely_sporadic$Tumor_Sample_Barcode)

# Get the number of samples with a mutation in each gene for likely LS
gene_mutation_counts_likely_LS <- mutation_data_likely_LS %>%
  group_by(Hugo_Symbol) %>%
  summarise(Unique_Samples = n_distinct(Tumor_Sample_Barcode)) %>%
  arrange(desc(Unique_Samples))

# Get the number of samples with a mutation in each gene for likely sporadic
gene_mutation_counts_likely_sporadic <- mutation_data_likely_sporadic %>%
  group_by(Hugo_Symbol) %>%
  summarise(Unique_Samples = n_distinct(Tumor_Sample_Barcode)) %>%
  arrange(desc(Unique_Samples))

# Merge the data frames on Hugo_Symbol
merged_mutation_counts <- merge(
  gene_mutation_counts_likely_LS, 
  gene_mutation_counts_likely_sporadic, 
  by = "Hugo_Symbol", 
  suffixes = c("_likely_LS", "_likely_Sporadic"), 
  all = TRUE
)

# Calculate the frequencies
merged_mutation_counts <- merged_mutation_counts %>%
  mutate(Frequency_Likely_LS = Unique_Samples_likely_LS / total_samples_likely_LS,
         Frequency_Likely_Sporadic = Unique_Samples_likely_Sporadic / total_samples_likely_sporadic)

# Replace NA with 0 for genes not present in one of the categories
merged_mutation_counts[is.na(merged_mutation_counts)] <- 0

# Ensure the column names are correct before merging
print(colnames(merged_mutation_counts))
print(colnames(results_table))

# Merge the frequencies into the results_table
results_table <- merge(results_table,
                       merged_mutation_counts[, c("Hugo_Symbol", "Frequency_Likely_LS", "Frequency_Likely_Sporadic")],
                       by.x = "Gene", by.y = "Hugo_Symbol",
                       all.x = TRUE)

# Replace NA with 0 for frequencies not present in the merged results_table
results_table[is.na(results_table)] <- 0

# Display the updated results_table
print(results_table)


# # # # # # GOOOOOD!!!!!

# Calculate adjusted p-values using Benjamini-Hochberg procedure
adjusted_p_values <- p.adjust(results_table$P_Value, method = "BH")

results_table$adjusted_p_value <- adjusted_p_values


# Filter results to include only significant genes after adjustment

significant_genes <- subset(results_table, adjusted_p_value < 0.05)

# View the filtered results
print(significant_genes)

# WOKING DO NOT CHANGE YAY

# To include which group the gene is enriched in we need log
# Calculate log2-based ratio for the p-values
# Ensure the result_table_with_frequencies is defined and contains the necessary frequency columns
# Calculating log2-based ratio for the frequencies
results_table$log2_ratio <- log2(results_table$Frequency_Likely_LS / results_table$Frequency_Likely_Sporadic)

# Filter the results_table to include only significant genes
significant_genes <- results_table %>%
  filter(adjusted_p_value < 0.05)

# Adding the Enriched_in column
significant_genes <- significant_genes %>%
  mutate(Enriched_in = case_when(
    log2_ratio > 0 ~ "Likely Lynch",
    log2_ratio <= 0 ~ "Likely Sporadic",
    TRUE ~ ""
  ))

# Filter out genes with less than 4 patients
filtered_genes <- significant_genes[significant_genes$Unique_Samples_likely_LS >= 4, ]

# Read in oncogenes & tumor suppressor
cancer_gene_list <- read.delim("C:/Users/jrkam/OneDrive - University of Cape Town/Honours/thesis/Finding genes/cancerGeneList.tsv")

# check their variables
unique(cancer_gene_list$Is.Oncogene)
unique(cancer_gene_list$Is.Tumor.Suppressor.Gene)

# Filter out genes that are not oncogenes and also tumor supressor genes 
cancer_gene_list_filtered <- cancer_gene_list[cancer_gene_list$Is.Oncogene == "Yes" | cancer_gene_list$Is.Tumor.Suppressor.Gene == "Yes", ]

final_filtered_data <- filtered_genes %>%
  filter(Gene %in% cancer_gene_list_filtered$Hugo.Symbol)

# Export final_filtered_genes data frame to a CSV file
write.csv(final_filtered_data, "final_filtered_data.csv", row.names = FALSE)

final_genes <- final_filtered_data$Gene

# Upload KEGG 2021 Human list from erichr
install.packages("enrichR")
library(enrichR)
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs <- "KEGG_2021_Human"
enriched <- enrichr(final_genes, dbs)

# Make a table for enrich data
df <- as.data.frame(enriched$KEGG_2021_Human)

# Plot the enrich data
plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
