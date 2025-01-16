
library(data.table)

# Load the CSV data
csv_file <- "vcf.csv"
vcf_data <- fread(csv_file)

# Check the structure of the CSV to ensure necessary columns exist
print(names(vcf_data))

# Function to identify somatic mutations
extract_somatic <- function(vcf_data) {
  # Parse the INFO field for 'SOMATIC' flag
  somatic_flag <- grepl("SOMATIC", vcf_data$INFO)
  
  # Split FORMAT and extract GT (Genotype) field for NORMAL and TUMOR
  format_fields <- strsplit(vcf_data$FORMAT, ":")
  normal_fields <- strsplit(vcf_data$NORMAL, ":")
  tumor_fields <- strsplit(vcf_data$TUMOR, ":")
  
  # Extract GT values
  normal_gt <- sapply(normal_fields, function(x) x[which(format_fields[[1]] == "GT")])
  tumor_gt <- sapply(tumor_fields, function(x) x[which(format_fields[[1]] == "GT")])
  
  # Filter for somatic mutations: SOMATIC flag present and NORMAL = 0/0, TUMOR = 1/1 or 0/1
  somatic_mutations <- vcf_data[somatic_flag & normal_gt == "0/0" & tumor_gt %in% c("1/1", "0/1")]
  
  return(somatic_mutations)
}

# Apply the function
somatic_mutations <- extract_somatic(vcf_data)

# Output the somatic mutations
print(somatic_mutations)

# save the results to a file
fwrite(somatic_mutations, "somatic_mutations_from_csv.tsv", sep = "\t")
