# Load libraries
library(argparse)    # For command line argument parsing
library(data.table)  # For fast data manipulation
library(tidyverse)   # For data wrangling and visualization

# Define command line arguments
parser <- argparse::ArgumentParser()

parser$add_argument(
                    "--input",
                    type = "character",
                    help = "Input file path")

parser$add_argument("--output",
                    type = "character",
                    help = "Output file path")

args <- parser$parse_args(commandArgs(TRUE))

# Read in the data using data.table::fread for fast reading of large files
vcf <- data.table::fread(args$input, nThread = 4)

# Reshape the data
vcf <- vcf %>%
  mutate(across(23:ncol(.), as.character)) %>% # convert columns to be reshaped to character type
  mutate(row_id = row_number()) %>%            # add a unique row identifier
  pivot_longer(
               cols = -c(CHROM:Errors, row_id),
               names_to = c("sample_name", "variable"),
               names_sep = "\\.") %>%
  pivot_wider(names_from = "variable", values_from = "value") %>%
  select(-row_id) %>% # drop the row_id column
  mutate(
    alt_allele_freq = if_else(str_detect(AD, ",0$"), 0,
                          as.numeric(str_extract(AD, "\\d+$")) / (as.numeric(str_extract(AD, "^[^,]+")) + as.numeric(str_extract(AD, "\\d+$")))),
    alt_allele_freq = round(alt_allele_freq, 3)
  ) %>%
  filter(! alt_allele_freq == 0) %>%
  arrange(sample_name, POS)

data.table::fwrite(vcf, args$output, nThread = 4)

write_tsv(vcf, args$output)