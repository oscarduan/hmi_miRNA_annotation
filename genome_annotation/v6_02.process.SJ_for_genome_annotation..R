setwd("C:/Users/duany/OneDrive/Desktop/Srivastava Lab/Results/PacBio/")
library(dplyr)
library(readr)
library(stringr)
library(data.table)

convert.SJ <- function(file_path){
  base <- str_replace(basename(file_path), ".SJ.out.tab", "")
  
  # Read the SJ.out.tab file
  data <- fread(file_path, select = c(1, 2, 3, 4), header = FALSE, col.names = c("Chr", "Start", "End", "Strand")) # nolint: line_length_linter.
  
  # Convert the numeric strand information to +, -, or .
  data[, Strand := fifelse(Strand == 1, "+", fifelse(Strand == 2, "-", "."))]
  
  # Write the output file
  output_path <- paste0("v6_annotation/first_STAR_SJ/", base, ".SJ1.txt")
  fwrite(data, file = output_path, sep = "\t", col.names = FALSE, quote = FALSE)
}

# Get a list of all SJ.out.tab files in the directory
input_dir <- "v6_annotation/first_STAR_SJ/"
files <- list.files(input_dir, pattern = "\\.SJ.out.tab$", full.names = TRUE)

# Process each file
lapply(files, convert.SJ)
