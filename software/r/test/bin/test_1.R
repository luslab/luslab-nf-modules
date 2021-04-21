#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set run location
if(length(commandArgs(trailingOnly = TRUE)) == 0){
  cat('No command line arguments provided, user defaults paths are set for running interactively in Rstudio on docker\n')
  opt$runtype = "user"
} else {
  if(is.null(opt$runtype)){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
  if(tolower(opt$runtype) != "user" & tolower(opt$runtype) != "nextflow"){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
}

# Set paths and load data
{
  if (opt$runtype == "user"){
      
    input_path = "./results/NF-downstream_analysis_test/test_data/"
    output_path = "./output/NF-downstream_analysis_test/test_1/"
    dir.create(output_path, recursive = T)


  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    input_path = "./input/"
    output_path = "./"

  }
  
}


test_data <- read.csv(paste0(input_path, 'test.csv'), sep = '\t')[1:2,]

write.csv(test_data, paste0(output_path, 'test_out.csv'))



