# Loading libraries -------------------------------------------------------

rm(list = ls())

source("R/libraries.R")

library(biomaRt)

# Inputs ------------------------------------------------------------------------------

#' Input parameters for TADiff part
#' 
#' @param dir_name Directory of input datasets containing feature counts and frequency tables
#' 
#' @param output_folder Folder name for printing output tables
#' 
#' @param meta meta-data file name used
#' 
#' @param names.meta meta data columns to process (names or indexes)
#' 
#' @param meth_data Parent index of methylation data. If no methylation is provided, place FALSE

dir_name <- "Datasets"

output_folder <- "results_bloodcancer"

# meta = "metaData_groups.csv"

expr_data <- 2

data.all <- fread(paste(output_folder, "/integrated-tad-table-methNorm.txt", sep = ""),
                 sep = "\t")

expression <- data.all[which(data.all$parent == expr_data), ]

ensembl.expression <- list()

if( any(str_detect(expression$ID, "ENSG")) ) {
  
  ensembl.expression[[1]] <- expression[which(str_detect(expression$ID, "ENSG")), ]
  expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]
  
}


if(nrow(expression) > 0) {
  
  ensembl <- useEnsembl(biomart = "genes", 
                       dataset = "hsapiens_gene_ensembl")
  
  new_gene_ids <- getBM(attributes = c('entrezgene_id', 'ensembl_gene_id'),
                       filters = 'entrezgene_id',
                       values = expression$ID, 
                       mart = ensembl)
  
  who <- match(new_gene_ids$entrezgene_id, expression$ID)
  expression[who,]$ID <- new_gene_ids$ensembl_gene_id
  
  ensembl.expression[[2]] <- expression[which(str_detect(expression$ID, "ENSG")), ]
  expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]
  
}



if(nrow(expression) > 0) {
  
  new_gene_ids <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'),
                       filters = 'external_gene_name',
                       values = expression$ID, 
                       mart = ensembl)
  
  who <- match(new_gene_ids$external_gene_name, expression$ID)
  expression[who,]$ID <- new_gene_ids$ensembl_gene_id
  
  ensembl.expression[[3]] <- expression[which(str_detect(expression$ID, "ENSG")), ]
  expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]
  
}




ensembl.expression <- rbindlist(ensembl.expression)
  

other <- data.all[which(data.all$parent != expr_data), ]

data.all <- rbind(expression, other)

write.table(data.all, 
            file = paste(output_folder, "/integrated-tad-table-methNorm-ensembl.txt", sep = ""), 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")
