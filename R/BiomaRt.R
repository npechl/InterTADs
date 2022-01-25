library(tidyverse)
library(biomaRt)
library(data.table)

rm(list = ls())

#output_folder <- "C:/Users/vasil/Desktop/Phd/InterTADs-data/results_bloodcancer/"
#expr_data <- 1
#tech <- "hg19"

ensembl_ids <- function(output_folder = NULL,
                        expr_data = NULL,
                        tech = NULL) {
  data.all <- fread(paste(output_folder, "integrated-tad-table-methNorm.txt",
                          sep = ""
  ),
  sep = "\t"
  )
  
  expression <- data.all[which(data.all$parent == expr_data), ]
  
  expression$names <- paste(expression$chromosome_name, expression$start_position, expression$end_position, sep = ":")
  
  ensembl.expression <- list()
  
  
  if (any(str_detect(expression$ID, "ENSG"))) {
    ensembl.expression[[1]] <- expression[which(str_detect(expression$ID, "ENSG")), ]
    
    expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]
  }
  
  if (nrow(expression) > 0) {
    if (tech == "hg19") {
      ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
      
      new_gene_ids <- getBM(
        attributes = c("chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
        mart = ensembl
      )
      
      new_gene_ids$names <- paste(new_gene_ids$chromosome_name, new_gene_ids$start_position, new_gene_ids$end_position, sep = ":")
      
      who <- match(expression$names, new_gene_ids$names)
      expression$ID <- new_gene_ids[who,]$ensembl_gene_id
      
      expression <- expression[, -c("names")]
      
      ensembl.expression[[2]] <- expression[which(str_detect(
        expression$ID,
        "ENSG"
      )), ]
      
      expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]
    } else {
      ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      
      new_gene_ids <- getBM(
        attributes = c("chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
        mart = ensembl
      )
      
      new_gene_ids$names <- paste(new_gene_ids$chromosome_name, new_gene_ids$start_position, new_gene_ids$end_position, sep = ":")
      
      who <- match(expression$names, new_gene_ids$names)
      expression$ID <- new_gene_ids[who,]$ensembl_gene_id
      
      expression <- expression[, -c("names")]
      
      ensembl.expression[[2]] <- expression[which(str_detect(
        expression$ID,
        "ENSG"
      )), ]
      
      expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]
    }
  }
  
  ensembl.expression <- rbindlist(ensembl.expression)
  
  other <- data.all[which(data.all$parent != expr_data), ]
  
  
  data.all <- rbind(expression, other, ensembl.expression)
  
  write.table(data.all,
              file = paste(output_folder,
                           "integrated-tad-table-methNorm-ensembl_1.txt",
                           sep = ""
              ),
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = "\t"
  )
  return(TRUE)
}


result_ensmbl <- ensembl_ids(
  output_folder = "C:/Users/vasil/Desktop/Phd/InterTADs-data/results_bloodcancer/",
  expr_data = 1,
  tech = "hg19"
)
print(result_ensmbl)
