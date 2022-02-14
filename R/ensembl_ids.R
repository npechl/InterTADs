# Loading libraries ------------------------------------------------------------


#
# source("R/libraries.R")
#
# library(biomaRt)


#' ensembl_ids
#'
#' @param dir_name Directory of input datasets containing feature counts and
#' frequency tables
#'
#' @param output_folder Folder name for printing output tables
#'
#' @param expr_data Parent index of expression data.
#' If no expression is provided, place FALSE
#' @import biomaRt
#' @import data.table
#' @import systemPipeR
#' @import data.table
#' @import tidyverse
#' @import org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import annotables
#' @import GenomicRanges
#' @import gplots
#' @import gghalves
#' @import limma
#' @importFrom  utils data  write.table
#'
#' @description
#'
#' @return
#' TRUE if function runs properly
#'
#'
#' @export
#'
#' @examples
#' result_ensmbl <- ensembl_ids(dir_name = "Datasets",
#' output_folder = 'results_bloodcancer',
#' expr_data = 2)
#' print(result_ensmbl)
#'




ensembl_ids <- function(dir_name = NULL,
                        output_folder = NULL,
                        expr_data = NULL){

    data.all <- fread(paste(output_folder,
                            "/integrated-tad-table-methNorm.txt",
                            sep = ""),
    sep = "\t")

    expression <- data.all[which(data.all$parent == expr_data), ]

    ensembl.expression <- list()


    if( any(str_detect(expression$ID, "ENSG")) ) {

        ensembl.expression[[1]] <- expression[which(str_detect(expression$ID,
                                                               "ENSG")), ]
        expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]

    }


    if(nrow(expression) > 0) {

        ensembl <- useEnsembl(biomart = "genes",
                            dataset = "hsapiens_gene_ensembl")

        new_gene_ids <- getBM(attributes = c('entrezgene_id','ensembl_gene_id'),
                                filters = 'entrezgene_id',
                                values = expression$ID,
                                mart = ensembl)

        who <- match(new_gene_ids$entrezgene_id, expression$ID)
        expression[who,]$ID <- new_gene_ids$ensembl_gene_id

        ensembl.expression[[2]] <- expression[which(str_detect(expression$ID,
                                                               "ENSG")), ]
        expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]

    }



    if(nrow(expression) > 0) {

        new_gene_ids <- getBM(attributes = c('external_gene_name',
                                             'ensembl_gene_id'),
                                filters = 'external_gene_name',
                                values = expression$ID,
                                mart = ensembl)

        who <- match(new_gene_ids$external_gene_name, expression$ID)
        expression[who,]$ID <- new_gene_ids$ensembl_gene_id

        ensembl.expression[[3]] <- expression[which(str_detect(expression$ID,
                                                               "ENSG")), ]
        expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]

    }




    ensembl.expression <- rbindlist(ensembl.expression)


    other <- data.all[which(data.all$parent != expr_data), ]


    data.all <- rbind(expression, other)

    write.table(data.all,
                file = paste(output_folder,
                             "/integrated-tad-table-methNorm-ensembl.txt",
                             sep = ""),
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
    if (length(ensembl.expression$chromosome_name)==0){
        message("None ensembl id was mached")
        return(FALSE)
    }
    else{
        return(TRUE)
    }

}

# result_ensmbl <- ensembl_ids(dir_name = "Datasets",
#                              output_folder = 'results_bloodcancer',
#                              expr_data = 2)


