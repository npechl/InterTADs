#' ensembl_ids
#'
#' @param expr_data 
#' Parent index of expression data.
#' 
#' @param integratedTADtable
#' IntegratedTADtable contains all data info about TADs. 
#' This is the output integrated table 
#' from data_integration function (or prepare_methylation_values function).
#'
#' @import biomaRt
#' @import data.table
#'
#' @description
#' Convert external gene IDs to ensembl IDs. 
#' In case the user has no gene IDs or do not wish to convert 
#' the provided gene IDs to ensembl IDs, is able to skip this function.
#' Note the the functional analysis will perform only with ensembl IDs.
#'
#' @return
#' A table with matched ensembl_ids.
#'
#'
#' @export
#'
#' @examples
#' result<- data_integration (
#' 
#'     counts_folder = system.file(
#'         "extdata", "Datasets", "counts", package = "InterTADs"
#'     ),
#' 
#'     freq_folder = system.file(
#'         "extdata", "Datasets", "freq", package = "InterTADs"
#'     ),
#' 
#'     mapping_file = system.file(
#'         "extdata", "Datasets", "meta-data.csv", package = "InterTADs"
#'     ),
#' 
#'     tad_file =system.file(
#'         "extdata", "Datasets",
#'         "hglft_genome_2dab_ec1330.bed", package = "InterTADs"
#'     ),
#' 
#'     tech = "hg19"
#' )
#' 
#' result_ensmbl <- ensembl_ids(
#'      integratedTADtable = result[[1]],
#'     expr_data = 3
#' )


ensembl_ids <- function(
    integratedTADtable = NULL,
    expr_data = NULL
) {

    if ((is.null(integratedTADtable) ) || ((is.null(expr_data)) )) {
      stop("Please provide all the input parameters.")
    }
  
    expression <- integratedTADtable[which(integratedTADtable$parent == expr_data), ]

    ensembl.expression <- list()


    if( any(str_detect(expression$ID, "ENSG")) ) {

        ensembl.expression[[1]] <- expression[which(
            str_detect(expression$ID, "ENSG")
        ), ]
        
        expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]

    }


    if(nrow(expression) > 0) {

        ensembl <- useEnsembl(
            biomart = "genes",
            dataset = "hsapiens_gene_ensembl"
        )

        new_gene_ids <- getBM(
            attributes = c('entrezgene_id','ensembl_gene_id'),
            filters = 'entrezgene_id',
            values = expression$ID,
            mart = ensembl
        )

        who <- match(new_gene_ids$entrezgene_id, expression$ID)
        
        expression[who,]$ID <- new_gene_ids$ensembl_gene_id

        ensembl.expression[[2]] <- expression[which(
            str_detect(expression$ID, "ENSG")
        ), ]
        
        expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]

    }



    if(nrow(expression) > 0) {

        new_gene_ids <- getBM(
            attributes = c('external_gene_name', 'ensembl_gene_id'),
            filters = 'external_gene_name',
            values = expression$ID,
            mart = ensembl
        )

        who <- match(new_gene_ids$external_gene_name, expression$ID)
        expression[who,]$ID <- new_gene_ids$ensembl_gene_id

        ensembl.expression[[3]] <- expression[which(
            str_detect(expression$ID, "ENSG")
        ), ]
        
        expression <- expression[which(!str_detect(expression$ID, "ENSG")), ]

    }




    ensembl.expression <- rbindlist(ensembl.expression)


    other <- integratedTADtable[which(integratedTADtable$parent != expr_data), ]


    integratedTADtable <- rbind(expression, other)


    if (length(ensembl.expression$chromosome_name) == 0){
        
        message("None ensembl IDs were matched")
        return(integratedTADtable)
        
    } else {

        return(integratedTADtable)
    }

}


