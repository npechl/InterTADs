
#' ensembl_ids
#'
#' @param expr_data Parent index of expression data.
#' @param input_file
#'
#' @import biomaRt
#' @import data.table
#'
#' @description
#'
#' @return
#' A table with matched ensembl_ids
#'
#'
#' @export
#'
#' @examples
#' result_ensmbl <- ensembl_ids(
#' input_file= methylo_result,
#' expr_data = 3)
#'
#'




ensembl_ids <- function(input_file,
                        expr_data = 3){

    # data.all <- fread(paste(output_folder,
    #                         input_file,
    #                         sep = ""),sep = "\t")

    #print(input_file)
    data.all <- input_file


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


    # write.table(data.all,
    #             file = paste(output_folder,
    #                          "/integrated-tad-table-methNorm-ensembl.txt",
    #                          sep = ""),
    #             col.names = TRUE,
    #             row.names = FALSE,
    #             quote = FALSE,
    #             sep = "\t")
    if (length(ensembl.expression$chromosome_name)==0){
        message("None ensembl id was mached")
        return(NULL)
    }
    else{

        return(data.all)
    }

}




