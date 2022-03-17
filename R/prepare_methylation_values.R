
#' prepare_methylation_values
#'
#' @param meth_data 
#' Parent index of methylation data, from meta-data.csv.
#' 
#' @param integratedTADtable 
#' IntegratedTADtable contains all data info about TADs. This is the output
#' integrated table from data_integration function.
#'
#' @param mapping_file 
#' A mapping file containing mapping betwen the columns of the input 
#' NGS datasets. The first column corresponds to the sample ID that 
#' will be used in the output file. The following columns correspond 
#' to the column names of the input datasets. 
#' The file also contains the related sample meta-data
#'
#' @import data.table
#'
#' @description
#' 
#' In case DNA methylation datasets are provided this function reverses 
#' methylation values (0% -> 100%, 100% -> 0%). The user can skip this 
#' function in case it is not suitable for the provided dataset. 
#'
#' @return
#'
#' An integrated TAD data table with the reversed methylation 
#' values is returned.
#'
#' @export
#'
#' @examples
#'
#' result <- data_integration (
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
#'        "extdata", "Datasets", "meta-data.csv", package = "InterTADs"
#'    ),
#'
#'     tad_file =system.file(
#'        "extdata", "Datasets",
#'        "hglft_genome_2dab_ec1330.bed", package = "InterTADs"
#'     ),
#' 
#'     tech = "hg19"
#' )
#' 
#' 
#' methylo_result <- prepare_methylation_values(
#'     integratedTADtable = result[[1]],
#'     mapping_file = system.file(
#'         "extdata", "Datasets", "meta-data.csv", package = "InterTADs"
#'     ),
#'     meth_data = 2
#' )



prepare_methylation_values <- function (
    integratedTADtable,
    mapping_file,
    meth_data
) {

    data.all <- integratedTADtable

    meta <- fread(mapping_file)


    sample.list <- meta[[1]]



    meth <- data.all[which(data.all$parent == meth_data),]

    meth.promoter <- meth[which(str_detect(meth$Gene_locus, "promoter")), ]

    meth.intergenic <- meth[which(str_detect(meth$Gene_locus, "intergenic")), ]

    meth.regulatory <- rbind(meth.promoter, meth.intergenic)
    meth.regulatory <- meth.regulatory[!duplicated(meth.regulatory$ID), ]

    for(i in sample.list){
        meth.regulatory[[i]] <- 100 - meth.regulatory[[i]]
    }

    other <- data.all[which(data.all$parent != meth_data), ]

    data.all <- rbind(meth.regulatory, other)

    return(data.all)
}





