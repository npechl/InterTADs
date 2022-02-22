
#' prepare_methylation_values
#'
#' @param meth_data Parent index of methylation data, from meta-data.csv.
#' If no methylation is provided, place FALSE
#' @param integratedTADtable IntegratedTADtable contains
#' all data info about TADs
#'
#' @param mapping_file A meta-data file
#'
#' @import data.table
#'
#' @description
#'
#' @return
#'
#' @export
#'
#' @examples
#' methylo_result <- prepare_methylation_values(
#' integratedTADtable = result[[1]],
#' mapping_file = system.file("extdata", "Datasets",
#'                          "meta-data.csv", package="InterTADs"),
#' meth_data = 2
#' )



prepare_methylation_values <- function (
    integratedTADtable,
    mapping_file,
    meth_data = 2
) {

    data.all <- integratedTADtable

    meta <- fread(mapping_file)


    sample.list <- meta$newNames


    meth <- data.all[which(data.all$parent == meth_data)]

    meth.promoter <- meth[which(str_detect(meth$Gene_locus, "promoter")) ]

    meth.intergenic <- meth[which(str_detect(meth$Gene_locus, "intergenic")) ]

    meth.regulatory <- rbind(meth.promoter, meth.intergenic)
    meth.regulatory <- meth.regulatory[!duplicated(meth.regulatory$ID) ]

    for(i in sample.list){
        meth.regulatory[[i]] <- 100 - meth.regulatory[[i]]
    }

    other <- data.all[which(data.all$parent != meth_data) ]

    data.all <- rbind(meth.regulatory, other)

    return(data.all)
}





