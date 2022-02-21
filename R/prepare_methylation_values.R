# Loading libraries ------------------------------------------------------------

# rm(list = ls())

# source("R/libraries.R")

# Inputs -----------------------------------------------------------------------

#' Input parameters for TADiff part
#'
#' @param dir_name Directory of input datasets containing feature counts
#' and frequency tables
#'
#' @param output_folder Folder name for printing output tables
#'
#' @param meta meta-data file name used
#'
#' @param names.meta meta data columns to process (names or indexes)
#'
#' @param meth_data Parent index of methylation data.
#' If no methylation is provided, place FALSE
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
#' prepare_methylation_values(
#' dir_name = system.file("extdata","Datasets",package='InterTADs'),
#' output_folder = system.file("extdata","results_bloodcancer",
#'                             package='InterTADs'),
#' meta = "meta-data.csv",
#' meth_data = 1)



prepare_methylation_values <- function (
    integratedTADtable, 
    mapping_file
) {

    data.all <- fread(integratedTADtable)

    meta <- fread(mapping_file)


    sample.list <- meta$newNames


    meth <- data.all[which(data.all$parent == meth_data), ]


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
