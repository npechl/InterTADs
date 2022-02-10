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
#' @description
#'
#' @return
#'
#' @export
#'
#' @examples



prepare_methylation_values <- function (dir_name = NULL,
                                        output_folder = NULL,
                                        meta = NULL,
                                        meth_data = NULL){

    data.all <- fread(paste(output_folder, "/integrated-tad-table.csv",
                        sep = ""),
                        sep = "\t")
    meta <- fread(paste(dir_name, meta, sep = "/"))
    who <- meta == ""
    who <- apply(who, 1, sum, na.rm = TRUE)
    meta <- meta[which(who == 0), ]


    sample.list <- meta$newNames


    meth <- data.all[which(data.all$parent == meth_data), ]


    meth.promoter <- meth[which(str_detect(meth$Gene_locus, "promoter")), ]

    meth.intergenic <- meth[which(str_detect(meth$Gene_locus, "intergenic")), ]

    meth.regulatory <- rbind(meth.promoter, meth.intergenic)
    meth.regulatory <- meth.regulatory[!duplicated(meth.regulatory$ID), ]

    rm(meth.intergenic, meth.promoter, meth)

    for(i in sample.list){
        meth.regulatory[[i]] <- 100 - meth.regulatory[[i]]
    }

    other <- data.all[which(data.all$parent != meth_data), ]

    data.all <- rbind(meth.regulatory, other)

    write.table(data.all,
                file = paste(output_folder,
                "/integrated-tad-table-methNorm.txt",
                sep = ""),
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")


 return(NULL)
}


# prepare_methylation_values(dir_name = "Datasets",
#                            output_folder = "results_bloodcancer",
#                            meta = "meta-data.csv",
#                            meth_data = 1)
