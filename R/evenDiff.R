# Loading libraries ------------------------------------------------------------


# source("R/libraries.R")
#
# start_tad_time <- Sys.time()

#' evenDiff
#'
#' @param names.meta meta data columns to process (names or indexes)
#' @param mapping_file A meta-data file
#' @param methylo_result
#'
#' @import data.table
#' @import limma
#'
#' @description
#'
#' @return
#'
#' @export
#'
#' @examples
#' result_evenDiff <- evenDiff(
#' mapping_file = system.file("extdata", "Datasets",
#' "meta-data.csv", package="InterTADs"),
#' methylo_result = methylo_result,
#' names.meta = c('group'))
#'
#' print(result_evenDiff)
#'
#'


evenDiff <- function(
                    mapping_file = NULL,
                    methylo_result,
                    names.meta = NULL){

    data.all <- methylo_result

    data.all$ID <- paste(data.all$tad_name, data.all$ID, sep = ";")

    meta <- fread(mapping_file)
    who <- meta == ""
    who <- apply(who, 1, sum, na.rm = TRUE)
    meta <- meta[which(who == 0), ]


    sample.list <- meta[[1]]

    who <- c("ID", "Gene_id", "parent")

    gene.parents <- data.all[,..who]

    colnames(gene.parents) <- c("ID", "gene_ID", "parents")


    diffevent <- matrix(data = "0",
                        nrow = length(names.meta),
                        ncol = (length(unique(data.all$parent)) + 1))

    colnames(diffevent) <- c("Analysis", paste("p_", unique(data.all$parent),
                                               sep = ""))

    rm(who)
    Diff_list <- list()
    for (z in 1:length(names.meta)){

        #print(names.meta[["IGHV"]])
        cat(c(names.meta[z], "\n"))

        analysis <- names.meta[z]
        groups <- as.character(meta[[analysis]])
        groups <- unique(groups)
        groups <- groups[which(groups != "")]
        groups <- groups[!is.na(groups)]

        group1 <- groups[1]
        group2 <- groups[2]


        meta.keep <- meta[which(meta[[analysis]] == group1 | meta[[analysis]] ==
                                  group2), ]


        sample.list <- meta.keep[[1]]

        df <- data.all[,..sample.list]

        df <- as.data.frame(df)
        row.names(df) <- data.all$ID

        pheno <- as.factor(meta.keep[[analysis]])

        phenoMat <- model.matrix(~pheno)
        colnames(phenoMat) <- sub("^pheno", "", colnames(phenoMat))

        fit <- limma::lmFit(object = df, design = phenoMat)

        gc()
        set.seed(6)
        fit <- limma::eBayes(fit)

        gc()
        top.rank <- limma::topTable(fit, number = nrow(df), adjust.method = "fdr",
                             sort.by = "p")

        sign.table <- top.rank[which(top.rank$adj.P.Val <= 0.01 &
                                       abs(top.rank$logFC) > 4), ]

        if (nrow(sign.table) == 0) {

            cat(c("No statistical significant events for:",
                    names.meta[z], "\n"))

            diffevent[z,1] <- analysis }

        else {

            # annotate sign.table

            sign.table$ID <- row.names(sign.table)

            sign.genes <- merge(gene.parents,
                                sign.table,
                                by.x = "ID",
                                by.y = "ID")


            sign.table <- merge(sign.table,
                                data.all,
                                by.x = "ID",
                                by.y = "ID")

            sign <- sign.genes %>% dplyr::count(parents)


            diffevent[z, 1] <- analysis
            diffevent[z, paste("p_", sign$parents, sep = "")] <- sign$n

            Diff_list[[analysis]] = sign.table

        }
    }

    diffevent <- as.data.frame(diffevent)

    for(i in 2:ncol(diffevent)){
        diffevent[,i] <- as.numeric(diffevent[,i])
        #print(diffevent[,i])
    }


    diffevent$`total events` <- rowSums(diffevent[,2:ncol(diffevent)]) #rowSums



    return(list(Diff_list,diffevent))

}


# result_evenDiff <- evenDiff(
# mapping_file = "/Users/aspaor/Downloads/bloodcancer/metaData_groups.csv",
# methylo_result = methylo_result,
# names.meta = c("IGHV","Gender"))
# mapping_file = "/Users/aspaor/Downloads/bloodcancer/metaData_groups.csv"
# methylo_result = methylo_result
# names.meta = c("IGHV","Gender")
# #
# print(result_evenDiff)
