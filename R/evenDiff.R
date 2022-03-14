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
#' result<- data_integration (
#' counts_folder = system.file("extdata", "Datasets",
#'                          "counts", package="InterTADs"),
#' counts_fls = NULL,
#' freq_folder = system.file("extdata", "Datasets",
#'                          "freq", package="InterTADs"),
#' freq_fls = NULL,
#' mapping_file = system.file("extdata", "Datasets",
#'                          "meta-data.csv", package="InterTADs"),
#'
#' tad_file =system.file("extdata", "Datasets",
#'                      "hglft_genome_2dab_ec1330.bed", package="InterTADs"),
#' tech = "hg38"
#' )
#'
#' methylo_result <- prepare_methylation_values(
#' integratedTADtable = result[[1]],
#' mapping_file = system.file("extdata", "Datasets",
#'                          "meta-data.csv", package="InterTADs"),
#' meth_data = 2
#' )
#'
#' result_evenDiff <- evenDiff(
#' mapping_file = system.file("extdata", "Datasets",
#' "meta-data.csv", package="InterTADs"),
#' methylo_result = methylo_result,
#' names.meta = c('group'),
#' adj.PVal = 0.01,
#' log_thr = 4)
#'
#' print(result_evenDiff)
#'
#'


evenDiff <- function(
                    mapping_file = NULL,
                    methylo_result,
                    names.meta = NULL,
                    adj.PVal = NULL,
                    log_thr = NULL){

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
        top.rank <- limma::topTable(fit, number = nrow(df),
                            adjust.method = "fdr",sort.by = "p")

        sign.table <- top.rank[which(top.rank$adj.P.Val <= adj.PVal &
                                       abs(top.rank$logFC) > log_thr), ]

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


