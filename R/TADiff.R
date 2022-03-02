
#' TADiff
#'
#' @param names.meta meta data columns to process (names or indexes)
#'
#' @param expr_data Parent index of expression data.
#' If no expression is provided, place FALSE
#' @param mapping_file A meta-data file
#' @param methylo_result
#'
#' @import data.table
#' @import dplyr
#' @description
#'
#' @return
#'
#' @export
#'
#' @examples
#' TADiff_result <- TADiff(
#'
#' mapping_file = system.file("extdata", "Datasets",
#'                      "meta-data.csv", package="InterTADs"),
#'                        methylo_result = methylo_result,
#'                        names.meta = c("group"),
#'                    expr_data = 3)
#' TADiff_result


TADiff<- function(
                  mapping_file = NULL,
                  methylo_result,
                  names.meta = NULL,
                  expr_data = NULL){


    data.all <- methylo_result

    data.all$ID <- paste(data.all$tad_name, data.all$ID, sep = ";")


    meta <- fread(mapping_file)

    who <- meta == ""
    who <- apply(who, 1, sum, na.rm = TRUE)
    meta <- meta[which(who == 0), ]


    sample.list <- meta[[1]]
    data.all$Mean <- rowMeans(data.all[,..sample.list])


    data.all <- data.all[which(round(data.all$Mean) > 10), ]


    sign_table <- matrix(0, nrow = length(names.meta), ncol = 2)

    rm(who)
    Diff_list <- list()
    for (z in 1:length(names.meta)){

        cat(c(names.meta[z], "\n"))

        analysis <- names.meta[z]
        groups <- as.character(meta[[analysis]])
        groups <- unique(groups)
        groups <- groups[which(groups != "")]
        groups <- groups[!is.na(groups)]

        group1 <- groups[1]
        group2 <- groups[2]

        meta.keep <- meta[which(meta[[analysis]] == group1 |
                                    meta[[analysis]] ==group2), ]

        sample.list <- meta.keep[[1]]

        df <- data.all[,..sample.list]

        df <- as.data.frame(df)
        row.names(df) <- data.all$ID

        pheno <- as.factor(meta.keep[[analysis]])


        phenoMat <- model.matrix(~pheno)
        colnames(phenoMat) <- sub("^pheno", "", colnames(phenoMat))

        fit <- lmFit(object = df, design = phenoMat)

        gc()
        set.seed(6)
        fit <- eBayes(fit)

        gc()
        top.rank <- topTable(fit, number = nrow(df), adjust.method = "fdr",
                            sort.by = "p")

        sign.table <- top.rank[which(top.rank$adj.P.Val <= 0.01 &
                                        abs(top.rank$logFC) > 2), ]

        if (nrow(sign.table) == 0) {

            cat(c("No statistical significant events for:",
                    names.meta[z],
                    "\n"))

            sign_table[z,1] <- analysis
            sign_table[z,2] <- "0"

            full.tads = list()
        }

        else {

            # annotate sign.table

            sign.table$ID <- row.names(sign.table)

            sign.table <- merge(sign.table,
                                data.all,
                                by.x = "ID",
                                by.y = "ID")

            sign.tad.info <- sign.table %>%
            dplyr::group_by(tad_name) %>%
            dplyr::summarise(count = n())



            # annotate tad.info table

            tad.info <- data.all %>%
            dplyr::group_by(tad_name) %>%
            dplyr::summarize(count = n())

            sign.tad.info <- merge(sign.tad.info,
                                    tad.info,
                                    by = "tad_name")

            sign.tad.info$freq <- sign.tad.info$count.x /
                sign.tad.info$count.y * 100

            sign.tad.info$pvalue <- 1

            for (i in 1:nrow(sign.tad.info)) {

                sign.tad.info$pvalue[i] <- 1 - phyper(sign.tad.info$count.x[i],
                                                    nrow(sign.table),
                                                    nrow(df) - nrow(sign.table),
                                                    sign.tad.info$count.y[i])

            }


            if(expr_data != FALSE){

                # get expression data

                expr <- data.all[which(data.all$parent == expr_data), ]

                df <- expr[,..sample.list]

                df <- as.data.frame(df)
                row.names(df) <- expr$ID

                # build model

                fit <- lmFit(object = df, design = phenoMat)

                gc()
                set.seed(6)
                fit <- eBayes(fit)

                # get top rank events

                gc()
                top.rank <- topTable(fit, number = nrow(df),
                                    adjust.method = "fdr")

                # any filtering ?

                # annotate sign.table (expression)

                sign.table.expr <- as.data.frame(top.rank)

                sign.table.expr$ID <- row.names(sign.table.expr)

                sign.table.expr <- merge(sign.table.expr,
                                        expr,
                                        by = "ID")

                sign.tad.expr <- sign.table.expr %>%
                dplyr::group_by(tad_name) %>%
                dplyr::summarise(mean_logFC = mean(abs(logFC)))



                # merge information

                tad.all.info <- merge(sign.tad.info,
                                        sign.tad.expr,
                                        by = "tad_name" )

                tad.all.info.f <- tad.all.info %>%
                    dplyr::filter(count.x > 4) %>%
                    dplyr::filter(pvalue < 0.01) %>%
                    dplyr::filter(freq > 15) %>%
                    dplyr::filter(mean_logFC > 2)

                sign_table[z,1] <- analysis
                sign_table[z,2] <- as.character(nrow(tad.all.info.f))



                # create output tables

                full.tads <- merge(tad.all.info.f,
                                    data.all,
                                    by = "tad_name")


                ################### CHECK IT   ##################################
                # if(nrow(full.tads) > 0) {
                Diff_list[[analysis]] = full.tads
                    # write.table(full.tads,
                    #             file = paste(output_folder, "/", analysis,
                    #                             "_TADiff.txt", sep = ""),
                    #             col.names = TRUE,
                    #             row.names = FALSE,
                    #             quote = FALSE,
                    #             sep = "\t")


                ##########################################################
            }

            # what's the value of `sign_table` if we don't have expression data
            # which are the outputs if we don't have expression data
        }
    }


    sign_table <- as.data.frame(sign_table)


    return(list(Diff_list,sign_table))
}


# TADiff_result <- TADiff(
#
# mapping_file = "/Users/aspaor/Downloads/bloodcancer/metaData_groups.csv",
#                        methylo_result = methylo_result,
#                        names.meta = c("IGHV","Gender"),
#                    expr_data = 2)
# TADiff_result


