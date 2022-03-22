
#' TADiff
#'
#' @param names.meta
#' meta data columns to process (names or indexes)
#' @param expr_data
#' Parent index of expression data.
#' @param integratedTADtable
#' IntegratedTADtable contains all data info about TADs.
#' This is the output integrated table
#' from data_integration, prepare_methylation or
#' exnsembl_ids functions.
#' @param adj.PVal
#' Significant events adjusted p value theshold. Defaults to 0.05.
#' @param log_thr
#' LogFC theshold. Defaults to 2
#' @param tad_event
#' Number of events for each TAD. Defaults to 4.
#' @param pval_thr
#' Significant TADs p value threshold. Defaults to 0.05.
#' @param freq_thr
#' Percentage theshold of significant events to all events for each TAD, Defaults to 20 (20%).
#' @param mean_logFC_thr
#' Mean logFC threshold for each TAD.
#' @param seed
#' A number used to initialize a pseudorandom number generator.
#' @param mapping_file
#' A mapping file containing mapping betwen the columns of the input
#' NGS datasets. The first column corresponds to the sample ID that
#' will be used in the output file. The following columns correspond
#' to the column names of the input datasets.
#' The file also contains the related sample meta-data
#'
#'
#' @import data.table
#' @import limma
#'
#' @description
#' TADiff performs differential analysis of events (expression, methylation etc.)
#' by taking into account the chromatin configuration of the genome,
#' i.e. the topologically associating domains (TADs).
#'
#' @return
#' A list of data tables containing significant TADs along with
#' a summary file
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
#'
#' result_ensmbl <- ensembl_ids(
#'     integratedTADtable = methylo_result,
#'     expr_data = 3
#' )
#'
#' TADiff_result <- TADiff(
#'
#'     integratedTADtable = result_ensmbl,
#'     mapping_file = system.file(
#'         "extdata", "Datasets", "meta-data.csv", package = "InterTADs"
#'     ),
#'     names.meta = c("group"),
#'     expr_data = 3,
#'     adj.PVal = 0.01,
#'     log_thr = 2,
#'     tad_event = 4,
#'     pval_thr = 0.01,
#'     freq_thr = 20,
#'     mean_logFC_thr = 2
#' )


TADiff <- function(
    integratedTADtable,
    mapping_file,
    names.meta,
    expr_data = NULL,
    adj.PVal = 0.05,
    log_thr = 2,
    tad_event = 4,
    pval_thr = 0.05,
    freq_thr = 15,
    mean_logFC_thr = 2,
    seed = 6
) {

    integratedTADtable$ID <- paste(
        integratedTADtable$tad_name, integratedTADtable$ID, sep = ";"
    )


    meta <- fread(mapping_file)

    sample.list <- meta[[1]]
    integratedTADtable$Mean <- rowMeans(
        integratedTADtable[ , sample.list, with = FALSE]
    )

    integratedTADtable <- integratedTADtable[which(
        round(integratedTADtable$Mean) > 10
    ), ]

    sign_table <- matrix(0, nrow = length(names.meta), ncol = 2)

    Diff_list <- list()

    for (z in 1:length(names.meta)){

        message(names.meta[z])

        analysis <- names.meta[z]

        groups <- unique(as.character(meta[[analysis]]))
        groups <- groups[which(groups != "")]
        groups <- groups[!is.na(groups)]

        meta.keep <- meta[which(
            meta[[analysis]] == groups[1] | meta[[analysis]] == groups[2]
        ), ]

        sample.list <- meta.keep[[1]]

        df <- integratedTADtable[, sample.list, with = FALSE]

        df <- setDF(df, rownames = integratedTADtable$ID)

        pheno <- as.factor(meta.keep[[analysis]])


        phenoMat <- model.matrix(~ pheno)

        colnames(phenoMat) <- sub(
            "^pheno", "", colnames(phenoMat)
        )

        fit <- lmFit(object = df, design = phenoMat)

        set.seed(6)
        fit <- eBayes(fit)

        top.rank <- topTable(
            fit,
            number = nrow(df),
            adjust.method = "fdr",
            sort.by = "p"
        )

        sign.table <- top.rank[which(
            top.rank$adj.P.Val <= adj.PVal & abs(top.rank$logFC) > log_thr
        ), ]

        if (nrow(sign.table) == 0) {

            message(c(
                "No statistical significant events for:", names.meta[z]
            ))

            sign_table[z,1] <- analysis
            sign_table[z,2] <- "0"

            full.tads = list()

        } else {

            sign.table$ID <- row.names(sign.table)

            sign.table <- merge(
                sign.table,
                integratedTADtable,
                by.x = "ID",
                by.y = "ID"
            )

            sign.tad.info <- sing.table[, by = tad_name, .(count = .N)]

            # annotate tad.info table

            tad.info <- integratedTADtable[, by = tad_name, .(count = .N)]

            sign.tad.info <- merge(
                sign.tad.info,
                tad.info,
                by = "tad_name"
            )

            sign.tad.info$freq <- sign.tad.info$count.x /
                sign.tad.info$count.y * 100

            sign.tad.info$pvalue <- 1

            for (i in 1:nrow(sign.tad.info)) {

                sign.tad.info$pvalue[i] <- 1 - phyper(
                    sign.tad.info$count.x[i],
                    nrow(sign.table),
                    nrow(df) - nrow(sign.table),
                    sign.tad.info$count.y[i]
                )

            }


            if(!is.null(expr_data)){

                # get expression data

                expr <- integratedTADtable[which(
                    integratedTADtable$parent == expr_data
                ), ]

                df <- expr[, sample.list, with = FALSE]

                df <- setDF(df, rownames = expr$ID)

                # build model

                fit <- lmFit(object = df, design = phenoMat)

                set.seed(6)
                fit <- eBayes(fit)

                # get top rank events

                top.rank <- topTable(
                    fit,
                    number = nrow(df),
                    adjust.method = "fdr"
                )

                # any filtering ?

                # annotate sign.table (expression)

                sign.table.expr <- as.data.frame(top.rank)

                sign.table.expr$ID <- row.names(sign.table.expr)

                sign.table.expr <- merge(
                    sign.table.expr,
                    expr,
                    by = "ID"
                )

                sign.tad.expr <- sign.table.expr[, by = tad_name, .(
                    mean_logFC = mean(abs(logFC))
                )]

                # merge information

                tad.all.info <- merge(
                    sign.tad.info,
                    sign.tad.expr,
                    by = "tad_name"
                )

                tad.all.info.f <- tad.all.info[which(
                    count.x > tad_event &
                        pvalue < pval_thr &
                        freq > freq_thr &
                        mean_logFC > mean_logFC_thr
                ), ]

                sign_table[z,1] <- analysis
                sign_table[z,2] <- as.character(nrow(tad.all.info.f))



                # create output tables

                full.tads <- merge(
                    tad.all.info.f,
                    integratedTADtable,
                    by = "tad_name"
                )

                Diff_list[[analysis]] = full.tads

            } else {

                Diff_list[[analysis]] = sign.tad.info

            }

            # what's the value of `sign_table` if we don't have expression data
            # which are the outputs if we don't have expression data
        }
    }


    sign_table <- as.data.frame(sign_table)


    return(
        list(
            Diff_list,
            sign_table
        )
    )
}





