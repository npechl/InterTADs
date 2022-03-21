
#' evenDiff
#'
#' @param names.meta
#' meta data columns to process (names or indexes).
#' @param integratedTADtable
#' IntegratedTADtable contains all data info about TADs.
#' This is the output integrated table
#' from data_integration, prepare_methylation or
#' exnsembl_ids functions.
#' @param adj.PVal
#' Significant events adjusted p value theshold. Defaults to 0.05.
#' @param log_thr
#' LogFC theshold. Defaults to 2.
#' @param seed
#' A number used to initialize a pseudorandom number generator.
#' @param mapping_file
#' A mapping file containing mapping betwen the columns of the input
#' NGS datasets. The first column corresponds to the sample ID that
#' will be used in the output file. The following columns correspond
#' to the column names of the input datasets.
#' The file also contains the related sample meta-data
#'
#' @import data.table
#' @import limma
#'
#' @description
#' evenDiff performs differential analysis of multi-omics
#' events (expression, methylation etc.).
#'
#' @return
#' A list of data tables containing significant events along with
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
#' result_evenDiff <- evenDiff(
#'     integratedTADtable = result_ensmbl,
#'     mapping_file = system.file(
#'         "extdata", "Datasets", "meta-data.csv", package = "InterTADs"
#'     ),
#'     names.meta = c('group'),
#'     adj.PVal = 0.01,
#'     log_thr = 4
#' )


evenDiff <- function(
    integratedTADtable,
    mapping_file,
    names.meta,
    adj.PVal = 0.05,
    log_thr = 2,
    seed = 6
) {

    integratedTADtable$ID <- paste(
        integratedTADtable$tad_name, integratedTADtable$ID, sep = ";"
    )

    meta <- fread(mapping_file)

    sample.list <- meta[[1]]

    gene.parents <- integratedTADtable[, c(
        "ID", "Gene_id", "parent"
    ), with = FALSE]

    colnames(gene.parents) <- c("ID", "gene_ID", "parents")


    diffevent <- matrix(
        data = "0",
        nrow = length(names.meta),
        ncol = (length(unique(integratedTADtable$parent)) + 1)
    )

    colnames(diffevent) <- c(
        "Analysis", paste0("p_", unique(integratedTADtable$parent))
    )

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
        colnames(phenoMat) <- sub("^pheno", "", colnames(phenoMat))

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
            top.rank$adj.P.Val <= adj.PVal &
                abs(top.rank$logFC) > log_thr
        ), ]

        if (nrow(sign.table) == 0) {

            message(c(
                "No statistical significant events for:",
                names.meta[z]
            ))

            diffevent[z,1] <- analysis

        } else {

            # annotate sign.table

            sign.table$ID <- row.names(sign.table)

            sign.genes <- merge(
                gene.parents,
                sign.table,
                by.x = "ID",
                by.y = "ID"
            )


            sign.table <- merge(
                sign.table,
                integratedTADtable,
                by.x = "ID",
                by.y = "ID"
            )

            #sign <- sign.genes %>% dplyr::count(parents)
            sign <- sign.genes[,.N,by = parents]


            diffevent[z, 1] <- analysis
            diffevent[z, paste("p_", sign$parents, sep = "")] <- sign$n

            Diff_list[[analysis]] = sign.table

        }
    }

    diffevent <- as.data.frame(diffevent)

    for(i in 2:ncol(diffevent)){

        diffevent[,i] <- as.numeric(diffevent[,i])

    }


    diffevent$`total events` <- rowSums(diffevent[,2:ncol(diffevent)])



    return(
        list(
            Diff_list,
            diffevent
        )
    )

}


