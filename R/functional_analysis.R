
#' Functional analysis.
#'
#' @description
#' Enrichment analysis with GO terms, KEGG pathways and motif enrichment
#' with TFs.
#' Optional step to be performed on the output files
#' of the `02c_evenDiff.R` and `02d_TADiff.R` scripts.
#'
#'
#' @param annotation_file
#' @param tech Human Genome Reference used
#' @param files_evenDiff
#' @param files_TADiff
#' @param names.meta meta data columns to process (names or indexes)
#' @param exp_parent number of the parent file of the expression data
#' @param dbs Databases used from EnrichR
#' @param type the prevously selected databases acronyms used for the names
#'             of the outputs files e.g. GO_MF for GO_Molecular_Function_2018
#' @param genes_cover gene coverage of the databases
#' @param p_adjust_method p adjustment method, the methods supported are:
#'                        `c("holm", "hochberg", "hommel", "bonferroni",`
#'                        `"BH", "BY", "fdr", "none")`
#' @param cut_off cut-off Enrichr enrichment (adjusted) p-value
#' @param criterio Enrichr result column selected as criterio:
#'                `"p-value"` or `"adjusted p-value"`
#' @param min_genes  min number of genes in over-represented terms
#' @param system RStudio supports different fonts for different
#' operating systems
#' @param dir_name name or filepath of the input folder
#' @param output_folder name or filepath of the output folder
#'
#' @import data.table
#' @import stats
#' @import systemPipeR
#' @import GenomicRanges
#' @import ape
#' @import enrichR
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
#' TADiff_result <- TADiff(
#' integratedTADtable = methylo_result,
#' mapping_file = system.file("extdata", "Datasets",
#' "meta-data.csv", package="InterTADs"),
#' names.meta = c("group"),
#' expr_data = 3,
#' adj.PVal = 0.01,
#' log_thr = 2,
#' tad_event = 4,
#' pval_thr = 0.01,
#' freq_thr = 15,
#' mean_logFC_thr = 2)
#'
#' result_evenDiff <- evenDiff(
#' integratedTADtable = methylo_result,
#' mapping_file = system.file("extdata", "Datasets",
#' "meta-data.csv", package="InterTADs"),
#' names.meta = c('group'),
#' adj.PVal = 0.01,
#' log_thr = 4)
#'
#' func_anal_res <- functional_analysis(
#' annotation_file = "/gencode.v36.annotation.gff3.gz",
#' tech = "hg38",
#' files_evenDiff = result_evenDiff[[1]],
#' files_TADiff = TADiff_result[[1]],
#' names.meta =NULL,
#' exp_parent = 3,
#' dbs = c("GO_Molecular_Function_2018",
#' "GO_Biological_Process_2018",
#' "KEGG_2019_Human"),
#' type = c("GO_MF", "GO_BP", "KEGG"),
#' genes_cover = c(11459, 14433, 7802),
#' p_adjust_method = "fdr",
#' cut_off = 0.05,
#' criterio = "adjusted p-value",
#' min_genes = 3,
#' system = "win",
#' dir_name =  system.file("extdata", "Datasets", package="InterTADs"),
#' output_folder = paste0("Outputs"))




functional_analysis<- function(annotation_file = NULL,
                               tech = NULL,
                               files_evenDiff = NULL,
                               files_TADiff = NULL,
                               names.meta =NULL,
                                exp_parent = NULL,
                                dbs = NULL,
                                type = NULL,
                                genes_cover = NULL,
                                p_adjust_method = NULL,
                                cut_off = NULL,
                                criterio = NULL,
                                min_genes = NULL,
                                system = NULL,
                                dir_name = NULL,
                                output_folder = NULL){


    # Set graph fonts --------------------
    # set_graph_fonts(system)
    #options(timeout=1000)


    if (tech == "hg19") {

        annot_check<- str_detect('/gencode.v19.annotation.gff3.gz'
                                 ,annotation_file)
        if (!isTRUE(annot_check)){ stop('Please check tech and annotation file
                                        to be compatible')}

    } else if (tech == "hg38") {

        annot_check <- str_detect('/gencode.v36.annotation.gff3.gz',
                                  annotation_file)
        if (!isTRUE(annot_check)){ stop('Please check tech and annotation file
                                        to be compatible')}

    }

    # Enrichment + Data Analysis -----------------

    dir.create(output_folder, showWarnings = FALSE)

    #file.create(paste0(output_folder, "/time_info.txt"))

    data_type <- data.table(type = type,
                            cover = as.numeric(genes_cover))



    # evenDiff files --------------

    if (length(files_evenDiff)!=0) { ##!rlang::is_empty

        for (each_file in 1:length(files_evenDiff)) {

            evenDiff_file = names.meta[[each_file]]

            folder <- paste0(output_folder, "/evenDiff_",evenDiff_file)
            dir.create(folder, showWarnings = FALSE)
            #out <- paste0(output_folder, "/evenDiff_", evenDiff_file)
            #folder <- create_folders(out,result_evenDiff[[1]])
            biodata <- files_evenDiff[[each_file]] ##fread(paste0(dir_name, "/", one_file))

            if (length(biodata$ID) ==0 ){
                next
            }

            # Enrichment all
            list_all <- enrich_all(biodata, dbs, cut_off, criterio, type)

            # Data analysis
            for (l in c(1:length(dbs))) {

                if (nrow(list_all[[l]]) > 0) {

                    name <- names(list_all)[l]
                    genes_coverage <- data_type[type == name, cover]
                    result <- data_analysis(list_all[[l]],
                                            name,
                                            list_all$data_with_genes,
                                            genes_coverage,
                                            p_adjust_method,
                                            min_genes)

                    if (!is.null(result)) assign(paste0('list_', name, '_all'),
                                                result)
                }

            }

            if (exists("list_GO_MF_all") & exists("list_GO_BP_all")) {

                # Join GO Molecular Function and Biological Process outputs
                data_all <- full_join(list_GO_MF_all$data_per_tad,
                                        list_GO_BP_all$data_per_tad,
                                        by = "TAD")
            }


            if (length(biodata[,which(biodata$parent == exp_parent)]) == 0){
                next
            }

            # Motif enrichment
            report_list <- motifs_enrich(biodata,
                                        folder ,
                                        dir_name,
                                        p_adjust_method,
                                        cut_off,
                                        tech,
                                        exp_parent,
                                        annotation_file)

            if (length(report_list)> 0) list_motif <- motif_outputs(report_list)

            # Output Files -------------

            # Enrichment all
            if (exists("data_all")) {
                fwrite(data_all,
                        paste0(folder ,
                            "/over-represented GO terms-enrichment all.txt"),
                            row.names = FALSE, sep = "\t", quote = FALSE)
            }

            if (exists("list_GO_MF_all")) {

                fwrite(list_GO_MF_all$data_per_term,
                        paste0(folder ,
                                "/GO MF Terms in different TADs.txt"),
                                row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_GO_MF_all$data_per_tad,
                        paste0(folder ,
                            "/over-represented GO MF terms-enrichment all.txt"),
                            row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_GO_BP_all")) {

                fwrite(list_GO_BP_all$data_per_term,
                        paste0(folder ,
                                "/GO BP Terms in different TADs.txt"),
                                row.names = FALSE, sep = "\t", quote = FALSE)


                fwrite(list_GO_BP_all$data_per_tad,
                        paste0(folder ,
                            "/over-represented GO BP terms-enrichment all.txt"),
                            row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_KEGG_all")) {

                fwrite(list_KEGG_all$data_per_tad,
                    paste0(folder ,
                        "/over-represented KEGG Pathways-enrichment all.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)


                fwrite(list_KEGG_all$data_per_term,
                        paste0(folder ,
                                "/KEGG Pathways in different TADs.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }


            # Motif enrichment
            if (length(report_list) > 0) {

                fwrite(list_motif$table_per_tad,
                        paste0(folder ,
                                "/over-represented TFs in each tad.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_motif$table_per_tfs,
                    paste0(folder , "/TFs in different TADs.txt"),
                    row.names = FALSE, sep = "\t", quote = FALSE)

                file.create(paste0(folder ,"/report MotifEA.txt"),
                            showWarnings = FALSE)

                dput(report_list, file = paste0(folder ,
                                                "/report MotifEA.txt"))
            }



        }


    }

    # TADiff files ----------------

    if (length(files_TADiff)!=0) { ##!is_empty

        for (each_file in 1:length(files_TADiff)) {

            each_file <- 2
            TADiff_file = names.meta[[each_file]]

            folder <- paste0(output_folder, "/TAD_", TADiff_file)
            biodata <- files_evenDiff[[each_file]]
            dir.create(folder,showWarnings = FALSE)

            biodata <- files_TADiff[[each_file]]
            # biodata <- biodata[which((biodata$significant == TRUE) &
            #(tad_name == "TAD2886")), ]


            if (length(biodata$ID) ==0 ){
                next
            }
            # Enrichment all
            list_all <- enrich_all(biodata, dbs, cut_off, criterio, type)

            # Data analysis
            for (l in c(1:length(dbs))) {

                if (nrow(list_all[[l]]) > 0) {

                    name <- names(list_all)[l]
                    genes_coverage <- data_type[type == name, cover]
                    result <- data_analysis(list_all[[l]],
                                            name,
                                            list_all$data_with_genes,
                                            genes_coverage,
                                            p_adjust_method,
                                            min_genes)

                    if (!is.null(result)) assign(paste0('list_', name, '_all'),
                                                    result)
                }

            }

            if (exists("list_GO_MF_all") & exists("list_GO_BP_all")) {

                # join GO Molecular Function and Biological Process outputs
                data_all <- full_join(list_GO_MF_all$data_per_tad,
                                        list_GO_BP_all$data_per_tad,
                                        by = "TAD")
            }


            # Enrichment per TAD
            list_per_tad <- enrich_per_tad(biodata, dbs, cut_off, criterio,
                                            type)

            for (l in c(1:length(dbs))) {

                if (nrow(list_per_tad[[l]]) > 0) {

                    name <- names(list_per_tad)[l]
                    genes_coverage <- data_type[type == name, cover]
                    result <- data_analysis(list_per_tad[[l]],
                                            name,
                                            list_per_tad$data_with_genes,
                                            genes_coverage,
                                            p_adjust_method,
                                            min_genes)

                    if (!is.null(result)) assign(paste0('list_', name,
                                                        '_per_tad'),result)
                }
            }


            if (exists("list_GO_MF_per_tad") & exists("list_GO_BP_per_tad")){

                # join GO Molecular Function and Biological Process outputs
                data_per_tad <- full_join(list_GO_MF_per_tad$data_per_tad,
                                            list_GO_BP_per_tad$data_per_tad,
                                            by = "TAD")
            }


            if (length(biodata[,which(biodata$parent == exp_parent)]) == 0){
                next
            }


            # Motif enrichment
            report_list <- motifs_enrich(biodata,
                                            folder,
                                            dir_name,
                                            p_adjust_method,
                                            cut_off,
                                            tech,
                                            exp_parent,
                                         annotation_file)

            if (length(report_list)> 0) list_motif <- motif_outputs(report_list)


            # Output Files --------------

            # Enrichment all
            if (exists("data_all")) {
                fwrite(data_all,
                        paste0(folder,
                            "/over-represented GO terms-enrichment all.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)
            }

            if (exists("list_GO_MF_all")) {
                fwrite(list_GO_MF_all$data_per_term,
                        paste0(folder,
                                "/GO MF Terms in different TADs.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_GO_MF_all$data_per_tad,
                        paste0(folder,
                            "/over-represented GO MF terms-enrichment all.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_GO_BP_all")) {

                fwrite(list_GO_BP_all$data_per_term,
                        paste0(folder,
                                "/GO BP Terms in different TADs.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)


                fwrite(list_GO_BP_all$data_per_tad,
                        paste0(folder,
                            "/over-represented GO BP terms-enrichment all.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_KEGG_all")) {

                fwrite(list_KEGG_all$data_per_tad,
                        paste0(folder,
                        "/over-represented KEGG Pathways-enrichment all.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)


                fwrite(list_KEGG_all$data_per_term,
                        paste0(folder,
                                "/KEGG Pathways in different TADs.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }

            # Enrichment per TAD
            if (exists("data_per_tad")) {
                fwrite(data_per_tad,
                        paste0(folder,
                        "/over-represented GO terms-enrichment per tad.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)
            }

            if (exists("list_GO_MF_per_tad")) {
                fwrite(list_GO_MF_per_tad$data_per_term,
                        paste0(folder,
                                "/GO MF Terms in different TADs.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_GO_MF_per_tad$data_per_tad,
                        paste0(folder,
                        "/over-represented GO MF terms-enrichment per TAD.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_GO_BP_per_tad")) {

                fwrite(list_GO_BP_per_tad$data_per_term,
                        paste0(folder,
                                "/GO BP Terms in different TADs.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_GO_BP_per_tad$data_per_tad,
                        paste0(folder,
                        "/over-represented GO BP terms-enrichment per TAD.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_KEGG_per_tad")) {

                fwrite(list_KEGG_per_tad$data_per_tad,
                    paste0(folder,
                    "/over-represented KEGG Pathways-enrichment per TAD.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_KEGG_per_tad$data_per_term,
                        paste0(folder,
                                "/KEGG Pathways in different TADs.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)
            }


            # Motif enrichment
            if (length(report_list) > 0) {

                fwrite(list_motif$table_per_tad,
                        paste0(folder,
                                "/over-represented TFs in each tad.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_motif$table_per_tfs,
                    paste0(folder, "/TFs in different TADs.txt"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                file.create(paste0(folder,"/report MotifEA.txt"),
                            showWarnings = FALSE)

                dput(report_list, file = paste0(folder,
                                                "/report MotifEA.txt"))
            }


            }

        }


    return(TRUE)
}






