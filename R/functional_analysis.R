# Load libraries and functions -------------------

#rm(list = ls())

# source("R/libraries_enrich.R")
# source("R/motif_enrich.R")
# source("R/go_pathway_enrich.R")
# source("R/visualization_enrich.R")

# Inputs -----------------


#' Functional analysis.
#'
#' @description
#' Enrichment analysis with GO terms, KEGG pathways and motif enrichment
#' with TFs.
#' Optional step to be performed on the output files
#' of the `02c_evenDiff.R` and `02d_TADiff.R` scripts.
#'
#'
#' @param tech Human Genome Reference used
#'
#' @param exp_parent number of the parent file of the expression data
#'
#' @param dbs Databases used from EnrichR
#'
#' @param type the prevously selected databases acronyms used for the names
#'             of the outputs files e.g. GO_MF for GO_Molecular_Function_2018
#'
#' @param genes_cover gene coverage of the databases
#'
#' @param p_adjust_method p adjustment method, the methods supported are:
#'                        `c("holm", "hochberg", "hommel", "bonferroni",`
#'                        `"BH", "BY", "fdr", "none")`
#'
#' @param criterio Enrichr result column selected as criterio:
#'                `"p-value"` or `"adjusted p-value"`
#'
#' @param cut_off cut-off Enrichr enrichment (adjusted) p-value
#'
#' @param min_genes min number of genes in over-represented terms
#'
#' @param system RStudio supports different fonts for different
#' operating systems
#'
#' @param dir_name name or filepath of the input folder
#'
#' @param output_folder name or filepath of the output folder
#'
#' @description
#'
#' @return
#'
#' @export
#'
#' @examples

functional_analysis<- function(tech = "hg19",
                                exp_parent = 2,
                                dbs = c("GO_Molecular_Function_2018",
                                        "GO_Biological_Process_2018",
                                        "KEGG_2019_Human"),
                                type = c("GO_MF", "GO_BP", "KEGG"),
                                genes_cover = c(11459, 14433, 7802),
                                p_adjust_method = "fdr",
                                cut_off = 0.05,
                                criterio = "adjusted p-value",
                                min_genes = 3,
                                system = "win",
                                dir_name = "Datasets"){
                                #output_folder = paste0("Outputs_whatever")){

    start_time <- proc.time()
    # Set graph fonts --------------------
    # set_graph_fonts(system)

    # Download hg gff file ----------------
    if (tech == "hg19") {

        download.file(url="ftp://ftp.ebi.ac.uk/pub/databases/gencode
                    /Gencode_human/release_19/gencode.v19.annotation.gff3.gz",
                    destfile=paste0(dir_name,'/gencode.v19.annotation.gff3.gz'),
                    method='auto')

    } else if (tech == "hg38") {

        download.file(url="ftp://ftp.ebi.ac.uk/pub/databases/gencode
                    /Gencode_human/release_36/gencode.v36.annotation.gff3.gz",
                    destfile=paste0(dir_name,'/gencode.v36.annotation.gff3.gz'),
                    method='auto')

    }

    # Enrichment + Data Analysis -----------------

    dir.create(output_folder, showWarnings = FALSE)

    file.create(paste0(output_folder, "/time_info.txt"))

    data_type <- data.table(type = type,
                            cover = as.numeric(genes_cover))

    files_evenDiff <- list.files(dir_name, pattern = c("evenDiff"))
    files_TADiff   <- list.files(dir_name, pattern = c("TADiff"))

    # evenDiff files --------------

    if (!is_empty(files_evenDiff)) {

        for (one_file in files_evenDiff) {

            file_name <- str_remove(one_file, ".txt")
            folder <- create_folders(paste0(output_folder, "/", file_name))
            biodata <- fread(paste0(dir_name, "/", one_file))

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


            # Motif enrichment
            report_list <- motifs_enrich(biodata,
                                        folder$motif_outputs,
                                        dir_name,
                                        p_adjust_method,
                                        cut_off,
                                        tech,
                                        exp_parent)

            if (length(report_list)> 0) list_motif <- motif_outputs(report_list)

            # Output Files -------------

            # Enrichment all
            if (exists("data_all")) {
                fwrite(data_all,
                        paste0(folder$go_all_outputs,
                            "/over-represented GO terms-enrichment all.csv"),
                            row.names = FALSE, sep = "\t", quote = FALSE)
            }

            if (exists("list_GO_MF_all")) {

                fwrite(list_GO_MF_all$data_per_term,
                        paste0(folder$go_all_outputs,
                                "/GO MF Terms in different TADs.csv"),
                                row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_GO_MF_all$data_per_tad,
                        paste0(folder$go_all_outputs,
                            "/over-represented GO MF terms-enrichment all.csv"),
                            row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_GO_BP_all")) {

                fwrite(list_GO_BP_all$data_per_term,
                        paste0(folder$go_all_outputs,
                                "/GO BP Terms in different TADs.csv"),
                                row.names = FALSE, sep = "\t", quote = FALSE)


                fwrite(list_GO_BP_all$data_per_tad,
                        paste0(folder$go_all_outputs,
                            "/over-represented GO BP terms-enrichment all.csv"),
                            row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_KEGG_all")) {

                fwrite(list_KEGG_all$data_per_tad,
                    paste0(folder$kegg_all_outputs,
                        "/over-represented KEGG Pathways-enrichment all.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)


                fwrite(list_KEGG_all$data_per_term,
                        paste0(folder$kegg_all_outputs,
                                "/KEGG Pathways in different TADs.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }


            # Motif enrichment
            if (length(report_list) > 0) {

                fwrite(list_motif$table_per_tad,
                        paste0(folder$motif_outputs,
                                "/over-represented TFs in each tad.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_motif$table_per_tfs,
                    paste0(folder$motif_outputs, "/TFs in different TADs.csv"),
                    row.names = FALSE, sep = "\t", quote = FALSE)

                file.create(paste0(folder$motif_outputs,"/report MotifEA.txt"),
                            showWarnings = FALSE)

                dput(report_list, file = paste0(folder$motif_outputs,
                                                "/report MotifEA.txt"))
            }

            # Visualization -------------

            # Enrich all visualization
            if (exists("list_GO_MF_all")) {

                enrichr_visual(folder$go_all_images,
                                "GO MF Terms",
                                list_GO_MF_all$data_visual,
                                criterio)
            }

            if (exists("list_GO_BP_all")) {

                enrichr_visual(folder$go_all_images,
                                "GO BP Terms",
                                list_GO_BP_all$data_visual,
                                criterio)
            }

            if (exists("list_KEGG_all")) {

                enrichr_visual(folder$kegg_all_images,
                                "KEGG Pathways",
                                list_KEGG_all$data_visual,
                                criterio)
            }

            # Motif enrichment analysis visualization
            if (length(report_list) > 0) {

                motif_visual(folder$motif_images,
                            folder$motif_outputs,
                            list_motif$data_visual,
                            report_list,
                            criterio)
            }

            rm(list_motif, list_all, list_GO_MF_all, list_GO_BP_all,
                list_KEGG_all,data_all, report_list)
        }




    # TADiff files ----------------

    if (!is_empty(files_TADiff)) {

        for (one_file in files_TADiff) {

            file_name <- str_remove(one_file, ".txt")
            folder <- create_folders(paste0(output_folder, "/", file_name))
            biodata <- fread(paste0(dir_name, "/", one_file))
            # biodata <- biodata[which((biodata$significant == TRUE) &
            #(tad_name == "TAD2886")), ]

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



            # Motif enrichment
            report_list <- motifs_enrich(biodata,
                                            folder$motif_outputs,
                                            dir_name,
                                            p_adjust_method,
                                            cut_off,
                                            tech,
                                            exp_parent)

            if (length(report_list)> 0) list_motif <- motif_outputs(report_list)


            # Output Files --------------

            # Enrichment all
            if (exists("data_all")) {
                fwrite(data_all,
                        paste0(folder$go_all_outputs,
                            "/over-represented GO terms-enrichment all.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)
            }

            if (exists("list_GO_MF_all")) {
                fwrite(list_GO_MF_all$data_per_term,
                        paste0(folder$go_all_outputs,
                                "/GO MF Terms in different TADs.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_GO_MF_all$data_per_tad,
                        paste0(folder$go_all_outputs,
                            "/over-represented GO MF terms-enrichment all.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_GO_BP_all")) {

                fwrite(list_GO_BP_all$data_per_term,
                        paste0(folder$go_all_outputs,
                                "/GO BP Terms in different TADs.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)


                fwrite(list_GO_BP_all$data_per_tad,
                        paste0(folder$go_all_outputs,
                            "/over-represented GO BP terms-enrichment all.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_KEGG_all")) {

                fwrite(list_KEGG_all$data_per_tad,
                        paste0(folder$kegg_all_outputs,
                        "/over-represented KEGG Pathways-enrichment all.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)


                fwrite(list_KEGG_all$data_per_term,
                        paste0(folder$kegg_all_outputs,
                                "/KEGG Pathways in different TADs.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }

            # Enrichment per TAD
            if (exists("data_per_tad")) {
                fwrite(data_per_tad,
                        paste0(folder$go_per_outputs,
                        "/over-represented GO terms-enrichment per tad.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)
            }

            if (exists("list_GO_MF_per_tad")) {
                fwrite(list_GO_MF_per_tad$data_per_term,
                        paste0(folder$go_per_outputs,
                                "/GO MF Terms in different TADs.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_GO_MF_per_tad$data_per_tad,
                        paste0(folder$go_per_outputs,
                        "/over-represented GO MF terms-enrichment per TAD.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_GO_BP_per_tad")) {

                fwrite(list_GO_BP_per_tad$data_per_term,
                        paste0(folder$go_per_outputs,
                                "/GO BP Terms in different TADs.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_GO_BP_per_tad$data_per_tad,
                        paste0(folder$go_per_outputs,
                        "/over-represented GO BP terms-enrichment per TAD.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

            }

            if (exists("list_KEGG_per_tad")) {

                fwrite(list_KEGG_per_tad$data_per_tad,
                    paste0(folder$kegg_per_outputs,
                    "/over-represented KEGG Pathways-enrichment per TAD.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_KEGG_per_tad$data_per_term,
                        paste0(folder$kegg_per_outputs,
                                "/KEGG Pathways in different TADs.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)
            }


            # Motif enrichment
            if (length(report_list) > 0) {

                fwrite(list_motif$table_per_tad,
                        paste0(folder$motif_outputs,
                                "/over-represented TFs in each tad.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                fwrite(list_motif$table_per_tfs,
                    paste0(folder$motif_outputs, "/TFs in different TADs.csv"),
                        row.names = FALSE, sep = "\t", quote = FALSE)

                file.create(paste0(folder$motif_outputs,"/report MotifEA.txt"),
                            showWarnings = FALSE)

                dput(report_list, file = paste0(folder$motif_outputs,
                                                "/report MotifEA.txt"))
            }

            # Visualization --------------

            # Enrich all visualization
            if (exists("list_GO_MF_all")) {
                enrichr_visual(folder$go_all_images,
                                "GO MF Terms",
                                list_GO_MF_all$data_visual,
                                criterio)
            }

            if (exists("list_GO_BP_all")) {
                enrichr_visual(folder$go_all_images,
                                "GO BP Terms",
                                list_GO_BP_all$data_visual,
                                criterio)
            }

            if (exists("list_KEGG_all")) {
                enrichr_visual(folder$kegg_all_images,
                                "KEGG Pathways",
                                list_KEGG_all$data_visual,
                                criterio)
            }


            # Enrich per TAD visualization
            if (exists("list_GO_MF_per_tad")) {
                enrichr_visual(folder$go_per_images,
                                "GO MF Terms",
                                list_GO_MF_per_tad$data_visual,
                                criterio)
            }

            if (exists("list_GO_BP_per_tad")) {
                enrichr_visual(folder$go_per_images,
                                "GO BP Terms",
                                list_GO_BP_per_tad$data_visual,
                                criterio)
            }

            if (exists("list_KEGG_per_tad")) {
                enrichr_visual(folder$kegg_per_images,
                                "KEGG Pathways",
                                list_KEGG_per_tad$data_visual,
                                criterio)
            }


            # Motif enrichment analysis visualization
            if (length(report_list) > 0) {
                motif_visual(folder$motif_images,
                            folder$motif_outputs,
                            list_motif$data_visual,
                            report_list,
                            criterio)
                }

            rm(list_all, list_GO_MF_all, list_GO_BP_all, list_KEGG_all,
                data_all,list_per_tad, list_GO_MF_per_tad, list_GO_BP_per_tad,
                list_KEGG_per_tad, data_per_tad, list_motif, report_list)
            }

        }

    }
    return(NULL)
}
# # Time measurements ----------
#
# end_time <- proc.time()
#
# time_diff <- end_time - start_time
#
# line <- paste0("Total time in secs")
# write(line, paste0(output_folder, "/time_info.txt"), append = T)
#
# line <- paste0("user: ", time_diff[1])
# write(line, paste0(output_folder,"/time_info.txt"), append = T)
#
# line <- paste0("system: ", time_diff[2])
# write(line, paste0(output_folder,"/time_info.txt"), append = T)
#
# line <- paste0("elapsed: ", time_diff[3])
# write(line, paste0(output_folder,"/time_info.txt"), append = T)
