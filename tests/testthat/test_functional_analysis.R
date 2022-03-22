
test_that('Test functional analysis runs properly',{

    result<- data_integration (
        counts_folder = system.file("extdata", "Datasets",
                                "counts", package="InterTADs"),
        freq_folder = system.file("extdata", "Datasets",
                                "freq", package="InterTADs"),
        mapping_file = system.file("extdata", "Datasets",
                                "meta-data.csv", package="InterTADs"),

        tad_file =system.file("extdata", "Datasets",
                                "hglft_genome_2dab_ec1330.bed",
                                package="InterTADs"),
        tech = "hg19"
    )

    methylo_result <- prepare_methylation_values(
        integratedTADtable = result[[1]],
        mapping_file = system.file("extdata", "Datasets",
                                    "meta-data.csv", package="InterTADs"),
        meth_data = 2
    )
    TADiff_result <- TADiff(
        integratedTADtable = methylo_result,
        mapping_file = system.file("extdata", "Datasets",
                                    "meta-data.csv", package="InterTADs"),
        names.meta = c("group"),
        expr_data = 3,
        adj.PVal = 0.01,
        log_thr = 2,
        tad_event = 4,
        pval_thr = 0.01,
        freq_thr = 15,
        mean_logFC_thr = 2)

    result_evenDiff <- evenDiff(
        integratedTADtable = methylo_result,
        mapping_file = system.file("extdata", "Datasets",
                                    "meta-data.csv", package="InterTADs"),
        names.meta = c('group'),
        adj.PVal = 0.01,
        log_thr = 4)

        expect_true(
            functional_analysis(
                annotation_file = "/gencode.v19.annotation.gff3.gz",
                tech = "hg19",
                files_evenDiff = result_evenDiff[[1]],
                files_TADiff = TADiff_result[[1]],
                names.meta = c('group'),
                exp_parent = 3,
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
                dir_name =  system.file("extdata", "Datasets",
                                        package="InterTADs"),
                output_folder = paste0("Outputs")))



})
