test_that('Test evendiff runs properly',{
    result <- data_integration (

        counts_folder = system.file(
            "extdata", "Datasets", "counts", package = "InterTADs"
        ),
    
        freq_folder = system.file(
            "extdata", "Datasets", "freq", package = "InterTADs"
        ),

        mapping_file = system.file(
            "extdata", "Datasets", "meta-data.csv", package = "InterTADs"
        ),

        tad_file =system.file(
            "extdata", "Datasets",
            "hglft_genome_2dab_ec1330.bed", package = "InterTADs"
        ),

        tech = "hg19")


    methylo_result <- prepare_methylation_values(
        integratedTADtable = result[[1]],
        mapping_file = system.file(
            "extdata", "Datasets", "meta-data.csv", package = "InterTADs"
        ),
        meth_data = 2
    )

    result_ensmbl <- ensembl_ids(
        integratedTADtable = methylo_result,
        expr_data = 3
    )

    expect_vector( evenDiff(
        integratedTADtable = result_ensmbl,
        mapping_file = system.file(
            "extdata", "Datasets", "meta-data.csv", package = "InterTADs"
        ),
        names.meta = c('group'),
        adj.PVal = 0.01,
        log_thr = 4,
        seed = 6
    ))
})
