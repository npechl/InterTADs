test_that('Test TADiff runs properly',{
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
        
        tech = "hg19"
    )


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

    expect_vector(TADiff(
  
        integratedTADtable = result_ensmbl,
        mapping_file = system.file(
            "extdata", "Datasets", "meta-data.csv", package = "InterTADs"),
        names.meta = c("group"),
        expr_data = 3,
        adj.PVal = 0.01,
        log_thr = 2,
        tad_event = 4,
        pval_thr = 0.01,
        freq_thr = 20,
        mean_logFC_thr = 2))
})
