test_that('Test prepare methylation values runs properly',{
    result <-data_integration (
    
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
  
    expect_vector(prepare_methylation_values(
            integratedTADtable = result[[1]],
            mapping_file = system.file(
                "extdata", "Datasets", "meta-data.csv", package = "InterTADs"),
            meth_data = 2
        )
    )
    
    expect_error(prepare_methylation_values(
            integratedTADtable = NULL,
            mapping_file = NULL,
            meth_data = NULL
      )
    )
  
  
})
