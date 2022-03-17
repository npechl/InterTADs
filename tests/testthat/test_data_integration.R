test_that('Test data integration runs properly',{
    expect_vector( data_integration (

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
    ))
  
    expect_error( data_integration (
    
        counts_folder = NULL,
    
        freq_folder = NULL,
    
        mapping_file =NULL,
    
        tad_file =NULL,
    
        tech =NULL
    ))
  
  
})
