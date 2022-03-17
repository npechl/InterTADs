#'
#' Data Integration part
#'
#' @param counts_folder 
#' A stringr-path to a repository containing files 
#' with NGS counts data (e.g. expression count table).
#' 
#' @param counts_fls 
#' A character vector of paths to files 
#' with NGS counts data (e.g. expression count table).
#' 
#' @param freq_folder 
#' A stringr-path to a repository containing files 
#' with NGS frequency data (e.g. DNA methylation table).
#' 
#' @param freq_fls 
#' A character vector of paths to files 
#' with NGS frequency data (e.g. DNA methylation table).
#' 
#' @param mapping_file 
#' A mapping file containing mapping betwen the columns of the input 
#' NGS datasets. The first column corresponds to the sample ID that 
#' will be used in the output file. The following columns correspond 
#' to the column names of the input datasets. 
#' The file also contains the related sample meta-data
#' 
#' @param tad_file 
#' A BED file containing information 
#' about the desired TADs' chromosomal location.
#' 
#' @param tech 
#' Human Genome Reference technology, either "hg19" or "hg38" can be used.
#'
#' @import data.table
#' @import systemPipeR
#' @import org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import annotables
#' @import IRanges 
#' @importFrom S4Vectors Rle queryHits subjectHits
#' @importFrom  GenomicRanges GRanges
#' @importFrom stringr str_detect str_replace str_to_lower str_remove
#'
#' @description
#'
#' @return 
#' A list of integratedTADtable and summary. 
#' IntegratedTADtable contains all data info about TADs and 
#' summary explains types of different events (methylation, expression etc.).
#'
#' @export
#'
#' @examples
# 
# # example 1
# result<- data_integration (
# 
#     counts_folder = system.file(
#         "extdata", "Datasets", "counts", package = "InterTADs"
#     ),
# 
#     freq_folder = system.file(
#         "extdata", "Datasets", "freq", package = "InterTADs"
#     ),
# 
#     mapping_file = system.file(
#         "extdata", "Datasets", "meta-data.csv", package = "InterTADs"
#     ),
# 
#     tad_file =system.file(
#         "extdata", "Datasets",
#         "hglft_genome_2dab_ec1330.bed", package = "InterTADs"
#     ),
# 
#     tech = "hg19"
# )
# 
# 
# # example 2
# result<- data_integration (
# 
#     counts_fls = list.files(
#         system.file(
#             "extdata", "Datasets", "counts", package = "InterTADs"
#         ),
# 
#         full.names = TRUE
#     ),
# 
#     freq_fls =  list.files(
#         system.file(
#             "extdata", "Datasets", "freq", package = "InterTADs"
#         ),
# 
#         full.names = TRUE
#     ),
# 
#     mapping_file = system.file(
#         "extdata", "Datasets", "meta-data.csv", package = "InterTADs"
#     ),
# 
#     tad_file =system.file(
#         "extdata", "Datasets",
#         "hglft_genome_2dab_ec1330.bed", package = "InterTADs"
#     ),
# 
#     tech = "hg19"
# )



data_integration <- function(
    counts_folder = NULL,
    counts_fls = NULL,

    freq_folder = NULL,
    freq_fls = NULL,

    mapping_file = mapping_file,
    tad_file = tad_file,
    tech = "hg38"
) {


    # Reading sample_metadata data file ---------------------
    sample_metadata <- fread(
        mapping_file, fill = TRUE
    )




    # Get input files ---------------------------

    if(is.null(freq_folder) & is.null(freq_fls)) {


    } else if(is.null(freq_fls)) {

        freq_fls = list.files(
            freq_folder,
            full.names = TRUE
        )

    } else {


    }



    if(is.null(counts_folder) & is.null(counts_fls)) {


    } else if(is.null(counts_fls)) {
        counts_fls = list.files(
            counts_folder,
            full.names = TRUE

        )

    } else {



    }


    biodata <- list()
    parent = 0

    # Reading input frequency tables ----------------------
    
    if(length(freq_fls) > 0) {

        for(i in 1:length(freq_fls)) {

            new <- fread(freq_fls[i], fill = TRUE)

            keep0 <- apply(
                sample_metadata[,2:ncol(sample_metadata)], 2, function(x, y) {

                    return(length(which(y %in% x)))

                }, colnames(new)
            )

            keep   <- names(which(keep0 == max(keep0))[1])
            keep   <- which(colnames(sample_metadata) == keep)
            # parent <- keep
            parent = parent + 1

            file_mapping <- sample_metadata[, c(1, keep), with = FALSE]

            file_mapping <- file_mapping[which(
                !(is.na( file_mapping[[2]] ))
            ), ]

            new <- cbind(new[,1:4], new[, file_mapping[[2]], with = FALSE])

            colnames(new) <- c(
                "ID",
                "chromosome_name",
                "start_position", "end_position",
                file_mapping[[1]]
            )

            new$chromosome_name <- str_to_lower(new$chromosome_name)
            new$chromosome_name <- str_remove(new$chromosome_name, "chr")

            if( max(new[,5:ncol(new)]) <= 1 ){

              new[, 5:ncol(new)] <- 100 * new[, 5:ncol(new)]

            }

            new$parent <- parent

            biodata[[i]] <- new

        }

    }

    # Reading input counts tables -----------------------
    if(length(counts_fls) > 0) {

        for(i in 1:length(counts_fls)){

            new <- fread(counts_fls[i], fill = TRUE)

            keep0 <- apply(
                sample_metadata[,2:ncol(sample_metadata)], 2, function(x, y) {

                    return(length(which(y %in% x)))

                }, colnames(new)
            )


            keep   <- names(which(keep0 == max(keep0))[1])
            keep   <- which(colnames(sample_metadata) == keep)
            # parent <- keep
            parent = parent + 1

            file_mapping <- sample_metadata[, c(1, keep), with = FALSE]

            file_mapping <- file_mapping[which(
                !(is.na( file_mapping[[2]] ))
            ), ]

            new <- cbind(new[,1:4], new[, file_mapping[[2]], with = FALSE])

            colnames(new) <- c(
                "ID",
                "chromosome_name",
                "start_position", "end_position",
                file_mapping[[1]]
            )

            new$chromosome_name <- str_to_lower(new$chromosome_name)
            new$chromosome_name <- str_remove(new$chromosome_name, "chr")

            temp    <- new[,5:ncol(new)]
            temp    <- log(temp + 1)
            col.max <- apply(temp, 2, max)
            col.max <- as.numeric(col.max)

            temp <- t(temp) * (100 / col.max)
            temp <- as.data.table(t(temp))

            new <- cbind(new[, 1:4], temp)

            new$parent <- parent

            biodata[[i + length(freq_fls)]] <- new
        }

    }

    biodata <- rbindlist(biodata, use.names = TRUE, fill = TRUE)


    # Filtering -----------------

    #  Keep chromosomes 1 - 22
    biodata <- biodata[which(
        biodata$chromosome_name %in% as.character(1:22)
    ), ]

    # Remove all-zeros records
    zeros   <- rowSums(biodata[, sample_metadata[[1]], with = FALSE])
    biodata <- biodata[which(zeros != 0), ]

    # Getting genomic features ---------------------------------


    # Get gene names (expressed as entrez ids) and locus
    # (e.g. exon, intron, cds etc.)
    # based on chromosomal location

    if(tech == "hg19"){

        feat <- genFeatures(
            TxDb.Hsapiens.UCSC.hg19.knownGene,
            reduce_ranges = FALSE,
            verbose = FALSE
        )

    } else {

        feat <- genFeatures(
            TxDb.Hsapiens.UCSC.hg38.knownGene,
            reduce_ranges = FALSE,
            verbose = FALSE
        )

    }


    feat <- unlist(feat)

    biodata_granges <- GRanges(
        seqnames = Rle(
            paste("chr", biodata$chromosome_name, sep = "")
        ),

        ranges = IRanges(
            start = as.numeric(biodata$start_position),
            end = as.numeric(biodata$end_position)
        )
    )

    overlaps <- findOverlaps(biodata_granges, feat)

    feat <- data.table(
        feature_by = as.character(feat$feature_by),
        featuretype = as.character(feat$featuretype)
    )

    feat <- cbind(
        biodata[from(overlaps), 1:4], feat[to(overlaps), ]
    )

    biodata <- cbind(
        feat, biodata[from(overlaps), 5:ncol(biodata)]
    )

    biodata <- unique(biodata)

    biodata[which(biodata$feature_by == "character(0)"), ]$feature_by <- NA
    biodata[which(biodata$feature_by == ""), ]$feature_by             <- NA

    # Converting entrez ids to hgnc symbols
    genes <- unique(biodata$feature_by)
    genes <- genes[which(genes != "NA")]
    genes <- genes[!str_detect(genes, "INTER")]

    mapping <- map_entrez_ids(entrez.ids = genes, tech = tech)

    map <- match(biodata$feature_by, mapping$entrez.id)
    biodata$feature_by <- mapping[map, ]$hgnc.symbol


    #
    # Collapsing genes and feature type
    # of same observation into same row
    #
    # e.g.
    #
    # index1 --- MYC, exon
    # index1 --- RUNX1T1, intron
    # ...........................
    # index1 --- MYC|RUNX1T1, exon|intron

    features <- biodata[, .(
        Gene_id = paste(feature_by, collapse = ","),
        Gene_feature = paste(featuretype, collapse = ",")
    ), by = ID]


    biodata <- biodata[which(!duplicated(biodata$ID)), ]

    who <- match(biodata$ID, features$ID)

    biodata$feature_by  <- features[who, ]$Gene_id
    biodata$featuretype <- features[who, ]$Gene_feature

    names <- colnames(biodata)
    names <- str_replace(names, "feature_by", "Gene_id")
    names <- str_replace(names, "featuretype", "Gene_locus")

    colnames(biodata) <- names

    # Annotate with TADs --------------------------------------

    TAD <- fread(
        tad_file,
        header = FALSE, sep = "\t"
    )

    # Reordering
    data.over <- biodata[, c(
        "ID", "chromosome_name", "start_position", "end_position"
    ), with = FALSE]

    data.over$chromosome_name <- paste0("chr", data.over$chromosome_name)
    data.over                 <- data.over[,c(2,3,4,1)]

    # Overlap of the TAD with events
    colnames(TAD) <- paste(c("chr", "start", "end", "name"))
    colnames(data.over) <- paste(c("chr", "start", "end", "name"))

    # Make GRanges object
    TAD_gr <- with(
        TAD, GRanges(
            chr, IRanges(start = start, end = end, names = name)
        )
    )

    biodata_gr <- with(
        data.over, GRanges(
            chr, IRanges(start = start, end = end, names = name)
        )
    )


    # Completely overlapping ----------------------------
    df0 <- findOverlaps(query = TAD_gr, subject = biodata_gr)
    df1 <- cbind(TAD[queryHits(df0),], data.over[subjectHits(df0),])
    df1 <- df1[,c(2,3,4,8)]

    colnames(df1) <- c("tad_start", "tad_end", "tad_name", "name.1")

    #' Overlapping with TADs' table
    #'
    full <- merge(df1, biodata, by.x = "name.1", by.y = "ID")

    colnames(full)[1] <- "ID"

    full <- full[, c(
        "chromosome_name",
        "tad_name",
        "tad_start",
        "tad_end",
        "ID",
        "start_position",
        "end_position",
        "Gene_id",
        "Gene_locus",
        "parent",
        sample_metadata[[1]]
    ), with = FALSE]

    #end_time <- Sys.time()

    # Generating outputs -----------------------------


    sum_table <- full[, by = parent, 
                      .(IDs = paste(sample(ID, 3), collapse = ", "))]


    return(
        list(
            "IntegratedTADtable" = full,
            "summary" = sum_table
        )
    )

}








