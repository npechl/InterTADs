#' This file contains functions used in the `03_enrichmentAnalysis.R`
#' script.



#' This function is called by the `functional_analysis.R` script.
#' motifs_enrich
#'
#' @param biodata
#' @param motif_output_folder
#' @param dir_name
#' @param p_adjust_method
#' @param cut_off
#' @param tech
#' @param exp_parent
#'
#' @description
#' It performs enrichment analysis using the PWMEnrich tool.
#' PWMEnrich input is the DNA sequences grouped per TAD.
#'
#'
#' @return
#'
#' @export
#'
#' @examples
motifs_enrich <- function(biodata,
                            motif_output_folder,
                            dir_name,
                            p_adjust_method,
                            cut_off,
                            tech,
                            exp_parent) {

    motif_data <- prepare_sequences(biodata,
                                    tech,
                                    motif_output_folder,
                                    exp_parent,
                                    dir_name)

    seq_tad_number <- get_dna_sequences(motif_data,
                                        motif_output_folder,
                                        tech)

    # Perform motif enrichment analysis using PWMEnrich.
    report_list <- list()

    # Load the pre-compiled lognormal background.
    data(PWMLogn.hg19.MotifDb.Hsap)

    k <- 1

    for (i in c(1:nrow(seq_tad_number))) {

        sequence <- readDNAStringSet(paste0(motif_output_folder,
                                            "/seq_perTADs.fasta"),
                                    format="fasta",
                                    skip = (seq_tad_number$start[i] - 1),
                                    nrec = (seq_tad_number$end[i] -
                                            seq_tad_number$start[i] + 1))

        res <- motifEnrichment(sequence, PWMLogn.hg19.MotifDb.Hsap)
        report <- groupReport(res, by.top.motifs = TRUE)

        report@d$adjusted.p.value <- p.adjust(report@d$p.value,
                                                method = p_adjust_method)

        report <- report[report$adjusted.p.value < cut_off]

        if (nrow(report@d) > 0) {

            report_list[[k]] <- report@d
            report_list[[k]]$tad <- seq_tad_number$tad[i]
            k <- k +1

        }
    }

    rm(motif_data, report, seq_tad_number, sequence, res)

    return(report_list)

}



#' This function is called by the `functional_analysis.R` script.
#' motif_outputs
#'
#' @param report_list
#'
#' @description
#' It manipulates the enriched data after the analysis and creates
#' three output data.tables to be used for the Output csv files and
#' the visualization.
#'
#' @return
#'
#' @export
#'
#' @examples
#'
#'

motif_outputs <- function(report_list) {

    csv_per_tad <- data.table(tad = character(),
                                tfs = character(),
                                motifs = character(),
                                adjusted_p_value = character(),
                                p_value = character())

    csv_per_tf <- data.table(tfs = character(),
                                tad = character(),
                                adjusted_p_value = character(),
                                motifs = character())

    top_motifs <- data.table(target = character(),
                            adjusted_p_value = numeric(),
                            id = character())

    for (i in c(1:length(report_list))) {

        temp <- report_list[[i]]

        # Exclude uncharacterized motifs
        temp <- temp[str_detect(temp$id, "UW.Motif.", negate = TRUE), ]

        per_tad <- temp %>%
            dplyr::select(tad, target, id, adjusted.p.value, p.value) %>%
            dplyr::summarise(tad = tad,
                            tfs = paste(target, collapse = "|"),
                            motifs = paste(id, collapse = "|"),
                            adjusted_p_value = paste(adjusted.p.value,
                                                    collapse = "|"),
                            p_value = paste(p.value, collapse = "|"), ) %>%
            as.data.table() %>%
            unique()

        csv_per_tad <- rbind(csv_per_tad, per_tad)

        per_tf_top_motif <- temp %>%
            dplyr::select(target, adjusted.p.value, id)

        colnames(per_tf_top_motif) <- c("target", "adjusted_p_value", "id")

        top_motifs <- rbind(top_motifs, per_tf_top_motif)

        per_tf <- temp %>%
            dplyr::select(target, tad, adjusted.p.value, id) %>%
            group_by(target,tad) %>%
            dplyr::summarise(target,
                            tad,
                            adjusted_p_value = paste(adjusted.p.value,
                                                    collapse = "|"),
                            motifs = paste(id, collapse = "|"), ) %>%
            as.data.table() %>%
            unique()

        colnames(per_tf) <- c("tfs", "tad", "adjusted_p_value", "motifs")
        csv_per_tf <- rbind(csv_per_tf, per_tf)
    }

    data_visual <- csv_per_tf %>%
        dplyr::select(tfs, tad, adjusted_p_value)

    csv_per_tf <- csv_per_tf %>%
        group_by(tfs) %>%
        dplyr::summarise(tfs,
                        tad = paste(tad,collapse ="|"),
                        adjusted_p_value = paste(adjusted_p_value,
                                                collapse = "|"),
                        motifs = paste(motifs, collapse = "|"), ) %>%
        as.data.table() %>%
        unique()

    for(i in c(1:nrow(csv_per_tf))) {

        temp <- csv_per_tf[i] %>%
            dplyr::select(tfs, motifs)
        temp <- temp %>%
            separate_rows(motifs, sep = "\\|") %>%
            as.data.table()
        temp <- unique(temp)

        temp <- temp %>%
            group_by(tfs) %>%
            dplyr::summarise(tfs,
                            motifs = paste(motifs, collapse = "|"), ) %>%
            as.data.table() %>%
            unique()

        csv_per_tf[i, 4] <- temp[1, 2]
      }

    table_tsf_motifs <- csv_per_tf %>%
        dplyr::select(tfs, motifs)

    data_visual <- left_join(data_visual, table_tsf_motifs)

    top_motifs <- top_motifs %>%
        group_by(target, id) %>%
        summarise(target, id, adjusted_p_value = mean(adjusted_p_value)) %>%
        unique()

    top_motifs <- top_motifs %>%
        group_by(target) %>%
        summarise(target,
                min_adjusted_p_value = min(adjusted_p_value),
                id, adjusted_p_value)

    who_top <- which(top_motifs$adjusted_p_value ==
                                    top_motifs$min_adjusted_p_value)

    top_motifs <- top_motifs[who_top, ]

    top_motifs <- dplyr::select(top_motifs, target, id)
    colnames(top_motifs) <- c("tfs", "top_motif")

    data_visual <- left_join(data_visual, top_motifs)
    csv_per_tf <- left_join(csv_per_tf, top_motifs)

    result <- list(table_per_tad = csv_per_tad,
                    table_per_tfs = csv_per_tf,
                    data_visual = data_visual)

    rm(csv_per_tf, csv_per_tad, data_visual, top_motifs, who_top,
        table_tsf_motifs, per_tf, temp)

    return(result)
}



#' This function is called by the `motifs_enrich` function.
#' prepare_sequences
#'
#' @param data
#' @param tech
#' @param motif_output_folder
#' @param exp_parent
#' @param dir_name
#'
#' @description
#' It is used to find the sequences, that correspond to TFBS
#' from the genomic coordinates of the events.
#'
#' @return
#'
#' @export
#'
#' @examples
prepare_sequences <- function(data,
                                tech,
                                motif_output_folder,
                                exp_parent,
                                dir_name) {

    # Filter for the sequences related to TFs location
    data <- data %>%
        dplyr::select(tad_name, start_position, end_position,
                    chromosome_name, parent, ID)

    data_cg <- data[data$parent != exp_parent, ]

    # Expand the CG sequences

    data_cg$start_position <- data_cg$start_position - 25
    data_cg$end_position <- data_cg$end_position + 25

    # Keep promoter of ENSG sequences
    if (tech == "hg19") {

        file_hg <- read.gff(paste0(dir_name, "/gencode.v19.annotation.gff3.gz"))

    } else if (tech == "hg38") {

        file_hg <- read.gff(paste0(dir_name, "/gencode.v36.annotation.gff3.gz"))

    }


    data_ensg <- data[data$parent == exp_parent, ]
    data_ensg$partID <- unlist(lapply(strsplit(data_ensg$ID,";"), '[[', 2))

    gr <- GRanges(seqnames = Rle(paste("chr", data_ensg$chromosome_name,
                                        sep = "")),
                ranges = IRanges(start = as.numeric(data_ensg$start_position),
                                    end = as.numeric(data_ensg$end_position)))

    hg <- GRanges(seqnames = Rle(file_hg$seqid),
                    ranges = IRanges(start = as.numeric(file_hg$start),
                                    end = as.numeric(file_hg$end)))

    overlaps <- findOverlaps(gr, hg)

    overlaps_from <- overlaps@from
    overlaps_to <- overlaps@to

    file_hg <- as.data.table(file_hg)
    file_hg <- file_hg[ ,c("strand", "attributes")]

    strand_data <- cbind(data_ensg[overlaps_from, 1:7],
                        file_hg[overlaps_to, ])

    strand_data <- strand_data[str_detect(strand_data$attributes,
                                            strand_data$partID), ]

    data_ensg <- strand_data[, -c("attributes")]

    data_ensg <- unique(data_ensg)

    who_plus <- which(data_ensg$strand == "+")
    who_minus <- which(data_ensg$strand == "-")

    data_ensg$end_position[who_plus] <- data_ensg$start_position[who_plus]
    data_ensg$start_position[who_plus] <-
        data_ensg$start_position[who_plus] - 2000

    data_ensg$start_position[who_minus] <- data_ensg$end_position[who_minus]
    data_ensg$end_position[who_minus] <-
        data_ensg$end_position[who_minus] + 2000

    data_ensg <- data_ensg %>%
        dplyr::select(tad_name, start_position, end_position,
                    chromosome_name, parent, ID)

    data <- rbind(data_cg, data_ensg)
    data <- data[order(data$tad_name,data$start_position, decreasing = FALSE)]
    data$AA <- c(1:nrow(data))

    keep_parent_ID <- unique(data[,c("AA", "ID", "parent")])
    keep_parent_ID <- keep_parent_ID[order(keep_parent_ID$AA,
                                            decreasing = TRUE),]

    new_data <- GRanges(seqnames = Rle(paste0(data$chromosome_name, "_",
                                                data$tad_name)),
                        ranges = IRanges(start = as.numeric(data$
                                                            start_position),
                                        end = as.numeric(data$end_position)),
                                        use.names = data$ID)

    new_data <- reduce(new_data, with.revmap = TRUE)

    new_data <- as.data.table(new_data)

    new_data <- new_data[, c("seqnames", "start" , "end" , "revmap")]
    new_data$revmap <- as.character(new_data$revmap)
    new_data$revmap <- str_replace_all(new_data$revmap, ":", " |")
    new_data$revmap <- paste0(new_data$revmap," ")
    new_data$parent <- ""

    new_data <- separate(new_data, "seqnames",
                        c("chromosome_name", "tad_name"),
                        sep = "_", remove = TRUE)

    for (k in keep_parent_ID$AA) {

        who <- which(str_detect(new_data$revmap, paste0(as.character(k), " ")))

        new_data$revmap[who] <- str_replace(new_data$revmap[who],
                                            paste0(as.character(k), " "),
                                            keep_parent_ID$ID[k])

        new_data$parent[who] <- paste0(keep_parent_ID$parent[k], "|",
                                        new_data$parent[who])
    }

    new_data$parent <- str_sub(new_data$parent,
                                start = 1, end = (nchar(new_data$parent) - 1))

    colnames(new_data) <- c("chromosome_name","tad_name",
                            "start_position", "end_position",
                            "merged_from", "parent")

    write.table(new_data,
                paste0(motif_output_folder,"/prepared sequences info.csv"),
                sep = "\t", row.names = FALSE)

    rm(data_cg, data_ensg, data, file_hg, gr, hg, overlaps,
        overlaps_from, overlaps_to, strand_data, who_minus,
        who_plus, keep_parent_ID, who)

    return(new_data)
}


#' This function is called by the `motifs_enrich` function.
#' get_dna_sequences
#'
#' @param input_data
#' @param outputs_folder
#' @param tech
#'
#' @description
#' It is used to query the Rest Ensembl API.
#' It gets the DNA sequences that correspond to the genomic
#' coordinates of the events.
#'
#' @return
#'
#' @export
#'
#' @examples
get_dna_sequences <- function(input_data,
                                outputs_folder,
                                tech) {

    if (tech == "hg19") {

        hg_version <- "GRCh37"

    } else if (tech == "hg38") {

        hg_version <- "GRCh38"
    }

    # Query Ensembl Rest Api to get the sequences
    # per TAD so as not to lose the TAD information.

    new_TADs <- input_data %>%
        group_by(tad_name)
    new_groups <- group_split(new_TADs)

    iterations <- c(1:length(new_groups))
    k <- length(new_groups) + 1
    seq_tad_number <- data.table(start = numeric(k),
                                end = numeric(k),
                                tad = character(k))

    seq_tad_number$start[1] <- 1

    for (i in iterations) {

        data <- new_groups[[i]]
        iter <- c(1:nrow(data))
        seq <- data.table(dna_seq = character(),
                          tad = character())

        for (j in iter) {

            start <- data$start_position[j]
            end <- data$end_position[j]
            chr <- data$chromosome_name[j]

            server <- "http://rest.ensembl.org"
            ext <- paste0("/sequence/region/human/", chr, ":", start, "..", end,
                            "?coord_system_version=", hg_version)

            r <- httr::GET(paste(server, ext, sep = ""),
                            content_type("text/plain"))

            seq <- rbind(seq, data.table(dna_seq = httr::content(r),
                                        tad = data$tad_name[j]))

        }

        if (!is.null(seq)) {

            write.fasta(sequences = as.list(seq$dna_seq), names = seq$tad,
                        file.out = paste0(outputs_folder,"/seq_perTADs.fasta"),
                        open = "a")

            seq_tad_number$end[i] <- (nrow(seq) + seq_tad_number$start[i] - 1)
            seq_tad_number$start[i + 1] <- (seq_tad_number$end[i] + 1)
            seq_tad_number$tad[i] <- data$tad_name[1]

        }
    }

    seq_tad_number <- seq_tad_number[which(seq_tad_number$end != 0), ]

    rm(new_groups, new_TADs, data, seq, r, iterations, iter, k)

    return(seq_tad_number)
}


