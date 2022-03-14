#' This file contains functions used in the `functional_analysis.R`
#' script.

#' This function is called by the `functional_analysis.R` script.
#' enrich_all
#'
#' @param biodata
#' @param dbs Databases used from EnrichR
#' @param cut_off cut-off Enrichr enrichment (adjusted) p-value
#' @param criterio Enrichr result column selected as criterio:
#'                `"p-value"` or `"adjusted p-value"`
#' @param type the prevously selected databases acronyms used for the names
#'             of the outputs files e.g. GO_MF for GO_Molecular_Function_2018
#' @import dplyr
#' @impot enrichR
#'
#' @description
#' It performs enrichment analysis using the Enrichr tool.
#' Enrichr input is all the genes of the dataset.
#'
#' @return
#'
#' @export
#'

enrich_all <- function(biodata,
                        dbs,
                        cut_off,
                        criterio,
                        type) {

    # Preparing data for enrichR
    data_selected <- biodata %>%
        separate_rows(Gene_id, sep = "\\|") %>%
        as.data.table()

    who_annotated <- which((data_selected$Gene_id != "NA") &
                            (data_selected$Gene_id != ""))
    data_selected <- data_selected[who_annotated, ]

    data_selected <- data_selected %>%
        dplyr::select(Gene_id,tad_name) %>%
        unique()

    data_for_enrich <- c(data_selected$Gene_id)
    data_for_enrich <- unique(data_for_enrich)

    # Enrichment with GO Molecular Function terms,
    # GO Biological Process terms and KEGG Pathways
    # using enrichr interface to connect to EnrichR.
    enriched <- enrichR::enrichr(data_for_enrich, dbs)

    result <- list()
    for (i in c(1:length(dbs))) {

        enriched_terms <- as.data.table(enriched[[dbs[i]]])

        if (criterio == "p-value") {

          enriched_terms <- subset(enriched_terms, P.value < cut_off)

        }
        else if (criterio == "adjusted p-value") {

            enriched_terms <- subset(enriched_terms,
                                    Adjusted.P.value < cut_off)
        }

        result[[i]] <- enriched_terms
    }

    result[[i+1]] <- data_selected
    names(result) <- c(type, "data_with_genes")

    rm(data_selected, data_for_enrich, who_annotated,
        enriched, enriched_terms)

    return(result)
}


#' This function is called by the `functional_analysis.R` script.
#' enrich_per_tad
#'
#' @param biodata
#' @param dbs Databases used from EnrichR
#' @param cut_off cut-off Enrichr enrichment (adjusted) p-value
#' @param criterio Enrichr result column selected as criterio:
#'                `"p-value"` or `"adjusted p-value"`
#' @param type the prevously selected databases acronyms used for the names
#'             of the outputs files e.g. GO_MF for GO_Molecular_Function_2018
#'
#' @import dplyr
#'
#' @description
#' It performs enrichment analysis using the Enrichr tool.
#' Enrichr input is the genes of the dataset grouped per TAD.
#'
#'
#' @return
#'
#' @export
#'

enrich_per_tad <- function(biodata,
                            dbs,
                            cut_off,
                            criterio,
                            type) {

    # Preparing data for enrichR
    full_tads <- biodata %>%
        dplyr::select(tad_name, Gene_id) %>%
        as.data.table()

    full_tads <- full_tads %>%
        separate_rows(Gene_id, sep = "\\|") %>%
        as.data.table()

    full_tads <- full_tads[which((full_tads$Gene_id != "NA") &
                                    (full_tads$Gene_id != "")), ]
    full_tads <- unique(full_tads)

    unique_tads <- unique(full_tads$tad_name)

    result <- list()

    for (i in c(1:length(unique_tads))) {

        data_per_tad <- full_tads[which(full_tads$tad_name == unique_tads[i]), ]
        data_for_enrich <- c(data_per_tad$Gene_id)
        data_for_enrich <- unique(data_for_enrich)

        # Enrichment with GO Molecular Function terms,
        # GO Biological Process terms abd KEGG Pathways
        # using enrichr interface to connect to EnrichR.
        enriched <- enrichr(data_for_enrich, dbs)

        for (j in c(1:length(dbs))) {

            enriched_terms <- as.data.table(enriched[[dbs[j]]])
            if (criterio == "p-value") {

                enriched_terms <- subset(enriched_terms, P.value < cut_off)

            }
            else if (criterio == "adjusted p-value") {

                enriched_terms <- subset(enriched_terms,
                                        Adjusted.P.value < cut_off)
            }

            if (i==1) {
                result[[j]] <- enriched_terms
            } else {
                result[[j]] <- rbind(result[[j]], enriched_terms)
            }

        }

    }

    result[[j+1]] <- full_tads
    names(result) <- c(type, "data_with_genes")

    rm(data_for_enrich, enriched, enriched_terms, data_per_tad,
        full_tads, unique_tads)

    return(result)

}



#' This function is called by the `functional_analysis.R` script.
#' data_analysis
#'
#' @param enriched_terms
#' @param type
#' @param data_selected
#' @param genes_coverage
#' @param p_adjust_method
#' @param min_genes
#' @import dplyr
#' @description
#' It manipulates the enriched data and performs hypergeometric test.
#'
#' @return
#'
#' @export
#'



data_analysis <- function(enriched_terms,
                            type,
                            data_selected,
                            genes_coverage,
                            p_adjust_method,
                            min_genes){

    enriched_terms <- enriched_terms %>%
      dplyr::select(Term, Overlap, Genes)

    # Calculate number of genes per term in database
    enriched_terms <- enriched_terms %>%
        separate(Overlap, c("numerator", "denominator"), sep = "\\/")
    enriched_terms$numerator <- as.numeric(as.character
                                            (enriched_terms$numerator))
    enriched_terms$denominator <- as.numeric(as.character
                                            (enriched_terms$denominator))
    enriched_terms <- enriched_terms[which
                                    (enriched_terms$numerator >= min_genes), ]

    if (nrow(enriched_terms) == 0)  return(NULL)

    enriched_terms <- enriched_terms %>%
        dplyr::select(Term, denominator, Genes) %>%
        group_by(Term) %>%
        summarise(Term,
                denominator = max(denominator),
                Genes = paste(Genes, collapse =";"), ) %>%
        as.data.table() %>%
        unique()


    data_extended <- enriched_terms %>%
        separate_rows(Genes, sep = ";", convert = TRUE) %>%
        as.data.table()

    data_extended <- left_join(data_extended,
                                data_selected,
                                by = c("Genes" = "Gene_id"))

    data_extended <- unique(data_extended)

    data_with_p <- calculate_pvalue(data_extended,
                                    data_selected,
                                    genes_coverage,
                                    p_adjust_method)

    results <- produce_outputs(data_with_p, type)

    rm(data_extended, data_selected, data_with_p, enriched_terms)

    return(results)

}



#' This function is called by the `data_analysis` function.
#' calculate_pvalue
#'
#' @param data_extended
#' @param data_selected
#' @param genes_coverage
#' @param p_adjust_method
#' @miport dplyr
#'
#' @description
#' It is used to calcualte the P value and the adjusted P value
#' for every Term per TAD.
#'
#' @return
#'
#' @export
#'

calculate_pvalue <- function(data_extended,
                            data_selected,
                            genes_coverage,
                            p_adjust_method) {

    data_extended <- data_extended[which(data_extended$tad_name != "NA"), ]
    tads <- unique(data_extended$tad_name)

    data_with_p <- data.table(Term = character(),
                              TAD = character(),
                              P.value = numeric())

    for (i in c(1:length(tads))) {

        tad <- tads[i]
        tad_terms <- data_extended[which(data_extended$tad_name == tad), ]
        terms_number <- dplyr::count(tad_terms, Term)
        iter <- c(1:nrow(terms_number))

        for (j in iter) {

            hit_in_sample <- terms_number$n[j]
            who <- which(tad_terms$Term == as.character(terms_number[j, 1]))
            hit_in_pop <- tad_terms[who, denominator]
            fail_in_pop <- genes_coverage - hit_in_pop
            tad_genes <- data_selected[which(data_selected$tad_name == tad), ]
            tad_genes <- tad_genes[which(tad_genes$Gene_id != ""), ]
            sample_size <- length(tad_genes$Gene_id)

            # Test for over-representation, enrichment
            p_value <- phyper(hit_in_sample - 1,
                                hit_in_pop,
                                fail_in_pop,
                                sample_size,
                                lower.tail = FALSE)

            data_with_p <- rbind(data_with_p,
                                data.table(Term =
                                    as.character(terms_number[j, 1]),
                                    TAD = as.character(tad),
                                    P.value = as.numeric(p_value[1])))

        }
    }

    data_with_p <- data_with_p[which(data_with_p$P.value != "NA"), ]

    # Adjust p values
    data_with_p$P.adjust <- p.adjust(data_with_p$P.value,
                                method = p_adjust_method)

    rm(data_extended, tads, tad, tad_terms, terms_number, iter,
        hit_in_sample, who, hit_in_pop, fail_in_pop, tad_genes,
        sample_size, p_value)

    return(data_with_p)
}


#' This function is called by the `data_analysis` function.
#' produce_outputs
#'
#' @param data_with_p
#' @param type
#'
#' @description
#' It manipulates the enriched data after the analysis and creates
#' three output tables to be used for the Output csv files and the
#' visualization.
#'
#' @return
#'
#' @export
#'
#'

produce_outputs <- function(data_with_p,
                            type) {

    if (str_detect(type,"GO")) {

        data_with_p$Term <- str_remove(data_with_p$Term, "\\)")
        data_with_p <- data_with_p %>%
            separate(Term, c("go_term", "go_number"), "\\(GO:" ) %>%
            as.data.table()

        data_visual <- data_with_p %>%
            dplyr::select(TAD, go_term, go_number, P.value, P.adjust)

        data_with_p <- data_visual %>%
            dplyr::group_by(TAD) %>%
            dplyr::summarise(TAD,
                            go_term = paste(go_term, collapse = "|"),
                            go_number = paste(go_number, collapse = "|"),
                            p_value = paste(P.value, collapse = "|"),
                            p_adjust = paste(P.adjust, collapse = "|"), ) %>%
            as.data.table() %>%
            unique()

        data_per_term <- data_visual %>%
            dplyr::group_by(go_term, go_number) %>%
            dplyr::summarise(go_term, go_number,
                            TAD = paste(TAD, collapse = "|"),
                            p_value = paste(P.value, collapse = "|"),
                            p_adjust = paste(P.adjust, collapse = "|"), ) %>%
            as.data.table() %>%
            unique()

        column_term <- paste0(type, "_Term")
        column_id   <- paste0(type, "_number")
        column_p    <- paste0(type, "_p_value")
        column_adj  <- paste0(type, "_p_adjust")
        colnames(data_with_p)   <- c("TAD", column_term, column_id, column_p,
                                    column_adj)
        colnames(data_visual)   <- c("TAD", "Term", "ID", "p_value", "p_adjust")
        colnames(data_per_term) <- c("go_term", "go_id", "TAD", "p_value",
                                    "p_adjust")

    }
    else{

        data_visual <- data_with_p %>%
            dplyr::select(TAD, Term, P.value, P.adjust)

        data_with_p <- data_visual %>%
            dplyr::group_by(TAD) %>%
            dplyr::summarise(TAD, Term = paste(Term, collapse = "|"),
                            p_value = paste(P.value, collapse = "|"),
                            p_adjust = paste(P.adjust, collapse = "|"), ) %>%
            as.data.table() %>%
            unique()

        data_per_term <- data_visual %>%
            dplyr::group_by(Term) %>%
            dplyr::summarise(Term,
                            TAD = paste(TAD, collapse = "|"),
                            p_value = paste(P.value, collapse = "|"),
                            p_adjust = paste(P.adjust, collapse = "|"), ) %>%
            as.data.table() %>%
            unique()

        colnames(data_with_p)   <- c("TAD", "Term", "p_value", "p_adjust")
        colnames(data_visual)   <- c("TAD", "Term", "p_value", "p_adjust")
        colnames(data_per_term) <- c("Term", "TAD", "p_value", "p_adjust")

    }


    new_list <- list(data_visual = data_visual,
                     data_per_term = data_per_term,
                     data_per_tad = data_with_p)

    rm(data_with_p, data_visual, data_per_term)

    return(new_list)
}


