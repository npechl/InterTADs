#' This file contains functions used in the `03_functional_analysis.R`
#' script.


#' This function is called by the `03_functional_analysis.R` script.
#' 
#' @description 
#' It creates the output folders.
create_folders <- function(output_folder) {
  
  go_output_folder <- paste0(output_folder, "/GO EA Outputs")
  kegg_output_folder <- paste0(output_folder, "/Pathways EA Outputs")
  go_all_outputs <- paste0(go_output_folder, "/GO enrich all")
  go_all_images <- paste0(go_all_outputs, "/Images")
  kegg_all_outputs <- paste0(kegg_output_folder, "/Path enrich all")
  kegg_all_images <- paste0(kegg_all_outputs, "/Images")
  motif_output_folder <- paste0(output_folder, "/Motif EA Outputs")
  image_output_folder <- paste0(motif_output_folder, "/Images")
  
  dir.create(output_folder, showWarnings = FALSE)
  dir.create(go_output_folder, showWarnings = FALSE)
  dir.create(go_all_outputs, showWarnings = FALSE)
  dir.create(go_all_images, showWarnings = FALSE)
  dir.create(kegg_output_folder, showWarnings = FALSE)
  dir.create(kegg_all_outputs, showWarnings = FALSE)
  dir.create(kegg_all_images, showWarnings = FALSE)
  dir.create(motif_output_folder, showWarnings = FALSE)
  dir.create(image_output_folder, showWarnings = FALSE)
  
  folder <- list(go_all_outputs = go_all_outputs,
                 go_all_images = go_all_images,
                 kegg_all_outputs = kegg_all_outputs,
                 kegg_all_images = kegg_all_images,
                 motif_outputs = motif_output_folder,
                 motif_images = image_output_folder)  
  
  if (str_detect(output_folder, "TADiff")) {
    
    go_per_outputs <- paste0(go_output_folder, "/GO enrich per TAD")
    go_per_images <- paste0(go_per_outputs, "/Images")
    kegg_per_outputs <- paste0(kegg_output_folder, "/Path enrich per TAD")
    kegg_per_images <- paste0(kegg_per_outputs, "/Images")
    
    dir.create(go_per_outputs, showWarnings = FALSE)
    dir.create(go_per_images, showWarnings = FALSE)
    dir.create(kegg_per_outputs, showWarnings = FALSE)
    dir.create(kegg_per_images, showWarnings = FALSE)
    
    folder <- append(folder, list(go_per_outputs = go_per_outputs,
                                  go_per_images = go_per_images,
                                  kegg_per_outputs = kegg_per_outputs,
                                  kegg_per_images = kegg_per_images))
  }
  
  rm(go_output_folder, kegg_output_folder, go_all_outputs, go_all_images,
     kegg_all_outputs, kegg_all_images, motif_output_folder, image_output_folder,
     go_per_outputs, go_per_images, kegg_per_outputs, kegg_per_images)
  
  return(folder)
}




#' This function is called by the `03_functional_analysis.R` script.
#' 
#' @description 
#' It performs enrichment analysis using the Enrichr tool.
#' Enrichr input is all the genes of the dataset.
enrich_all <- function(biodata, 
                       dbs,
                       threshold, 
                       criterio, 
                       type) {
  
  # Preparing data for enrichR
  data_selected <- biodata %>% 
    separate_rows(Gene_id, sep = "\\|") %>%
    as.data.table()
  
  who_annotated <- which((data_selected$Gene_id != "NA") & (data_selected$Gene_id != ""))
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
      
      enriched_terms <- subset(enriched_terms, P.value < threshold)
      
    } else if (criterio == "adjusted p-value") {
      
      enriched_terms <- subset(enriched_terms, Adjusted.P.value < threshold)
    }
    
    result[[i]] <- enriched_terms
  }
  
  result[[i+1]] <- data_selected 
  names(result) <- c(type, "data_with_genes")
  
  rm(data_selected, data_for_enrich, who_annotated, 
     enriched, enriched_terms)
  
  return(result)
}


#' This function is called by the `03_functional_analysis.R` script.
#' 
#' @description 
#' It performs enrichment analysis using the Enrichr tool.
#' Enrichr input is the genes of the dataset grouped per TAD.
enrich_per_tad <- function(biodata, 
                           dbs,
                           threshold, 
                           criterio,
                           type) {
  
  # Preparing data for enrichR
  full_tads <- biodata %>%
    dplyr::select(tad_name, Gene_id) %>%
    as.data.table()
  
  full_tads <- full_tads %>% 
    separate_rows(Gene_id, sep = "\\|") %>%
    as.data.table() 
  
  full_tads <- full_tads[which((full_tads$Gene_id != "NA") & (full_tads$Gene_id != "")), ]
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
        
        enriched_terms <- subset(enriched_terms, P.value < threshold)
        
      } else if (criterio == "adjusted p-value") {
        
        enriched_terms <- subset(enriched_terms, Adjusted.P.value < threshold)
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



#' This function is called by the `03_functional_analysis.R` script.
#' 
#' @description 
#' It manipulates the enriched data and performs hypergeometric test.
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
  enriched_terms$numerator <- as.numeric(as.character(enriched_terms$numerator))
  enriched_terms$denominator <- as.numeric(as.character(enriched_terms$denominator))
  enriched_terms <- enriched_terms[which(enriched_terms$numerator >= min_genes), ]
  
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
#' 
#' @description 
#' It is used to calcualte the P value and the adjusted P value 
#' for every Term per TAD.
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
      
      data_with_p <- rbind(data_with_p, data.table(Term = as.character(terms_number[j, 1]),
                                                   TAD = as.character(tad),
                                                   P.value = as.numeric(p_value[1])))
      
    }
  }
  
  data_with_p <- data_with_p[which(data_with_p$P.value != "NA"), ]
  
  # Adjust p values
  data_with_p$P.adjust <- p.adjust(data_with_p$P.value, method = p_adjust_method)
  
  rm(data_extended, tads, tad, tad_terms, terms_number, iter, 
     hit_in_sample, who, hit_in_pop, fail_in_pop, tad_genes, 
     sample_size, p_value)
  
  return(data_with_p) 
}


#' This function is called by the `data_analysis` function.
#' 
#' @description 
#' It manipulates the enriched data after the analysis and creates 
#' three output tables to be used for the Output csv files and the 
#' visualization.
produce_outputs <- function(data_with_p,
                            type) {
  
  if (str_detect(type,"GO")) {
    
    data_with_p$Term <- str_remove(data_with_p$Term, "\\)")
    data_with_p <- data_with_p %>% 
      separate(Term, c("go_term", "go_number"), "\\(GO:" ) %>%
      as.data.table()
    
    data_visual <- data_with_p %>%
      dplyr::select(TAD, go_term, go_number, P.value, P.adjust)
    
    data_with_p<- data_visual %>%
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
    colnames(data_with_p)   <- c("TAD", column_term, column_id, column_p, column_adj)
    colnames(data_visual)   <- c("TAD", "Term", "ID", "p_value", "p_adjust")
    colnames(data_per_term) <- c("go_term", "go_id", "TAD", "p_value", "p_adjust")
    
  }else{
    
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


# 
# 
# # This function is called by the "enrichmentAnalysis.R" script
# # It is used to query the KEGG Pathway DB about the Kegg ids of the pathways returned from the enrichment analysis
# getKEGGIds <- function(enriched_KEGG, genes_diff){
#   
#   enriched_KEGG <- enriched_KEGG %>%
#     dplyr::select(Term, Genes)
#   
#   enriched_KEGG <- enriched_KEGG %>% 
#     group_by(Term) %>%
#     summarise(Term,
#               Genes = paste(Genes, collapse = ";"),) %>%
#     as.data.table()
#   
#   enriched_KEGG <- unique(enriched_KEGG)
#   
#   enriched_KEGG$ID <- "NA"
#   iter <- c(1:nrow(enriched_KEGG))
#   for (i in iter){
#     if (str_detect(enriched_KEGG$Term[i],"\\(")){
#       split_term <- str_split(enriched_KEGG$Term[i],"\\(",simplify = T)
#       enriched_KEGG$Term[i] <- split_term[1]
#     }
#     enriched_KEGG$Term[i] <- str_replace(enriched_KEGG$Term[i],"/","")
#     k <- keggFind("pathway",c(enriched_KEGG$Term[i]))
#     if (!is_empty(k)){
#       enriched_KEGG$ID[i] <- str_replace(names(k[1]),"path:map","hsa")
#     }
#   }
#   
#   enriched_KEGG <- enriched_KEGG[which(enriched_KEGG$ID != "NA"),]
#   
#   genes_diff <- genes_diff %>% 
#     separate_rows(Gene_id, sep = "\\|") %>%
#     as.data.table()
#   genes_diff <- genes_diff[which((genes_diff$Gene_id != "NA") &(genes_diff$Gene_id !="")), ]
#   
#   output.csv <- data.table(Term = character(),
#                            ID = character(),
#                            Genes = character(),
#                            diff = character())
#   pathview.input <- data.table(Term = character(),
#                                ID = character(),
#                                Genes = character(),
#                                diff = character())
#   
#   loops <- c(1:nrow(enriched_KEGG))
#   for (l in loops){
#     
#     path <- enriched_KEGG[l]
#     path <- path %>%
#       separate_rows(Genes, sep = ";", convert = TRUE) %>%
#       unique() %>%
#       as.data.table()
#     path <- merge(path, genes_diff,by.x = "Genes", by.y = "Gene_id")
#     path <- group_by(path,Genes) %>%
#       dplyr::summarise(ID,Term,Genes,diff = mean(diff)) %>%
#       as.data.table() %>%
#       unique()
#     
#     output.path <- path %>%
#       group_by(Term, ID) %>%
#       summarise(Term, ID, Genes = paste0(Genes, collapse = "|"),diff = paste0(diff, collapse ="|")) %>%
#       unique()
#     
#     output.csv <- rbind(output.csv, output.path)
#     
#     pathview.input <- rbind(pathview.input, path)
#     
#   }
#   
#   newlist <- list(pathview.input = pathview.input, output.csv = output.csv )
#   return(newlist)
# }
# 
