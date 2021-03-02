#This file contains functions used in the "enrichmentAnalysis.R" script


#This function is called by the "dataAnalysis" function
#It is used to calcualte the P value and the adjusted P value for every Term per TAD
calculatePvalue <- function(data.extended ,data_selected ,genes.coverage ,p.adjust.method ) {
  
  data.extended <- data.extended[which(data.extended$tad_name != "NA"),]
  tads <- unique(data.extended$tad_name)
  
  data.with.p <- data.table(Term = character(),
                            TAD = character(),
                            P.value = numeric(),
                            stringsAsFactors = F)
  
  iterations <- c(1:length(tads))
  
  for (i in iterations){
    
    tad <- tads[i]
    tad.terms <- data.extended[which(data.extended$tad_name == tad),]
    terms.number <- dplyr::count(tad.terms,Term)
    iter <- c(1:nrow(terms.number))
    
    for (j in iter){
      
      hitInSample <- terms.number$n[j]
      hitInPop <- tad.terms[which(tad.terms$Term == as.character(terms.number[j,1])),denominator]
      failInPop <- genes.coverage -  hitInPop
      tad.genes <- data_selected[which(data_selected$tad_name == tad),]
      tad.genes <- tad.genes[which(tad.genes$Gene_id != ""),]
      sampleSize <- length(tad.genes$Gene_id)
      
      #test for over-representation, enrichment
      p.value <- phyper(hitInSample-1,hitInPop, failInPop,sampleSize, lower.tail=FALSE)
      
      temp <- data.table(Term = as.character(terms.number[j,1]),
                         TAD = as.character(tad),
                         P.value = as.numeric(p.value[1]),
                         stringsAsFactors = F)
      
      data.with.p <- rbind(data.with.p ,temp)
      
    }
  }
  
  data.with.p <- data.with.p[which(data.with.p$P.value != "NA"),]
  #data.with.p <- unique(data.with.p)
  
  #p adjust
  data.with.p$P.adjust <- p.adjust(data.with.p$P.value, method = p.adjust.method)
  
  return(data.with.p) 
}


#This function is called by the "enrichmentAnalysis.R" script
#It performs enrichment analysis using the Enrichr tool
#Enrichr input is all the genes of the dataset
enrichAll <- function(biodata, dbs, threshold, criterio, type){
  
  #preparing data for enrichR
  data_selected <- biodata %>% 
    separate_rows(Gene_id, sep = "\\|") %>%
    as.data.table()
  data_selected <- data_selected[which((data_selected$Gene_id != "NA") & (data_selected$Gene_id != "")), ]
  data_selected <- data_selected %>%
    dplyr::select(Gene_id,tad_name) %>%
    unique()
  dataForEnrich <-c(data_selected$Gene_id)
  dataForEnrich <- unique(dataForEnrich) 
  
  #Enrichment with GO Molecular Function terms, GO Biological Process terms and KEGG Pathways
  #using enrichr interface to connect to EnrichR
  enriched <- enrichR::enrichr(dataForEnrich, dbs)
  
  result <- list()
  for (i in c(1:length(dbs))){
    
    enriched_terms <- as.data.table(enriched[[dbs[i]]])
    if (criterio == "P.value"){
      enriched_terms <- subset(enriched_terms, P.value < threshold)
    }
    
    if (criterio == "Adjusted.P.value"){
      enriched_terms <- subset(enriched_terms, Adjusted.P.value < threshold)
    }
    
    result[[i]] <- enriched_terms
  }
  
  result[[i+1]] <- data_selected 
  names(result) <- c(type,"data.with.genes")
  
   return(result)
}


#This function is called by the "enrichmentAnalysis.R" script
#It performs enrichment analysis using the Enrichr tool
#Enrichr input is the genes of the dataset grouped per TAD
enrichPerTAD <- function(biodata,dbs, threshold, criterio, type){
  
  #preparing data for enrichR
  full.tads <- biodata %>%
    dplyr::select(tad_name,Gene_id) %>%
    as.data.table()
  
  data_ <- full.tads %>% 
    separate_rows(Gene_id, sep = "\\|") %>%
    as.data.table()
  data_ <- data_[which((data_$Gene_id != "NA") & (data_$Gene_id != "")), ]
  data_ <- unique(data_)
  
  unique.tads <- unique(data_$tad_name)

  #iterations <- c(1:10)
  iterations <- c(1:length(unique.tads))
  
  result <- list()
  for (i in iterations){
    
    data.per.tad <- data_[which(data_$tad_name == unique.tads[i]),]
    dataForEnrich <-c(data.per.tad$Gene_id)
    dataForEnrich <- unique(dataForEnrich)
    
    #Enrichment with GO Molecular Function terms, GO Biological Process terms abd KEGG Pathways
    #using enrichr interface to connect to EnrichR
    enriched <- enrichr(dataForEnrich, dbs)
 
    
    for (j in c(1:length(dbs))){
          
      enriched_terms <- as.data.table(enriched[[dbs[j]]])
      if (criterio == "P.value"){
          enriched_terms <- subset(enriched_terms, P.value < threshold)
      }
          
      if (criterio == "Adjusted.P.value"){
          enriched_terms <- subset(enriched_terms, Adjusted.P.value < threshold)
      }
      if (i==1) {result[[j]] <- enriched_terms
      
      }else{ result[[j]] <- rbind(result[[j]],enriched_terms)}    
      
    }
        
  }

  result[[j+1]] <- data_
  names(result) <- c(type,"data.with.genes")
  
  return(result)
  
}


#This function is called by the "dataAnalysis" function
#It manipulates the enriched data after the analysis and creates three output data.tables 
#to be used for the Output csv files and the visualization
produceOutputs <- function(data.with.p,type){
  if ( str_detect(type,"GO")){
    
    data.with.p$Term <- str_remove(data.with.p$Term, "\\)")
    data.with.p <- data.with.p %>% 
      separate(Term, c("GO.term", "GO.number"), "\\(GO:" ) %>%
      as.data.table()
    
    data.visual <- data.with.p%>%
      dplyr::select(TAD, GO.term, GO.number, P.value, P.adjust)
    
    data.with.p<- data.visual %>%
      dplyr::group_by(TAD) %>%
      dplyr::summarise(TAD,GO.term = paste(GO.term, collapse = "|"),
                       GO.number = paste(GO.number, collapse = "|"),
                       p.value = paste(P.value, collapse = "|"),
                       P.adjust = paste(P.adjust, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    
    data.per.term <- data.visual %>%
      dplyr::group_by(GO.term, GO.number) %>%
      dplyr::summarise(GO.term, GO.number, 
                       TAD = paste(TAD, collapse = "|"),
                       p.value = paste(P.value, collapse = "|"),
                       P.adjust = paste(P.adjust, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    
    column.term <- paste0(type,".Term")
    column.id <- paste0(type,".number")
    column.p <- paste0(type,".P.value")
    column.adj <- paste0(type, ".P.adjust")
    colnames(data.with.p) <- c("TAD",column.term,column.id,column.p,column.adj)
    colnames(data.visual) <- c("TAD","Term","ID","P.value","P.adjust")
    colnames(data.per.term) <- c("GO.Term","GO.ID","TAD","P.value","P.adjust")
    
  }else{
    
    data.visual <- data.with.p%>%
      dplyr::select(TAD, Term, P.value, P.adjust)
    
    data.with.p<- data.visual %>%
      dplyr::group_by(TAD) %>%
      dplyr::summarise(TAD,Term = paste(Term, collapse = "|"),
                       p.value = paste(P.value, collapse = "|"),
                       P.adjust = paste(P.adjust, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    
    data.per.term <- data.visual %>%
      dplyr::group_by(Term) %>%
      dplyr::summarise(Term , 
                       TAD = paste(TAD, collapse = "|"),
                       p.value = paste(P.value, collapse = "|"),
                       P.adjust = paste(P.adjust, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    
    colnames(data.with.p) <- c("TAD","Term","P.value","P.adjust")
    colnames(data.visual) <- c("TAD","Term","P.value","P.adjust")
    colnames(data.per.term) <- c("Term","TAD","P.value","P.adjust")
    
  }
  
  
  newList <- list(data.visual = data.visual, data.perTerm = data.per.term, data.perTAD = data.with.p)
  return(newList)
}


#This function is called by the "enrichmentAnalysis.R" script
#It manipulates the enriched data and performs hypergeometric test
dataAnalysis <- function(enriched_terms, type,data_selected, genes.coverage, p.adjust.method, min_genes){
      
      enriched_terms <- enriched_terms %>%
        dplyr::select(Term,Overlap,Genes)
      
      #calculate number of genes per term in database 
      enriched_terms <- enriched_terms %>%
        separate(Overlap, c("numerator", "denominator"), sep = "\\/")
      enriched_terms$numerator <- as.numeric(as.character(enriched_terms$numerator))
      enriched_terms$denominator <- as.numeric(as.character(enriched_terms$denominator))
      enriched_terms <- enriched_terms[which(enriched_terms$numerator >= min_genes),]
      if (nrow(enriched_terms) == 0){
        return()
      }
      
      enriched_terms <- enriched_terms %>%
          dplyr::select(Term, denominator, Genes)%>%
          group_by(Term) %>%
          summarise(Term, 
                    denominator = max(denominator),
                    Genes = paste(Genes, collapse =";"),) %>%
          as.data.table() %>%
          unique()
      
      
      data.extended <- enriched_terms %>%
        separate_rows(Genes, sep = ";", convert = TRUE) %>%
        as.data.table()
      
      data.extended <- left_join(data.extended,data_selected, by =c("Genes" = "Gene_id"))
      
      data.extended <- unique(data.extended) 
      
      data.with.p <- calculatePvalue(data.extended,data_selected, genes.coverage, p.adjust.method)
      
      resultsList <- produceOutputs(data.with.p,type)
      
      return(resultsList)
  
}



#This function is called by the "enrichmentAnalysis.R" script
#It is used to query the KEGG Pathway DB about the Kegg ids of the pathways returned from the enrichment analysis
getKEGGIds <- function(enriched_KEGG, genes_diff){
  
  enriched_KEGG <- enriched_KEGG %>%
    dplyr::select(Term, Genes)
  
  enriched_KEGG <- enriched_KEGG %>% 
    group_by(Term) %>%
    summarise(Term,
              Genes = paste(Genes, collapse = ";"),) %>%
    as.data.table()
  
  enriched_KEGG <- unique(enriched_KEGG)
  
  enriched_KEGG$ID <- "NA"
  iter <- c(1:nrow(enriched_KEGG))
  for (i in iter){
    if (str_detect(enriched_KEGG$Term[i],"\\(")){
      split_term <- str_split(enriched_KEGG$Term[i],"\\(",simplify = T)
      enriched_KEGG$Term[i] <- split_term[1]
    }
    enriched_KEGG$Term[i] <- str_replace(enriched_KEGG$Term[i],"/","")
    k <- keggFind("pathway",c(enriched_KEGG$Term[i]))
    if (!is_empty(k)){
      enriched_KEGG$ID[i] <- str_replace(names(k[1]),"path:map","hsa")
    }
  }
  
  enriched_KEGG <- enriched_KEGG[which(enriched_KEGG$ID != "NA"),]
  
  genes_diff <- genes_diff %>% 
    separate_rows(Gene_id, sep = "\\|") %>%
    as.data.table()
  genes_diff <- genes_diff[which((genes_diff$Gene_id != "NA") &(genes_diff$Gene_id !="")), ]
  
  output.csv <- data.table(Term = character(),
                           ID = character(),
                           Genes = character(),
                           diff = character())
  pathview.input <- data.table(Term = character(),
                               ID = character(),
                               Genes = character(),
                               diff = character())
  
  loops <- c(1:nrow(enriched_KEGG))
  for (l in loops){
    
    path <- enriched_KEGG[l]
    path <- path %>%
      separate_rows(Genes, sep = ";", convert = TRUE) %>%
      unique() %>%
      as.data.table()
    path <- merge(path, genes_diff,by.x = "Genes", by.y = "Gene_id")
    path <- group_by(path,Genes) %>%
      dplyr::summarise(ID,Term,Genes,diff = mean(diff)) %>%
      as.data.table() %>%
      unique()
    
    output.path <- path %>%
      group_by(Term, ID) %>%
      summarise(Term, ID, Genes = paste0(Genes, collapse = "|"),diff = paste0(diff, collapse ="|")) %>%
      unique()
    
    output.csv <- rbind(output.csv, output.path)
    
    pathview.input <- rbind(pathview.input, path)
    
  }
  
  newList <- list(pathview.input = pathview.input, output.csv = output.csv )
  return(newList)
}


#This function is called by the "enrichmentAnalysis.R" script
#It creates the output folders
createFolders <- function(output_folder){

  go_output_folder = paste0(output_folder,"/GO EA Outputs")
  kegg_output_folder = paste0(output_folder,"/Pathways EA Outputs")
  go_all_outputs = paste0(go_output_folder,"/GO enrich all")
  go_all_images = paste0(go_all_outputs,"/Images")
  kegg_all_outputs = paste0(kegg_output_folder,"/Path enrich all")
  kegg_all_images = paste0(kegg_all_outputs,"/Images")
  motif_output_folder = paste0(output_folder,"/Motif EA Outputs")
  image_output_folder = paste0(motif_output_folder,"/Images")
  
  dir.create(output_folder, showWarnings = FALSE)
  dir.create(go_output_folder, showWarnings = FALSE)
  dir.create(go_all_outputs, showWarnings = FALSE)
  dir.create(go_all_images, showWarnings = FALSE)
  dir.create(kegg_output_folder, showWarnings = FALSE)
  dir.create(kegg_all_outputs, showWarnings = FALSE)
  dir.create(kegg_all_images, showWarnings = FALSE)
  dir.create(motif_output_folder, showWarnings = FALSE)
  dir.create(image_output_folder, showWarnings = FALSE)
  
  folder <- list(goAllOutputs = go_all_outputs,
                 goAllImages = go_all_images,
                 keggAllOutputs = kegg_all_outputs,
                 keggAllImages = kegg_all_images,
                 motifOutputsFolder = motif_output_folder,
                 motifImageOutputs = image_output_folder)  
  
  if (str_detect(output_folder,"TADiff")){
    go_per_outputs = paste0(go_output_folder,"/GO enrich per TAD")
    go_per_images = paste0(go_per_outputs,"/Images")
    kegg_per_outputs = paste0(kegg_output_folder,"/Path enrich per TAD")
    kegg_per_images = paste0(kegg_per_outputs,"/Images")
    
    dir.create(go_per_outputs, showWarnings = FALSE)
    dir.create(go_per_images, showWarnings = FALSE)
    dir.create(kegg_per_outputs, showWarnings = FALSE)
    dir.create(kegg_per_images, showWarnings = FALSE)
    
    folder <- append(folder, list(goPerOutputs = go_per_outputs,
                                  goPerImages = go_per_images,
                                  keggPerOutputs = kegg_per_outputs,
                                  keggPerImages = kegg_per_images))
  }
 
  return(folder)
}

