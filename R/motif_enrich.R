#This file contains functions used in the "enrichmentAnalysis.R" script


#This function is called by the "motifEnrich" function
#It is used to find the sequences, that correspond to TFBS from the genomic coordinates of the events
prepareSequences <- function(biodata, tech, motif_output_folder, exp.parent, dir_name){
  
  #filter for the sequences related to TFs location
  data <- biodata %>%
    dplyr::select(ID,tad_name,chromosome_name,start_position,end_position, parent)
  
  #expression data parent == 1
  
  data.cg <- data[data$parent != exp.parent,] 
  
  #expand the CG sequences
  
  data.cg$start_position <- data.cg$start_position - 25
  data.cg$end_position <- data.cg$end_position + 25
  
  
  data.cg <- data.cg %>% dplyr::select(tad_name, start_position, end_position, chromosome_name, parent, ID)
  
  #keep promoter of ENSG sequences
  if (tech == "hg19"){
    
    download.file(url="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz",
                  destfile=paste0(dir_name,'/gencode.v19.annotation.gff3.gz'), method='curl')
    file.hg <- read.gff(paste0(dir_name,"/gencode.v19.annotation.gff3.gz"))
    
  }else if (tech == "hg38"){
    
    download.file(url="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gff3.gz",
                  destfile=paste0(dir_name,'/gencode.v36.annotation.gff3.gz'), method='curl')
    file.hg <- read.gff(paste0(dir_name,"/gencode.v36.annotation.gff3.gz"))
    
  }
  
  
  data.ensg <- data[data$parent == exp.parent,] 
  data.ensg$partID <- unlist(lapply(strsplit(data.ensg$ID,";"),'[[',2))
  
  gr = GRanges(seqnames = Rle(paste("chr", data.ensg$chromosome_name, sep = "")),
               ranges = IRanges(start = as.numeric(data.ensg$start_position),
                                end = as.numeric(data.ensg$end_position)))
  
  hg = GRanges(seqnames = Rle(file.hg$seqid),
                 ranges = IRanges(start = as.numeric(file.hg$start),
                                  end = as.numeric(file.hg$end)))
  
  overlaps = findOverlaps(gr, hg)
  
  overlaps.from = overlaps@from
  overlaps.to = overlaps@to
  
  file.hg = as.data.table(file.hg)
  file.hg = file.hg[,c("strand", "attributes")]
  
  strand.data = cbind(data.ensg[overlaps.from,1:7], file.hg[overlaps.to, ])
  strand.data <- strand.data[str_detect(strand.data$attributes, strand.data$partID),]

  data.ensg = strand.data[,-c("attributes")]
  
  data.ensg <- unique(data.ensg)
  
  who.plus <- which(data.ensg$strand == "+")
  who.minus <- which(data.ensg$strand == "-")
  
  data.ensg$end_position[who.plus] <- data.ensg$start_position[who.plus]
  data.ensg$start_position[who.plus] <- data.ensg$start_position[who.plus] - 2000
  
  data.ensg$start_position[who.minus] <- data.ensg$end_position[who.minus] 
  data.ensg$end_position[who.minus] <- data.ensg$end_position[who.minus] + 2000
  
  data.ensg <- data.ensg %>% dplyr::select(tad_name, start_position, end_position, chromosome_name, parent, ID)
  
  data <- rbind(data.cg, data.ensg)
  data <- data[order(data$tad_name,data$start_position,decreasing = F)]
  data$AA <- c(1:nrow(data))
  
  keep.parent.ID <- unique(data[,c("AA", "ID","parent")])
  keep.parent.ID <- keep.parent.ID[order(keep.parent.ID$AA, decreasing = T),]
  
  new.data = GRanges(seqnames = Rle(paste0(data$chromosome_name,"_",data$tad_name)),
               ranges = IRanges(start = as.numeric(data$start_position),
                                end = as.numeric(data$end_position)), use.names = data$ID)
  
  new.data <- reduce(new.data, with.revmap = T)
  
  new.data <- as.data.table(new.data)
  
  new.data <- new.data[,c("seqnames", "start" , "end" ,"revmap")]
  new.data$revmap <- as.character(new.data$revmap)
  new.data$revmap <- str_replace_all(new.data$revmap,":"," |")
  new.data$revmap <- paste0(new.data$revmap," ")
  new.data$parent <- ""
  
  new.data <- separate(new.data,"seqnames",c("chromosome_name", "tad_name"),sep = "_", remove = T)
  
  for (k in keep.parent.ID$AA){
    
    who <- which(str_detect(new.data$revmap,paste0(as.character(k)," ")))
    new.data$revmap[who] <- str_replace(new.data$revmap[who],paste0(as.character(k)," "),keep.parent.ID$ID[k])
    new.data$parent[who] <- paste0(keep.parent.ID$parent[k],"|",new.data$parent[who]) 
  } 
  
  new.data$parent <- str_sub(new.data$parent, start = 1, end = (nchar(new.data$parent)-1) )
  
  colnames(new.data) <- c("chromosome_name","tad_name","start_position", "end_position", "merged.from", "parent")
  
  write.table(new.data, paste0(motif_output_folder,"/prepared sequences info.csv"), sep = "\t")
  return(new.data)
}


#This function is called by the "motifEnrich" function
#It is used to query the Rest Ensembl API  
#It gets the DNA sequences that correspond to the genomic coordinates of the events
getDNASequences <- function(input.data, outputs_folder, tech){
  
  if (tech == "hg19"){
    hg.version = "GRCh37"
  }else if (tech == "hg38"){
    hg.version = "GRCh38"
  }
  
  #query Ensembl Rest Api to get the sequences 
  #per TAD so as not to lose the TAD information
  new.TADs <- input.data %>%
    group_by(tad_name)
  new.groups <- group_split(new.TADs)
  
  #iterations <- c(1:5)
  #k <- 5+1
  iterations <- c(1:length(new.groups))
  k <- length(new.groups)+1
  seq.tad.number <- data.table(start = numeric(k),
                               end = numeric(k),
                               tad = character(k),
                               stringsAsFactors = F)
  seq.tad.number$start[1] <- 1
  
  for (i in iterations){
    
    data <- new.groups[[i]]
    iter <- c(1:nrow(data))  
    seq <- data.table(dna.seq = character(),
                      tad = character(),
                      stringsAsFactors = FALSE)
    
    for (j in iter){
      start <- data$start_position[j]
      end <- data$end_position[j]
      chr <- data$chromosome_name[j]
      server <- "http://rest.ensembl.org"
      ext <- paste0("/sequence/region/human/",chr,":",start, "..", end, "?coord_system_version=",hg.version )
      r <- httr::GET(paste(server, ext, sep = ""), content_type("text/plain"))
      
     # if (httr::status_code(r) == 503){
    #    for (a in c(1:5)){
     #     r <- httr::GET(paste(server, ext, sep = ""), content_type("text/plain"))
      #  }
    #  }
      
      #stop_for_status(r)
      x <- data.table(dna.seq = httr::content(r),
                      tad = data$tad_name[j], 
                      stringsAsFactors = FALSE)
      seq <- rbind(seq,x)
    } 
    
    if (!is.null(seq)){
      
      write.fasta(sequences = as.list(seq$dna.seq) , names = seq$tad, file.out = paste0(outputs_folder,"/seq_perTADs.fasta"), open = "a")
      seq.tad.number$end[i] <- (nrow(seq) + seq.tad.number$start[i] -1)
      seq.tad.number$start[i+1] <- (seq.tad.number$end[i] +1)
      seq.tad.number$tad[i] <- data$tad_name[1]
      
    }
  }
  
  seq.tad.number <- seq.tad.number[which(seq.tad.number$end != 0),]
  return(seq.tad.number)
}


#This function is called by the "enrichmentAnalysis.R" script 
#It manipulates the enriched data after the analysis and creates three output data.tables 
#to be used for the Output csv files and the visualization
motifOutputs <- function(report.list){
  
  csv_perTAD <- data.table(TAD = character(),
                           TFs = character(),
                           Motifs = character(),
                           Adjusted.P.value = character(),
                           P.value = character())
  
  csv_perTF <- data.table(TFs = character(),
                          TAD = character(),
                          Adjusted.P.value = character(),
                          Motifs = character())
  
  topMotifs <- data.table(target = character(),
                          adjusted.p.value = numeric(),
                          id = character())
  
  iterations <- c(1: length(report.list))
  for (i in iterations){
    temp <- report.list[[i]]@d
    
    #exclude uncharacterized motifs
    temp <- temp[str_detect(temp$id,"UW.Motif.",negate = TRUE),]
    perTAD <- temp %>%
      dplyr::select(tad,target,id,adjusted.p.value, p.value) %>%
      dplyr::summarise(TAD = tad,TFs = paste(target, collapse = "|"),
                       Motifs = paste(id, collapse = "|"),
                       Adjusted.P.value = paste(adjusted.p.value, collapse = "|"),
                       P.value = paste(p.value, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    csv_perTAD <- rbind(csv_perTAD,perTAD)
    
    perTF.topMotif <- temp %>%
      dplyr::select(target,adjusted.p.value,id)
    
    topMotifs <- rbind(topMotifs, perTF.topMotif)
    
    perTF <- temp %>%
      dplyr::select(target,tad,adjusted.p.value,id) %>%
      group_by(target,tad) %>%
      dplyr::summarise(target, tad,
                       Adjusted.P.value = paste(adjusted.p.value, collapse = "|"),
                       Motifs = paste(id, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    colnames(perTF) <- c("TFs", "TAD", "Adjusted.P.value", "Motifs")
    csv_perTF <- rbind(csv_perTF,perTF)
  }
  
  data.visual <- csv_perTF %>%
    dplyr::select(TFs,TAD,Adjusted.P.value)
  
  csv_perTF <- csv_perTF %>%
    group_by(TFs) %>%
    dplyr::summarise(TFs,TAD = paste(TAD,collapse ="|"),
                     Adjusted.P.value = paste(Adjusted.P.value, collapse = "|"),
                     Motifs = paste(Motifs, collapse = "|"),) %>%
    as.data.table()%>%
    unique()
  
  iterations <- c(1:nrow(csv_perTF))
  for(i in iterations){
    
    temp <- csv_perTF[i] %>%
      dplyr::select(TFs,Motifs)
    temp <- temp %>%
      separate_rows(Motifs, sep = "\\|") %>%
      as.data.table()
    temp <- unique(temp)
    
    temp <- temp %>%
      group_by(TFs) %>%
      dplyr::summarise(TFs,
                       Motifs = paste(Motifs, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    csv_perTF[i,4] <- temp[1,2]
  }
  
  table.TFs.Motifs <- csv_perTF %>%
    dplyr::select(TFs,Motifs)
  data.visual <- left_join(data.visual,table.TFs.Motifs)
  
  topMotifs <- topMotifs %>%
    group_by(target, id) %>%
    summarise(target, id, adjusted.p.value = mean (adjusted.p.value)) %>%
    unique()
  
  topMotifs <- topMotifs %>%
    group_by(target) %>%
    summarise(target, Adjusted.P.value = min(adjusted.p.value), id, adjusted.p.value)
  
  topMotifs <- topMotifs[which(topMotifs$adjusted.p.value == topMotifs$Adjusted.P.value),]
  
  topMotifs <- dplyr::select(topMotifs, target, id)
  colnames(topMotifs) <- c("TFs","top.motif")
  
  data.visual <- left_join(data.visual,topMotifs)
  csv_perTF <- left_join(csv_perTF,topMotifs)
  
  newList <- list(table_perTAD = csv_perTAD,table_perTFs = csv_perTF, data.visual = data.visual)
  return(newList)
}


#This function is called by the "enrichmentAnalysis.R" script
#It performs enrichment analysis using the PWMEnrich tool
#PWMEnrich input is the DNA sequences grouped per TAD
motifEnrich <- function(biodata, motif_output_folder, dir_name,
                        p.adjust.method, cut.off, tech, exp.parent){
  
  #number of cores available for motif enrichment analysis
  #N <- 3
  
  #speed up execution
  #registerCoresPWMEnrich(N)
  #useBigMemoryPWMEnrich(TRUE)
  
  motif.data <- prepareSequences(biodata, tech, motif_output_folder, exp.parent, dir_name)
 
  seq.tad.number <- getDNASequences(motif.data,motif_output_folder, tech)
 
  #perform motif enrichment analysis using PWMEnrich
  report.list <- list()
  #report.list.adj <- list()
  
  # load the pre-compiled lognormal background
  data(PWMLogn.hg19.MotifDb.Hsap)
  #l <- 1
  k <- 1
  #iterations <- c(1:10)
  iterations <- c(1:nrow(seq.tad.number))
  for (i in iterations){
    
    sequence = readDNAStringSet(paste0(motif_output_folder,"/seq_perTADs.fasta"), format="fasta", skip = (seq.tad.number$start[i]-1), nrec = (seq.tad.number$end[i]-seq.tad.number$start[i]+1) )
    res = motifEnrichment(sequence, PWMLogn.hg19.MotifDb.Hsap)
    report = groupReport(res, by.top.motifs = TRUE)
    
    report@d$adjusted.p.value <- p.adjust(report@d$p.value, method = p.adjust.method)
    #report.p <- report[report$p.value < cut.off]
    report <- report[report$adjusted.p.value < cut.off]
    
    if (nrow(report@d)>0){
      report.list[[k]]<-report
      report.list[[k]]@d$tad <- seq.tad.number$tad[i]
      k <- k +1
    }
    
    #if (nrow(report.p@d)>0){
    #  report.list.p[[l]]<-report.p
    #  report.list.p[[l]]@d$tad <- seq.tad.number$tad[i]
     # l <- l+1
    #}
  }
 
  #registerCoresPWMEnrich(NULL)
  #useBigMemoryPWMEnrich(FALSE)
  #filename <- paste0(motif_output_folder,"/report MotifEA P value.txt")
  #file.create(filename, showWarnings = FALSE)
  #dput(report.list.p, file = filename)
  
  #filename <- paste0(motif_output_folder,"/report MotifEA adjusted P value.txt")
  #file.create(filename, showWarnings = FALSE)
  #dput(report.list.adj, file = filename)
  return(report.list)
  
}