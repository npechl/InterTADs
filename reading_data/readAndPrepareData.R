readAndPrepateData = function(pathsToData, colCmb = NULL, scale = 100, txdb = NULL){
  #------------------------------ Start time of algorithm ------------------------------
  # start = Sys.time()
  
  #------------------------------ Source required functions and libraries ------------------------------
  # print("Loading functions and libraries")
  # 
  # source("reading_data/check_CHROM_POS.R")
  # source("reading_data/check_CHROM_POS_ALT.R")
  # source("reading_data/readCSV.R")
  # source("reading_data/readTXT.R")
  # source("reading_data/readVCF.R")
  # source("reading_data/readVCF_snp_indel.R")
  # source("reading_data/readVCFindel.R")
  # source("reading_data/readVCFsnp.R")
  # source("preparing_data/createScale.R")
  # source("preparing_data/getAttConnection.R")
  # source("preparing_data/getBioLocation.R")
  # source("preparing_data/getGenomicFeatures.R")
  # 
  # library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  # library(GenomicFeatures)
  # library(systemPipeR)
  # library(data.table)
  # library(stringr)
  # library(stringi)
  # library(biomaRt)
  # library(gtools)
  # library(vcfR)
  # 
  # #------------------------------ Reading data ------------------------------
  # print("Reading data")
  # 
  # tablet = list()
  # dataVCF = list()
  # 
  # for(i in 1:length(pathsToData)){
  #   temp = pathsToData[[i]]
  #   
  #   if(str_detect(as.character(temp), ".csv")){
  #     tablet = c(tablet, readCSV(temp))
  #   } else if(str_detect(as.character(temp), ".txt")) {
  #     tablet = c(tablet, lapply(temp, readTXT))
  #   } else {
  #     dataVCF = c(dataVCF, lapply(temp, readVCF))
  #   }
  # }
  # 
  # biodata = tablet
  # 
  # rm(tablet)
  # 
  # # out = biodata
  # 
  # #------------------------------ Preparing data for integrating process ------------------------------
  # 
  # print("Column connection")
  # 
  # numMat = 1:length(biodata)
  # 
  # biodata = lapply(numMat, getAttConnection, biodata, colCmb, 0, scale)
  # 
  # numMat = 1:length(dataVCF) + length(biodata)
  # 
  # dataVCF = lapply(numMat, getAttConnection, dataVCF, colCmb, length(biodata))
  # 
  # print("Creating scale")
  # 
  # ret = createScale(biodata, dataVCF, scale)
  # biodata = ret[[1]]
  # dataVCF = ret[[2]]
  # 
  # if(is.null(colCmb)){
  #   biodata = rbindlist(biodata, use.names = FALSE)
  #   dataVCF = rbindlist(dataVCF, use.names = FALSE)
  #   colnames(biodata)[1] = "ID"
  # } else {
  #   biodata = rbindlist(biodata, use.names = TRUE)
  # }
  # 
  # rm(colCmb)
  # 
  # print("Getting location data")
  # 
  # bioLocation = getBioLocation(biodata$ID)
  # 
  # biodata = biodata[biodata$ID %in% bioLocation$ID, ]
  # 
  # biodata = biodata[order(as.character(biodata$ID)), ]
  # bioLocation = bioLocation[order(as.character(bioLocation$ID)), ]
  # 
  # biodata = cbind(bioLocation, biodata[ ,2:ncol(biodata)])
  # 
  # colnames(dataVCF) = colnames(biodata)
  # 
  # biodata = rbind(biodata, dataVCF)
  # 
  # rm(bioLocation)
  
  # print("Choosing desired chromosomes")
  #
  # if(is.null(chromosomes)){
  #   chr = as.factor(biodata$chromosome_name)
  #   chr = levels(chr)
  # } else {
  #   chr = as.factor(chromosomes)
  #   chr = levels(chr)
  # }
  #
  # biodata = biodata[which(biodata$chromosome_name %in% chr), ]
  
  biodata = pathsToData

  # biodata = biodata[sample(1:nrow(biodata), size = 10), ]
  
  all_chr = unique(biodata$chromosome_name)
  
  out = list()
  all_times = list()
  
  feat = genFeatures(txdb, reduce_ranges = FALSE)
  feat = unlist(feat)
  feat = feat[!(unlist(feat$featuretype) == "intergenic"), ]
  
  for(i in all_chr){
    time_elapsed = system.time({
      print(c(which(all_chr == i), i))
      
      temp_biodata = biodata[which(biodata$chromosome_name == i), ]
      
      gr = GRanges(seqnames = Rle(as.character(i)),
                   ranges = IRanges(start = as.numeric(temp_biodata$start_position),
                                    end = as.numeric(temp_biodata$end_position)))
      
      overlaps = findOverlaps(feat, gr)
      
      indexes = c(1:nrow(temp_biodata))
      
      genomic_features = lapply(indexes, getGenomicFeatures, overlaps, feat)
      genomic_features = rbindlist(genomic_features)
      
      temp = cbind(temp_biodata[ ,c(1:5)], genomic_features)
      temp_biodata = cbind(temp, temp_biodata[ ,6:ncol(biodata)])
      
      out = c(out, list(temp_biodata))
    })
    
    time_elapsed = as.data.frame(data.matrix(time_elapsed))
    time_elapsed = as.data.table(t(time_elapsed))
    
    time_elapsed = cbind(data.table(chromosome_name = i), time_elapsed)
    
    all_times = c(all_times, list(time_elapsed))
  }
  
  out = rbindlist(out)
  
  write.table(out, file = "all_pre_data_with_features.csv", row.names=FALSE, sep=";")
  
  all_times = rbindlist(all_times)
  
  write.table(all_times, file = "times_for_getting_features.csv", row.names = FALSE, sep = ";")
  
  return(all_times)
}