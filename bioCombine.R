bioCombine = function(biodata, colCmb = NULL, scale = NULL, chromosomes = NULL, txdb){
  #------------------------------ Start time of algorithm ------------------------------
  start = Sys.time()
  
  #------------------------------ Source required functions and libraries ------------------------------
  print("Loading functions and libraries")

  source("reading_data/check_CHROM_POS.R")
  source("reading_data/check_CHROM_POS_ALT.R")
  source("reading_data/readCSV.R")
  source("reading_data/readTXT.R")
  source("reading_data/readVCF.R")
  source("reading_data/readVCF_snp_indel.R")
  source("reading_data/readVCFindel.R")
  source("reading_data/readVCFsnp.R")
  source("preparing_data/createScale.R")
  source("preparing_data/getAttConnection.R")
  source("preparing_data/getBioLocation.R")
  source("preparing_data/getGenomicFeatures.R")
  source("integrating_process/chrProcessing.R")
  source("integrating_process/rangeProcessing.R")
  
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(GenomicFeatures)
  library(systemPipeR)
  library(data.table)
  library(stringr)
  library(stringi)
  library(biomaRt)
  library(gtools)
  library(vcfR)
  
  #------------------------------ Reading data ------------------------------
  print("Reading data")

  tablet = list()
  dataVCF = list()
  
  for(i in 1:length(biodata)){
    temp = biodata[[i]]
    
    if(str_detect(paste(temp, collapse = ""), ".csv")){
      tablet = c(tablet, readCSV(temp))
    } else if(str_detect(paste(temp, collapse = ""), ".txt")) {
      tablet = c(tablet, readTXT(temp))
    } else {
      dataVCF = c(dataVCF, readVCF(temp))
    }
  }
  
  biodata = tablet
  
  rm(tablet)
  
  # out = biodata
  
  #------------------------------ Preparing data for integrating process ------------------------------

  print("Column connection")
  
  numMat = 1:length(biodata)

  biodata = lapply(numMat, getAttConnection, biodata, colCmb, 0)
  
  numMat = 1:length(dataVCF) + length(biodata)
  
  dataVCF = lapply(numMat, getAttConnection, dataVCF, colCmb, length(biodata))
  
  print("Creating scale")

  if(!is.null(scale)){
    ret = createScale(biodata, dataVCF, scale)
    biodata = ret[[1]]
    dataVCF = ret[[2]]
    rm(ret)
  }

  if(is.null(colCmb)){
    biodata = rbindlist(biodata, use.names = FALSE)
    dataVCF = rbindlist(dataVCF, use.names = FALSE)
    colnames(biodata)[1] = "ID"
  } else {
    biodata = rbindlist(biodata, use.names = TRUE)
  }

  rm(colCmb)

  print("Getting location data")

  bioLocation = getBioLocation(biodata$ID)

  biodata = biodata[biodata$ID %in% bioLocation$ID, ]

  biodata = biodata[order(as.character(biodata$ID)), ]
  bioLocation = bioLocation[order(as.character(bioLocation$ID)), ]

  biodata = cbind(bioLocation, biodata[ ,2:ncol(biodata)])
  
  colnames(dataVCF) = colnames(biodata)
  
  biodata = rbind(biodata, dataVCF)
  
  rm(bioLocation)

  print("Choosing desired chromosomes")

  if(is.null(chromosomes)){
    chr = as.factor(biodata$chromosome_name)
    chr = levels(chr)
  } else {
    chr = as.factor(chromosomes)
    chr = levels(chr)
  }

  biodata = biodata[which(biodata$chromosome_name %in% chr), ]
  
  print("Getting genomic features")
  
  gr = GRanges(seqnames = Rle(as.character(biodata$chromosome_name)), 
               ranges = IRanges(start = as.numeric(biodata$start_position), 
                                end = as.numeric(biodata$end_position)))
  feat = genFeatures(txdb, reduce_ranges = FALSE)
  feat = unlist(feat)
  feat = feat[!(unlist(feat$featuretype) == "intergenic"), ]
  
  overlaps = findOverlaps(feat, gr)
  
  indexes = c(1:nrow(biodata))
  
  genomic_features = lapply(indexes, getGenomicFeatures, overlaps, feat)
  genomic_features = rbindlist(genomic_features)
  
  temp = cbind(biodata[ ,c(1:5)], genomic_features)
  biodata = cbind(temp, biodata[ ,6:ncol(biodata)])
  
  rm(gr, feat, overlaps, indexes, genomic_features, temp)
  
  
  write.table(biodata, file = paste("biodata_pre_integration_", paste(chromosomes, collapse = "."), ".csv", sep = ""), row.names=FALSE, sep=";")

  #------------------------------ Start integrating process ------------------------------

  print("Start of Integration process")

  out = lapply(chr, chrProcessing, biodata)
  out = rbindlist(out)
  out = as.data.table(out)
  write.table(out, file = paste("biodata_post_integration_", chromosomes, ".csv", sep = ""), row.names=FALSE, sep=";")

  print("End of Integration process")

  #------------------------------ End of algorithm ------------------------------

  end = Sys.time()
  print(end - start)

  return(out)
}