########## Loading libraries ##########
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(systemPipeR)
library(data.table)
library(stringr)
# library(biomaRt)
library(dplyr)

library(org.Hs.eg.db)
# library(EnsDb.Hsapiens.v79)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("helpers.R")

start_time = Sys.time()

############ Inputs ############
dir_name = "Data Integration"

tech = "h19"

meta = fread(paste(dir_name, "/meta-data.csv", sep = ""))
files = colnames(meta)
files = files[which(!(files %in% c("groups", "newNames")))]

counts = list.files(paste(dir_name, "/counts", sep = ""))
freq = list.files(paste(dir_name, "/freq", sep = ""))

names = meta$newNames

############ Reading files ############ 

biodata = list()

for(i in 1:length(freq)){
  new = fread(paste(dir_name, "/freq/", freq[i], sep = ""))
  
  keep = which(str_detect(freq[i], files))
  keep = meta[,..keep]
  keep = as.character(unlist(keep))
  
  new = cbind(new[,1:4], new[,..keep])
  colnames(new) = c("ID", "chromosome_name", "start_position", "end_position", names)
  
  new$chromosome_name = str_remove(new$chromosome_name, "chr")
  new$chromosome_name = str_remove(new$chromosome_name, "Chr")
  
  if(max(new[,5:ncol(new)]) == 1){
    new[,5:ncol(new)] = 100 * new[,5:ncol(new)]
  }
  
  biodata[[i]] = new
}


for(i in 1:length(counts)){
  new = fread(paste(dir_name, "/counts/", counts[i], sep = ""))
  
  keep = which(str_detect(counts[i], files))
  keep = meta[,..keep]
  keep = as.character(unlist(keep))
  
  new = cbind(new[,1:4], new[,..keep])
  colnames(new) = c("ID", "chromosome_name", "start_position", "end_position", names)
  
  new$chromosome_name = str_remove(new$chromosome_name, "chr")
  new$chromosome_name = str_remove(new$chromosome_name, "Chr")
  
  temp = new[,5:ncol(new)]
  
  temp = log(temp + 1)
  col.max = apply(temp, 2, max)
  col.max = as.numeric(col.max)
  
  temp = t(temp) * (100 / col.max)
  temp = as.data.table(t(temp))
  
  new = cbind(new[,1:4], temp)
  
  biodata[[i + length(freq)]] = new
}

biodata = rbindlist(biodata, use.names = FALSE)

one.run <- function(biodata){
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
  feat = genFeatures(txdb, reduce_ranges = FALSE, verbose = FALSE)
  
  gr = GRanges(seqnames = Rle(paste("chr", biodata$chromosome_name, sep = "")),
               ranges = IRanges(start = as.numeric(biodata$start_position),
                                end = as.numeric(biodata$end_position)))
  
  overlaps = findOverlaps(gr, feat@unlistData)
  
  feat = feat@unlistData
  feat = as.data.table(feat)
  feat = feat[,c("feature_by", "featuretype")]
  feat$feature_by = as.character(feat$feature_by)
  
  feat = cbind(biodata[overlaps@from, 1:4], feat[overlaps@to, ])
  biodata = cbind(feat, biodata[overlaps@from, 5:ncol(biodata)])
  
  biodata = unique(biodata)
  biodata[which(biodata$feature_by == "character(0)"), ]$feature_by = "NA"
  biodata[which(biodata$feature_by == ""), ]$feature_by = "NA"
  
  genes = unique(biodata$feature_by)
  genes = genes[which(genes != "NA")]
  genes = genes[!str_detect(genes, "INTER")]
  
  mapping = map_entrez_ids(entrez.ids = genes)
  
  map = BiocGenerics::match(biodata$feature_by, mapping$entrez.id)
  biodata$feature_by = mapping[map, ]$hgnc.symbol
  
  features = biodata %>% group_by(ID) %>% dplyr::select(ID, feature_by, featuretype) %>% summarise(Gene_id = paste(feature_by, collapse = "|"),
                                                                                                   Gene_feature = paste(featuretype, collapse = "|"))
  
  features = as.data.table(features)
  
  who = which(!duplicated(biodata$ID))
  biodata = biodata[who, ]
  
  map = match(features$ID, biodata$ID)
  biodata[map, ]$feature_by = features$Gene_id
  biodata[map, ]$featuretype = features$Gene_feature
  
  names = colnames(biodata)
  names = str_replace(names, "feature_by", "Gene_id")
  names = str_replace(names, "featuretype", "Gene_locus")
  
  colnames(biodata) = names
  
  genes = str_split(biodata$Gene_id, "\\|")
  genes = unlist(genes)
  
  return()
}

results = list()

keep = c(10000, 100000, 150000, 200000, 300000, 500000, seq(1000000, 4500000, by=500000))

for(i in keep){
  print(i)
  mydata = biodata[sample(1:nrow(biodata), size = i, replace = TRUE), ]
  count.time = system.time(one.run(biodata = mydata))
  count.time = as.data.table(t(data.matrix(count.time)))[,1:3]
  count.time$num_rows = nrow(mydata)
  
  results[[i]] = count.time
}

results = rbindlist(results)

plot(x = results$num_rows, 
     y = results$user.self, 
     type = "l", ylab = "sec", xlab = "number of rows", 
     col = "gold3")

lines(x = results$num_rows, 
      y = results$sys.self,  
      type = "l", col = "blue")

lines(x = results$num_rows, 
      y = results$elapsed,  
      type = "l", col = "red")

legend("topleft", inset=0.05, legend=c("Elapsed", "System time", "User time"),
       col=c("red", "blue", "gold3"), lwd = 1, cex=0.8)
