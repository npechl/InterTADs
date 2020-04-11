########## Loading libraries ##########
library(systemPipeR)
library(data.table)
library(stringr)
library(dplyr)

library(tidyr)

library(org.Hs.eg.db)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("helpers.R")

start_time = Sys.time()

############ Inputs ############
dir_name = "Datasets"
output_folder = "output_tables"

tech = "hg19"

meta = fread(paste(dir_name, "/meta-data.csv", sep = ""))
who = meta == ""
who = apply(who, 1, sum)
meta = meta[which(who == 0), ]

counts = list.files(paste(dir_name, "/counts", sep = ""))
freq = list.files(paste(dir_name, "/freq", sep = ""))

files = colnames(meta)
files = files[which(!(files %in% c("groups", "newNames")))]

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
  
  if(max(new[,5:ncol(new)]) <= 1){
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

########### Filtering ########### 

biodata = biodata[which(biodata$chromosome_name %in% as.character(1:22)), ]

zero.ids = rep(0, length(names))
zero.ids = paste(zero.ids, collapse = "")

data.num.ids = unite(data = biodata[,..names], col = ids, sep = "")
biodata = biodata[which(data.num.ids != zero.ids), ]

############ Getting genomic features ############ 

if(tech == "hg19"){
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
} else {
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
}

feat = genFeatures(txdb, reduce_ranges = FALSE, verbose = FALSE)
feat = feat@unlistData

gr = GRanges(seqnames = Rle(paste("chr", biodata$chromosome_name, sep = "")),
             ranges = IRanges(start = as.numeric(biodata$start_position),
                              end = as.numeric(biodata$end_position)))

overlaps = findOverlaps(gr, feat)

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

mapping = map_entrez_ids(entrez.ids = genes, tech = tech)

map = base::match(biodata$feature_by, mapping$entrez.id)
biodata$feature_by = mapping[map, ]$hgnc.symbol

unique.ids = unique(biodata$ID)
iterations = seq(from = 0, to = (length(unique.ids) - 1), by = 200000)
iterations = c(iterations, length(unique.ids))

total = list()

for(i in 1:(length(iterations) - 1)){
  
  temp = unique.ids[(iterations[i] + 1):(iterations[i + 1])]
  
  temp = biodata[which(biodata$ID %in% temp), ]
  
  features = temp %>% 
             dplyr::group_by(ID) %>% 
             dplyr::select(ID, feature_by, featuretype) %>% 
             dplyr::summarise(Gene_id = paste(feature_by, collapse = "|"),
                              Gene_feature = paste(featuretype, collapse = "|"))
  
  features = as.data.table(features)
  
  who = which(!duplicated(temp$ID))
  temp = temp[who, ]
  
  map = base::match(features$ID, temp$ID)
  temp[map, ]$feature_by = features$Gene_id
  temp[map, ]$featuretype = features$Gene_feature
  
  total[[i]] = temp
}

biodata = rbindlist(total)

rm(total)

names = colnames(biodata)
names = str_replace(names, "feature_by", "Gene_id")
names = str_replace(names, "featuretype", "Gene_locus")

colnames(biodata) = names

rm(list=setdiff(ls(), c("biodata", "meta", "start_time", "dir_name", "output_folder")))

end_time = Sys.time()

############ Generating outputs ############ 

dir.create(output_folder, showWarnings = FALSE)

write.table(biodata, paste(output_folder, "/integrated_table.csv", sep = ""), 
            row.names = FALSE, sep = "\t", quote = FALSE)
