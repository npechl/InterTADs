# Loading libraries ------------------------------------------------------------

rm(list = ls())

source("R/libraries.R")
source("R/helpers.R")

start_time = Sys.time()

# Inputs -----------------------------------------------------------------------

#'
#' Input parameters for Data Integration part
#'
#' @param dir_name Directory of input datasets containing feature counts and frequency tables
#'
#' @param output_folder Folder name for printing output tables
#'
#' @param tech Human Genome Reference used
#'
#' @param meta meta-data file name used
#'
#' @param counts_dir Directory name of counts NGS data
#'
#' @param freq_dir Directory name of freq NGS data
#'
#' @param tad_file BED file containing information about TADs
#'

dir_name = "Datasets_bloodcancer/"

output_folder = "results_bloodcancer/"

tech = "hg19" # or "hg38"

meta = "metaData_groups.csv"

counts_dir = "counts"

freq_dir = "freq"

tad_file = "hglft_genome_2dab_ec1330.bed"

# Reading files ----------------------------------------------------------------

#'
#' Reading meta data file
#'

meta = fread(paste(dir_name, meta, sep = "/"))
who = meta == ""
who = apply(who, 1, sum, na.rm = TRUE)
meta = meta[which(who == 0), ]

#'
#' Getting input files
#'

counts = list.files(paste(dir_name, counts_dir, sep = "/"))
freq = list.files(paste(dir_name, freq_dir, sep = "/"))

files = colnames(meta)
# files = files[which(!(files %in% c("groups", "newNames")))]

names = meta$newNames

biodata = list()

#'
#' Reading input frequency tables
#'

if(length(freq) > 0){
  for(i in 1:length(freq)){
    new = fread(paste(dir_name, freq_dir, freq[i], sep = "/"))

    keep = which(str_detect(freq[i], files))
    parent = keep
    keep = meta[,..keep]
    keep = as.character(unlist(keep))

    new = cbind(new[,1:4], new[,..keep])
    colnames(new) = c("ID", "chromosome_name", "start_position", "end_position", names)

    new$chromosome_name = str_to_lower(new$chromosome_name)
    new$chromosome_name = str_remove(new$chromosome_name, "chr")

    if(max(new[,5:ncol(new)]) <= 1){
      new[,5:ncol(new)] = 100 * new[,5:ncol(new)]
    }

    new$parent = parent

    biodata[[i]] = new
  }
}

#'
#' Reading input counts tables
#'

if(length(counts) > 0){
  for(i in 1:length(counts)){
    new = fread(paste(dir_name, counts_dir, counts[i], sep = "/"))

    keep = which(str_detect(counts[i], files))
    parent = keep
    keep = meta[,..keep]
    keep = as.character(unlist(keep))

    new = cbind(new[,1:4], new[,..keep])
    colnames(new) = c("ID", "chromosome_name", "start_position", "end_position", names)

    new$chromosome_name = str_to_lower(new$chromosome_name)
    new$chromosome_name = str_remove(new$chromosome_name, "chr")

    temp = new[,5:ncol(new)]

    temp = log(temp + 1)
    col.max = apply(temp, 2, max)
    col.max = as.numeric(col.max)

    temp = t(temp) * (100 / col.max)
    temp = as.data.table(t(temp))

    new = cbind(new[,1:4], temp)

    new$parent = parent

    biodata[[i + length(freq)]] = new
  }
}

biodata = rbindlist(biodata, use.names = FALSE)

# Filtering --------------------------------------------------------------------

#'
#'  Keep chromosomes 1 - 22
#'  
biodata = biodata[which(biodata$chromosome_name %in% as.character(1:22)), ]

#'
#' Remove all-zeros records
#' 
# zero.ids = rep(0, length(names))
# zero.ids = paste(zero.ids, collapse = "")

data.num.ids = base::rowSums(biodata[, meta$newNames, with = FALSE]) # unite(data = biodata[,..names], col = ids, sep = "")
biodata = biodata[which(data.num.ids != 0), ]

# Getting genomic features -----------------------------------------------------

#'
#' Getting gene names (expressed as entrez ids) and locus (e.g. exon, intron, cds etc.)
#' based on chromosomal location
#' 
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

overlaps.from = overlaps@from
overlaps.to = overlaps@to

rm(overlaps)

feat = as.data.table(feat)
feat = feat[,c("feature_by", "featuretype")]
feat$feature_by = as.character(feat$feature_by)


rm(gr, temp, new, data.num.ids, txdb)

feat = cbind(biodata[overlaps.from, 1:4], feat[overlaps.to, ])
biodata = cbind(feat, biodata[overlaps.from, 5:ncol(biodata)])

rm(feat, overlaps.from, overlaps.to, col.max, keep, names)

biodata = unique(biodata)
biodata[which(biodata$feature_by == "character(0)"), ]$feature_by = "NA"
biodata[which(biodata$feature_by == ""), ]$feature_by = "NA"

#'
#' Converting entrez ids to hgnc symbols
#' 

genes = unique(biodata$feature_by)
genes = genes[which(genes != "NA")]
genes = genes[!str_detect(genes, "INTER")]

mapping = map_entrez_ids(entrez.ids = genes, tech = tech)

map = base::match(biodata$feature_by, mapping$entrez.id)
biodata$feature_by = mapping[map, ]$hgnc.symbol

rm(map, mapping, genes)

#'
#' Collapsing genes and feature type of same observation into same row
#' 
#' e.g.
#' 
#' index1 --- MYC, exon
#' index1 --- RUNX1T1, intron
#' --------------------------
#' index1 --- MYC|RUNX1T1, exon|intron
#' 

features = biodata[, .(Gene_id = paste(feature_by, collapse = "|"),
                       Gene_feature = paste(featuretype, collapse = "|")), by = ID]

biodata = biodata[which(!base::duplicated(biodata$ID)), ]

who = base::match(biodata$ID, features$ID)

biodata$feature_by = features[who, ]$Gene_id
biodata$featuretype = features[who, ]$Gene_feature


# unique.ids = unique(biodata$ID)
# iterations = seq(from = 0, to = (length(unique.ids) - 1), by = 200000)
# iterations = c(iterations, length(unique.ids))
# 
# total = list()
# 
# for(i in 1:(length(iterations) - 1)) {
#   
#   temp = unique.ids[(iterations[i] + 1):(iterations[i + 1])]
#   
#   temp = biodata[which(biodata$ID %in% temp), ]
#   
#   features = temp %>% 
#              dplyr::group_by(ID) %>% 
#              dplyr::select(ID, feature_by, featuretype) %>% 
#              dplyr::summarise(Gene_id = paste(feature_by, collapse = "|"),
#                               Gene_feature = paste(featuretype, collapse = "|"))
#   
#   features = as.data.table(features)
#   
#   who = which(!duplicated(temp$ID))
#   temp = temp[who, ]
#   
#   map = base::match(features$ID, temp$ID)
#   temp[map, ]$feature_by = features$Gene_id
#   temp[map, ]$featuretype = features$Gene_feature
#   
#   total[[i]] = temp
#   
# }
# 
# biodata = rbindlist(total)

rm(features)

names = colnames(biodata)
names = str_replace(names, "feature_by", "Gene_id")
names = str_replace(names, "featuretype", "Gene_locus")

colnames(biodata) = names

rm(list = setdiff(ls(), c("biodata", "meta", "start_time", "dir_name", "output_folder", "x", "res", "tad_file")))

# Annotate with TADs -----------------------------------------------------------

TAD = fread(paste(dir_name, tad_file, sep = "/"), 
            header = F, 
            sep = "\t", 
            check.names = FALSE)

#'
#' Create a barcode
#'  
# biodata$name = paste0(biodata$ID, biodata$Gene_id)

#' Reordering
#' 
data.over = biodata[,c("ID", "chromosome_name", "start_position", "end_position")]
data.over$chromosome_name = paste("chr", data.over$chromosome_name, sep = "")
data.over = data.over[,c(2,3,4,1)]

#'
#' Overlap of the TAD with events
#' 
colnames(TAD) = paste(c("chr", "start", "end", "name"))
colnames(data.over) = paste(c("chr", "start", "end", "name"))

#'
#' Make GRanges object
#' 
gr1 = with(TAD, GRanges(chr, IRanges(start = start, end = end, names = name)))
gr2 = with(data.over, GRanges(chr, IRanges(start = start, end = end, names = name)))

#'
#' Completely overlapping
#' 
type1 = findOverlaps(query = gr1, subject = gr2)

type1.df = cbind(TAD[queryHits(type1),], data.over[subjectHits(type1),])
type1.df = type1.df[,c(2,3,4,8)]
colnames(type1.df) = c("tad_start", "tad_end", "tad_name", "name.1")

rm(type1, gr1, gr2)

#' Overlapping with TADs' table
#' 
full = base::merge(type1.df, biodata, by.x = "name.1", by.y = "ID")

colnames(full)[1] = "ID"

keep = c("chromosome_name", 
         "tad_name", 
         "tad_start", 
         "tad_end", 
         "ID", 
         "start_position", 
         "end_position", 
         "Gene_id", 
         "Gene_locus", 
         "parent",
         meta$newNames)

full = full[,..keep]

rm(list = setdiff(ls(), c("biodata", "full", "meta", "start_time", "dir_name", "output_folder", "x", "res")))

end_time = Sys.time()

# Generating outputs -----------------------------------------------------------

dir.create(output_folder, showWarnings = FALSE)

write.table(biodata, paste(output_folder, "/integrated-table.csv", sep = ""),
            row.names = FALSE, sep = "\t", quote = FALSE)

write.table(full, paste(output_folder, "/integrated-tad-table.csv", sep = ""),
            row.names = FALSE, sep = "\t", quote = FALSE)

for(i in unique(biodata$parent)){

  temp = biodata[which(biodata$parent == i), ]$ID

  temp = temp[sample(1:length(temp), 3)]

  temp = paste(temp, collapse = ", ")

  cat(c("File", i, ":", temp, "...", "\n\n"),
      file = paste(output_folder, "/summary.txt", sep = ""),
      append = TRUE)



}
