########## Loading libraries ########## 

rm(list = ls())

source("libraries.R")

########### Inputs ###########

#' Input parameters for TADiff part
#' 
#' @param dir_name Directory of input datasets containing feature counts and frequency tables
#' 
#' @param output_folder Folder name for printing output tables
#' 
#' @param meta meta-data file name used
#' 
#' @param names.meta meta data columns to process (names or indexes)
#' 
#' @param meth_data Parent index of methylation data. If no methylation is provided, place FALSE

dir_name = "Datasets_bloodcancer"

output_folder = "results_bloodcancer"

meta = "metaData_groups.csv"

meth_data = 2

data.all = fread(paste(output_folder, "/integrated-tad-table.csv", sep = ""),
                 sep = "\t")

meta = fread(paste(dir_name, meta, sep = "/"))
who = meta == ""
who = apply(who, 1, sum, na.rm = TRUE)
meta = meta[which(who == 0), ]


sample.list = meta$newNames


meth = data.all[which(data.all$parent == meth_data), ]


meth.promoter = meth[which(str_detect(meth$Gene_locus, "promoter")), ]

meth.intergenic = meth[which(str_detect(meth$Gene_locus, "intergenic")), ]

meth.regulatory = rbind(meth.promoter, meth.intergenic)
meth.regulatory = meth.regulatory[!duplicated(meth.regulatory$ID), ]

rm(meth.intergenic, meth.promoter, meth)

for(i in sample.list){
  meth.regulatory[[i]] = 100 - meth.regulatory[[i]]
}

other = data.all[which(data.all$parent != meth_data), ]

data.all = rbind(meth.regulatory, other)

write.table(data.all, 
            file = paste(output_folder, "/integrated-tad-table-methNorm.txt", sep = ""), 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

#' # overlap of the TAD with events
#' TAD = fread(paste("D:\\BACKUP\\Data_Integration", "/hglft_genome_2dab_ec1330.bed", sep = ""), 
#'             header = F, sep = "\t", check.names = FALSE)
#' 
#' data.all$name = paste0(data.all$ID, data.all$Gene_id)
#' 
#' #' Reordering
#' #' 
#' data.over = data.all[,c("name", "chromosome_name", "start_position", "end_position")]
#' data.over$chromosome_name = paste("chr", data.over$chromosome_name, sep = "")
#' data.over = data.over[,c(2,3,4,1)]
#' 
#' #'
#' #' Overlap of the TAD with events
#' #' 
#' colnames(TAD) = paste(c("chr", "start", "end", "name"))
#' colnames(data.over) = paste(c("chr", "start", "end", "name"))
#' 
#' #'
#' #' Make GRanges object
#' #' 
#' gr1 = with(TAD, GRanges(chr, IRanges(start = start, end = end, names = name)))
#' gr2 = with(data.over, GRanges(chr, IRanges(start = start, end = end, names = name)))
#' 
#' #'
#' #' Completely overlapping
#' #' 
#' type1 = findOverlaps(query = gr1, subject = gr2)
#' 
#' type1.df = cbind(TAD[queryHits(type1),], data.over[subjectHits(type1),])
#' type1.df = type1.df[,c(2,3,4,8)]
#' colnames(type1.df) = c("tad_start", "tad_end", "tad_name", "name.1")
#' 
#' #'
#' #' Generating output overlapping with
#' #' TADs table
#' #' 
#' full = base::merge(type1.df, data.all, by.x = "name.1", by.y = "name")
#' 
#' 
#' #write.table(full,file=paste("integrated_table_with_TADs_regulatory.csv"), col.names=T, row.names=F, quote=F, sep="\t")
#' #####################################################
#' ##mean table
#' colnames(full)
#' full<- full[,-1]
#' 
#' expr<- full %>% 
#'   filter(str_detect(ID, 'ENSG')) 
#' colnames(expr)
#' 
#' t.1 = aggregate(expr[,10:144], list(expr$tad_name), mean)
#' head(t.1)
#' 
#' write.table(t.1,file=paste("exprmean_in_tad_regulatory.txt"), col.names=T, row.names=F, quote=F, sep="\t")
#' #events_mean_in_tad.txt", sep = ""
