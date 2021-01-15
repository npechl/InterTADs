########## Loading libraries ########## 

rm(list = ls())

source("libraries.R")

start_tad_time = Sys.time()

########### Inputs ###########

#' Input parameters for TADiff part
#' 
#' @param dir_name Directory of input datasets containing feature counts and frequency tables
#' 
#' @param output_folder Folder name for printing output tables
#' 
#' @param image_output_folder Folder name for print output images
#' 
#' @param meta meta-data file name used
#' 
#' @param names.meta meta data columns to process (names or indexes)
#' 
#' @param expr_data Parent index of expression data. If no expression is provided, place FALSE

dir_name = "Datasets"

output_folder = "results"

image_output_folder = "results-vis"

meta = "meta-data.csv"

names.meta = c("groups")

expr_data = 2 



data.all = fread(paste(output_folder, "/integrated-tad-table.csv", sep = ""),
                 sep = "\t")

data.all$ID = paste(data.all$tad_name, data.all$ID, sep = ";")

meta = fread(paste(dir_name, meta, sep = "/"))
who = meta == ""
who = apply(who, 1, sum, na.rm = TRUE)
meta = meta[which(who == 0), ]


sample.list = meta$newNames

who = c("ID", "Gene_id", "parent")

gene.parents = data.all[,..who]

colnames(gene.parents) = c("ID", "gene_ID", "parents")


diffevent = matrix(data = "0", 
                   nrow = length(names.meta), 
                   ncol = (length(unique(data.all$parent)) + 1))

colnames(diffevent) = c("Analysis", paste("p_", unique(data.all$parent), sep = ""))

rm(who)

for (z in 1:length(names.meta)){
  
  analysis = names.meta[z]
  groups = as.character(meta[[analysis]])
  groups = unique(groups)
  groups = groups[which(groups != "")]
  groups = groups[!is.na(groups)]
  
  group1 = groups[1]
  group2 = groups[2]
  
  meta.keep = meta[which(meta[[analysis]] == group1 | meta[[analysis]] == group2), ]
  
  sample.list = meta.keep$newNames
  
  df = data.all[,..sample.list]
  
  df = as.data.frame(df)
  row.names(df) = data.all$ID
  
  pheno = as.factor(meta.keep[[analysis]])
  
  phenoMat = model.matrix(~pheno)
  colnames(phenoMat) = sub("^pheno", "", colnames(phenoMat))
  
  fit = lmFit(object = df, design = phenoMat)
  
  gc()
  set.seed(6)
  fit = eBayes(fit)
  
  gc()
  top.rank = topTable(fit, number = nrow(df), adjust.method = "fdr", sort.by = "p")
  
  sign.table = top.rank[which(top.rank$adj.P.Val <= 0.01 & abs(top.rank$logFC) > 4), ]
  
  if (nrow(sign.tad.info) == 0) {
    
    cat(c("No statistical significant events for:", names.meta[z], "\n"))
    
    diffevent[z,1] = analysis 
    
  } else {
    
    # annotate sign.table
    
    sign.table$ID = row.names(sign.table)
    
    sign.genes = merge(gene.parents, 
                       sign.table, 
                       by.x = "ID",
                       by.y = "ID")
    
    sign.table = merge(sign.table, 
                       data.all, 
                       by.x = "ID", 
                       by.y = "ID")
    
    sign = sign.genes %>% count(parents)
    
    diffevent[z, 1] = analysis
    diffevent[z, paste("p_", sign$parents, sep = "")] = sign$n
    
    # final file - export
    write.table(sign.table,
                file = paste(output_folder, "/", analysis, "_evenDiff.txt", sep = ""), 
                col.names = TRUE, 
                row.names = FALSE, 
                quote = FALSE, 
                sep = "\t")
  } 
}

diffevent = as.data.frame(diffevent)

for(i in 2:ncol(diffevent)){
  diffevent[,i] = as.numeric(diffevent[,i])
}

diffevent$`total events` = rowSums(diffevent[,2:ncol(diffevent)])

write.table(diffevent, 
            file = paste(output_folder, "/Summary_evenDiff.txt", sep = ""), 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")






