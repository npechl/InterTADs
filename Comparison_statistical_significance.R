library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(gplots)
library(data.table)

data.all = fread("output_tables/integrated_table.csv")
meta = fread("Data_Integration/meta-data.csv")
tad_sum = fread("output_tables/tad_statistics.csv")

expr = data.all[which(str_detect(ID, 'XLOC')), ]

meth = data.all[which(str_detect(ID, 'cg')), ]

var = data.all[which(!str_detect(ID, 'cg') & !str_detect(ID, 'XLOC'))]

groups = unique(meta$groups)

list1 = meta[which(meta$groups == groups[1]), ]$newNames
list2 = meta[which(meta$groups == groups[2]), ]$newNames

new.names = meta$newNames

one.run <- function(index){
  
  statistics <- wilcox.test(unlist(tad.1[index, ]), 
                            unlist(tad.2[index, ]), 
                            na.rm=TRUE, 
                            correct = FALSE)
  
  
  return(statistics$p.value) 
}

########## expr ########## 

tad.1 = expr[,..list1]
tad.2 = expr[,..list2]

indexes = 1:nrow(expr)
p.values = lapply(indexes, one.run)
p.values = unlist(p.values)

expr$pvalue = p.values

data.b = as.data.frame(as.numeric(expr$pvalue))
expr$FDR = p.adjust(unlist(data.b),
                    method = c("BH"),
                    n = nrow(data.b))

########## meth ########## 

tad.1 = meth[,..list1]
tad.2 = meth[,..list2]

indexes = 1:nrow(meth)
p.values = lapply(indexes, one.run)
p.values = unlist(p.values)

meth$pvalue = p.values

data.b = as.data.frame(as.numeric(meth$pvalue))
meth$FDR = p.adjust(unlist(data.b),
                    method = c("BH"),
                    n = nrow(data.b))

########## var ########## 

tad.1 = var[,..list1]
tad.2 = var[,..list2]
tad.all = cbind(tad.1, tad.2)

ids = unite(tad.all, col = "id", sep = "")
tad.all = cbind(tad.all, ids)

uniq.tad.all = unique(tad.all)
tad.1 = uniq.tad.all[,..list1]
tad.2 = uniq.tad.all[,..list2]

indexes = 1:nrow(uniq.tad.all)

p.values = lapply(indexes, one.run)
p.values = unlist(p.values)
uniq.tad.all$p.values = p.values

map = base::match(tad.all$id, uniq.tad.all$id)

tad.all$p.values = uniq.tad.all[map, ]$p.values
var$pvalues = tad.all$p.values

data.b = as.data.frame(var$pvalue)
var$FDR = p.adjust(unlist(data.b),
                   method = c("BH"), 
                   n = nrow(data.b))



meth$category<- paste("DNA methylation")
expr$category<- paste("gene expression")
var$category<- paste("variants")
tad_sum$category<- paste("TADs")

########## Generating output image ########## 

table.FDR = rbind(expr[,c("FDR", "category")],
                  meth[,c("FDR", "category")],
                  var[,c("FDR", "category")],
                  tad_sum[,c("FDR", "category")])

table.FDR = table.FDR[which(!is.na(table.FDR$FDR)), ]

gr  = ggplot(table.FDR, aes(y = category, x = FDR, color = category)) +
      scale_color_brewer(palette="Dark2") +
      geom_boxplot() +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


