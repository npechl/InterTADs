library(data.table)
library(gtools)

source("integrating_process/chrProcessing.R")
source("integrating_process/rangeProcessing.R")

pre_data = fread("all_pre_data/all_pre_data_with_features.csv", sep = ";")

# pre_data = pre_data[sample(1:nrow(pre_data), size = 10), ]

all_chr = unique(pre_data$chromosome_name)

all_chr = mixedsort(all_chr)

rowsPerchr = data.table(numOfRows = numeric(), chromosome_name = character())

for(i in all_chr){
  print(c(which(all_chr == i), i))
  temp = nrow(pre_data[which(pre_data$chromosome_name == i), ])
  
  rowsPerchr = rbind(rowsPerchr, data.table(numOfRows = temp, chromosome_name = i))
}

write.table(rowsPerchr, file = "rowsPerChr.csv", sep = ";")