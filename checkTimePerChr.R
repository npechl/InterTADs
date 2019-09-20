library(data.table)
library(gtools)

source("integrating_process/chrProcessing.R")
source("integrating_process/rangeProcessing.R")

pre_data = fread("all_pre_data/all_pre_data_with_features.csv", sep = ";")

# pre_data = pre_data[sample(1:nrow(pre_data), size = 10), ]

all_chr = unique(pre_data$chromosome_name)

rm(pre_data)

all_chr = mixedsort(all_chr)

all_times = list()

for(i in all_chr){
  print(c(i, which(all_chr == i)))
  
  current_integration = read.csv(file = paste("all_pre_data/pre_integration_chr_", i, ".csv", sep = ""), sep = ";")
  current_integration = as.data.table(current_integration)
  
  time_for_chr = system.time({
    current_integration = chrProcessing(i, current_integration)
  })
  
  rm(current_integration)

  time_for_chr = as.data.frame(data.matrix(time_for_chr))
  time_for_chr = as.data.table(t(time_for_chr))

  time_for_chr = cbind(data.table(chromosome_name = i), time_for_chr)

  all_times = c(all_times, list(time_for_chr))
}

all_times = rbindlist(all_times)
write.table(all_times, file = "time_pre_chr_integration.csv", sep = ";")