bioIntegrate = function(data, chromosomes){
  source("integrating_process/chrProcessing.R")
  source("integrating_process/rangeProcessing.R")
  
  library(data.table)
  
  out = lapply(chromosomes, chrProcessing, biodata)
  out = rbindlist(out)
  out = as.data.table(out)
  
  return(out)
}