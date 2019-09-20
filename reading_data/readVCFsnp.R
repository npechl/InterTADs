readVCFsnp = function(name){
  library(vcfR)
  library(data.table)
  library(stringi)
  library(stringr)
  
  f = read.vcfR(as.character(name), verbose = FALSE)
  out = f@fix[,c(1:3, 4:5)]
  out = as.data.table(out)
  out[,2] = as.numeric(out$POS)
  
  data = f@gt[,ncol(f@gt)]
  data = str_split(data, ":")
  data = unlist(data)
  data = data[str_detect(data, "%")]
  
  data = str_split(data, ",")
  data = unlist(data)
  data = data[str_detect(data, "%")]
  
  data = str_remove(data, "%")
  data = as.numeric(data)
  data = as.data.table(data)
  
  uniq = unique(out[,1:3])
  uniq = transpose(uniq)
  
  out = cbind(out, data)
  
  out = lapply(uniq, check_CHROM_POS_ALT, out)
  out = rbindlist(out)
  
  temp = cbind(out[,1:3], data.table(end_position = (out$start_position + 1)))
  out = cbind(temp, out[,4:ncol(out)])
  
  return(out)
}