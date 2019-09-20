chrProcessing = function(chr, exout){
  
  exout = exout[which(exout$chromosome_name == chr), ]

  library(data.table)

  seq = sort(c(exout$start_position, exout$end_position))
  seq = unique(seq)
  bioranges = transpose(data.table(start = seq[1:(length(seq) - 1)], end = seq[2:length(seq)]))
  # bioranges = transpose(bioranges)
  
  out = lapply(bioranges, rangeProcessing, exout)
  out = rbindlist(out)
  out = as.data.table(out)

  return(out)
}