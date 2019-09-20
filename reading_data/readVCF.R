readVCF = function(loc){
  library(data.table)
  library(stringi)
  library(stringr)
  
  out = matrix(nrow = 0, ncol = 5)
  out = as.data.table(out)
  colnames(out) = c("ID", "chromosome_name", "start_position", "end_position", "VarAnnotation")
  
  for(i in seq(1, length(loc), 2)){
    temp = readVCF_snp_indel(loc[i:(i + 1)])
    
    index_out = paste(out$ID, out$chromosome_name, out$start_position, out$end_position, out$VarAnnotation, sep = "")
    index_temp = paste(temp$ID, temp$chromosome_name, temp$start_position, temp$end_position, temp$VarAnnotation, sep = "")
    
    sameOut = out[which(index_out %in% index_temp), ]
    sameTemp = temp[which(index_temp %in% index_out), ]
    
    sameOut = sameOut[order(sameOut$chromosome_name, sameOut$start_position, sameOut$end_position, sameOut$VarAnnotation), ]
    sameTemp = sameTemp[order(sameTemp$chromosome_name, sameTemp$start_position, sameTemp$end_position, sameTemp$VarAnnotation), ]
    
    otherOut = out[!(index_out %in% index_temp), ]
    otherTemp = temp[!(index_temp %in% index_out), ]
    
    out = cbind(sameOut, sameTemp[ ,6])
    
    nas = matrix(data = numeric(), nrow = nrow(otherTemp), ncol = (ncol(otherOut) - 5))
    nas = as.data.table(nas)
    temp2 = cbind(otherTemp[,1:5], nas)
    otherTemp = cbind(temp2, otherTemp[,6])
    colnames(otherTemp) = colnames(out)
    
    nas = matrix(data = numeric(), nrow = nrow(otherOut), ncol = 1)
    nas = as.data.table(nas)
    otherOut = cbind(otherOut, nas)
    colnames(otherOut) = colnames(out)
    
    out = rbind(out, otherOut)
    out = rbind(out, otherTemp)
  }
  
  nas = out[,6:ncol(out)]
  nas = is.na(nas)
  nas = as.data.table(nas)

  nas = rowSums(nas)
  nas = (nas == ncol(out[,6:ncol(out)]))

  out = out[!nas, ]
  
  out$chromosome_name = str_remove(out$chromosome_name, "chr")
  
  out = list(out)
  
  return(out)
}