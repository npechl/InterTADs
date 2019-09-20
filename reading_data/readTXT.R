readTXT = function(biodata){
  library(data.table)
  
  out = read.table(as.character(biodata[1]), colClasses = c("vector"))
  
  out = out[,c(1,ncol(out))]
  colnames(out) = c("ID", as.character(out[1,2]))
  # out[,2:ncol(out)] = as.numeric(out[,2:ncol(out)])
  out = as.data.table(out)
  out = out[2:nrow(out), ]
  
  len = length(biodata)
  
  for(i in c(2:len)){
    
    newt = read.table(as.character(biodata[i]), colClasses = c("vector"))
    
    newt = newt[,c(1,ncol(newt))]
    colnames(newt) = c("ID", as.character(newt[1,2]))
    # newt[,2:ncol(newt)] = as.numeric(newt[,2:ncol(newt)])
    newt = as.data.table(newt)
    newt = newt[2:nrow(newt), ]
    
    sameIds1 = out[which(out$ID %in% newt$ID), ]
    otherIds1 = out[!(out$ID %in% newt$ID), ]
    
    sameIds2 = newt[which(newt$ID %in% out$ID), ]
    otherIds2 = newt[!(newt$ID %in% out$ID), ]
    
    sameIds1 = sameIds1[order(sameIds1$ID), ]
    sameIds2 = sameIds2[order(sameIds2$ID), ]
    
    out = cbind(sameIds1, sameIds2[ ,2])
    
    nas = matrix(nrow = nrow(otherIds1), ncol = 1)
    nas = as.data.table(nas)
    
    otherIds1 = cbind(otherIds1, nas)
    colnames(otherIds1) = colnames(out)
    
    nas = matrix(nrow = nrow(otherIds2), ncol = (ncol(out) - 2))
    nas = as.data.table(nas)
    
    temp = cbind(otherIds2[,1], nas)
    otherIds2 = cbind(temp, otherIds2[,2])
    colnames(otherIds2) = colnames(out)
    
    out = rbind(out, otherIds1)
    out = rbind(out, otherIds2)
  }
  
  for(i in 2:ncol(out)){
    out[, i] = as.numeric(out[,get(names(out)[i])])
  }
  
  out = list(out)
  
  return(out)
}