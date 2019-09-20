rangeProcessing = function(biorange, exout){
  connectors = exout[which(biorange[1] >= exout$start_position), ]
  connectors = connectors[which(biorange[2] <= connectors$end_position), ]
  temp = connectors[ ,8:ncol(connectors)]
  
  numOfRows = nrow(connectors)
  
  out = data.table(ID = character(numOfRows), 
                   chromosome_name = character(numOfRows), 
                   start_position = numeric(numOfRows), 
                   end_position = numeric(numOfRows),
                   VarAnnotation = character(numOfRows),
                   Gene_id = character(numOfRows),
                   Genomic_feature = character(numOfRows))
  
  if(numOfRows != 0){
    
    out$ID = connectors$ID
    out$chromosome_name = connectors$chromosome_name
    out$start_position = biorange[1]
    out$end_position = biorange[2]
    out$VarAnnotation = connectors$VarAnnotation
    out$Gene_id = connectors$Gene_id
    out$Genomic_feature = connectors$Genomic_feature
    
    scl = (biorange[2] - biorange[1]) / (connectors$end_position - connectors$start_position)
    
    # out = rbind(out, data.table(ID = connectors$ID, chromosome_name = connectors$chromosome_name,
    #                             start_position = biorange[1], end_position = biorange[2]))
    
    temp = scl * temp
    
    # out = data.table(ID = connectors$ID, chromosome_name = connectors$chromosome_name,
    #                  start_position = biorange[1], end_position = biorange[2])
  }
  
  out = cbind(out, temp)
  
  return(out)
}