createScale = function(data, dataVCF, scale){
  # dataVCF has by default maximum value = 100 
  # matrixScales = c()
  # 
  # for(temp in data){
  #   matrixScales = cbind(matrixScales, max(colSums(temp[,2:ncol(temp)])))
  # }
  # 
  # matrixScales = c(matrixScales, 100)
  # 
  # maxScale = max(matrixScales)
  
  i = 0
  
  for(temp in data){
    i = i + 1
    # tempScale = max(temp[,2:ncol(temp)], na.rm = TRUE)
    # 
    # numericData = temp[,2:ncol(temp)]
    # numericData = numericData * (scale / tempScale)
    # temp = cbind(temp[,1], numericData)
    colScales = sapply(temp[,2:ncol(temp)], max, na.rm = TRUE)
    colScales = as.vector(colScales)
    numericData = temp[,2:ncol(temp)]
    numericData = mapply("*", numericData, (scale / colScales))
    temp = cbind(temp[,1], numericData)
    data[[i]] = temp
  }
  
  i = 0
  
  for(temp in dataVCF){
    i = i + 1

    # numericData = temp[,6:ncol(temp)]
    # numericData = numericData * (scale / 100)
    # temp = cbind(temp[,1:5], numericData)
    
    colScales = sapply(temp[,6:ncol(out)], max, na.rm = TRUE)
    colScales = as.vector(colScales)
    numericData = temp[,6:ncol(temp)]
    numericData = mapply("*", numericData, (scale / colScales))
    temp = cbind(temp[,1:5], numericData)
    dataVCF[[i]] = temp
  }
  
  out = list(data, dataVCF)
  
  return(out)
}