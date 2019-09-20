getAttConnection = function(matIndex, data, bioCmb, vcfmatrix = 0){
  library(data.table)
  library(stringr)
  
  data = data[[matIndex - vcfmatrix]]
  
  if(!is.null(bioCmb)){
    bioCmb = bioCmb[which(bioCmb$matrix == matIndex),]
    
    out = data.table(ID = data[,1])
    colnames(out) = "ID"
    
    for(i in 1:nrow(bioCmb)){
      columns = bioCmb[i, "cols"]
      operation = bioCmb[i, "operation"]
      weights = bioCmb[i, "weights"]
      newNames = bioCmb[i, "newNames"]
      
      columns = strsplit(as.character(columns), "; ")
      columns = unlist(columns)
      columns = str_trim(columns, side = "both")
      
      operation = as.character(operation)
      operation = str_trim(operation, side = "both")
      if(is.na(operation)){
        operation = " "
      }
      
      if(!is.na(weights)){
        weights = strsplit(as.character(weights), "; ")
        weights = unlist(weights)
        weights = str_trim(weights, side = "both")
        weights = as.numeric(weights)
      }
      
      newNames = strsplit(as.character(newNames), "; ")
      newNames = unlist(newNames)
      newNames = str_trim(newNames, side = "both")
      
      if(is.numeric(weights) && operation != "divide"){
        newData = transpose(weights * transpose(data[ ,..columns]))
      } else if(is.numeric(weights) && operation == "divide"){
        newData = data.table(col1 = (data[ ,..columns] * weights[1]), col2 = (data[ ,..columns] * weights[2]))
      } else if(!(is.numeric(weights)) && operation == "divide"){
        newData = data.table(col1 = (data[ ,..columns] * 0.5), col2 = (data[ ,..columns] * 0.5))
      } else {
        newData = data[ ,..columns]
      }
      
      if(operation == "avg"){
        newData = data.table(col = rowMeans(newData, na.rm = TRUE))
        colnames(newData) = newNames
        out = cbind(out, newData)
      } else if(operation == "sum"){
        newData = data.table(col = rowSums(newData, na.rm = TRUE))
        colnames(newData) = newNames
        out = cbind(out, newData)
      } else if(operation == "divide"){
        colnames(newData) = newNames
        out = cbind(out, newData)
      } else if(operation == "omit"){
        next
      } else {
        newData = data.table(col = newData)
        colnames(newData) = newNames
        out = cbind(out, newData)
      }
    }
  } else {
    out = data
  }
  
  # if(vcfmatrix == 0){
  #   # colScales = colSums(out[,2:ncol(out)])
  #   colScales = sapply(out[,2:ncol(out)], max, na.rm = TRUE)
  #   colScales = as.vector(colScales)
  #   # maxScale = max(colScales)
  #   numericData = out[,2:ncol(out)]
  #   numericData = mapply("*", numericData, (scale / colScales))
  #   out = cbind(out[,1], numericData)
  # }
  
  return(out)
}