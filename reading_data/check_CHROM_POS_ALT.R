check_CHROM_POS_ALT = function(vec, matr){
  matr = matr[which(matr$CHROM == as.character(vec[1])), ]
  matr = matr[which(matr$POS == as.numeric(vec[2])), ]
  
  temp = matrix(nrow = 1, ncol = 4)
  temp = as.data.table(temp)
  colnames(temp) = c("A", "G", "T", "C")
  
  for(i in 1:nrow(matr)){
    if(matr[i, ALT] == "T"){
      temp$T = max(temp$T, matr[i, data], na.rm = TRUE)
    } else if(matr[i, ALT] == "G"){
      temp$G = max(temp$G, matr[i, data], na.rm = TRUE)
    } else if(matr[i, ALT] == "C"){
      temp$C = max(temp$C, matr[i, data], na.rm = TRUE)
    } else {
      temp$A = max(temp$A, matr[i, data], na.rm = TRUE)
    }
  }
  
  output =  data.table(ID = vec[3], chromosome_name = as.character(vec[1]), start_position = as.numeric(vec[2]))
  output = cbind(output, temp)
  
  return(output)
}