check_CHROM_POS = function(vec, matr){
  matr = matr[which(matr$CHROM == as.character(vec[1])), ]
  matr = matr[which(matr$POS == as.numeric(vec[2])), ]
  data = max(matr[, data])
  
  output =  data.table(ID = vec[3], chromosome_name = as.character(vec[1]), start_position = as.numeric(vec[2]), INDEL = data)
  
  return(output)
}