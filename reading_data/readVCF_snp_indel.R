readVCF_snp_indel = function(names){
  library(data.table)
  library(stringi)
  library(stringr)
  
  indel = readVCFindel(names[1])
  snp = readVCFsnp(names[2])
  
  name = names[1]
  name = str_split(name, "/")
  name = unlist(name)
  name = name[str_detect(name, ".vcf")]
  name = str_remove(name, ".vcf")
  name = str_split(name, "indel")
  name = unlist(name)
  name = paste(name, collapse = "")
  name = str_split(name, "snp")
  name = unlist(name)
  name = paste(name, collapse = "")
  name = str_remove(name, "[.]")
  
  
  chrpos_indel = paste(indel$chromosome_name, indel$start_position, indel$end_position, sep = "")
  chrpos_snp   = paste(snp$chromosome_name, snp$start_position, snp$end_position, sep = "")
  
  sameindel = indel[which(chrpos_indel %in% chrpos_snp), ]
  samesnp   = snp[which(chrpos_snp %in% chrpos_indel), ]
  
  sameindel = sameindel[order(sameindel$chromosome_name, sameindel$start_position, sameindel$end_position), ]
  samesnp   = samesnp[order(samesnp$chromosome_name, samesnp$start_position, samesnp$end_position), ]

  otherindel = indel[!(chrpos_indel %in% chrpos_snp), ]
  othersnp   = snp[!(chrpos_snp %in% chrpos_indel), ]
  
  out = cbind(samesnp, sameindel[ ,5])
  # t   = c("A", "G", "T", "C", "INDEL")
  # colnames(out)[5:ncol(out)] = t
  
  nas = matrix(nrow = nrow(otherindel), ncol = 4)
  nas = as.data.table(nas)
  
  temp = cbind(otherindel[,1:4], nas)
  otherindel = cbind(temp, otherindel[ ,5])
  colnames(otherindel) = colnames(out)
  
  nas = matrix(nrow = nrow(othersnp), ncol = 1)
  nas = as.data.table(nas)
  
  othersnp = cbind(othersnp, nas)
  colnames(othersnp) = colnames(out) 
  
  out = rbind(out, otherindel)
  out = rbind(out, othersnp)
  
  temp  = matrix(data = c("A", "G", "T", "C", "INDEL"), nrow = (5 * nrow(out)), ncol = 1)
  temp  = as.data.table(temp)
  temp2 = matrix(nrow = (5 * nrow(out)), ncol = 1)
  temp2 = as.data.table(temp2)
  temp  = cbind(temp, temp2)
  rm(temp2)
  
  colnames(temp) = c("VarAnnotation", "second")
  
  temp$second = as.numeric(temp$second)
  temp[which(temp$VarAnnotation == "A"), ]$second = out$A
  temp[which(temp$VarAnnotation == "G"), ]$second = out$G
  temp[which(temp$VarAnnotation == "T"), ]$second = out$T
  temp[which(temp$VarAnnotation == "C"), ]$second = out$C
  temp[which(temp$VarAnnotation == "INDEL"), ]$second = out$INDEL
  
  colnames(temp) = c("VarAnnotation", name)
  
  out = out[,1:4]
  out = out[rep(seq_len(nrow(out)), each=5), ]
  
  out = cbind(out, temp)
  
  return(out)
}