getBioLocation = function(Biovector, database = "ensembl", dataSet="hsapiens_gene_ensembl"){
  
  #------------------------------ Libraries ------------------------------
  
  if(!require("biomaRt")){
    print("biomaRt package is required")
    return(Biovector)
  } else {
    library(biomaRt)
    biomart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
  }
  
  if(!require("IlluminaHumanMethylation450kanno.ilmn12.hg19")){
    print("IlluminaHumanMethylation450kanno.ilmn12.hg19 is required")
    return(Biovector)
  } else {
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    data(Locations)
  }
  
  library(stringi)
  library(stringr)
  library(data.table)
  
  #------------------------------ Creating chromo data table ------------------------------
  
  Biovector = as.character(Biovector)
  
  output = data.table(ID = character(), chromosome_name = character(), start_position = numeric(), end_position = numeric(), VarAnnotation = character())
  
  teV = Biovector[str_detect(Biovector, "cg")]
  
  if(length(teV) != 0){
    # Methylation data
    
    chromo = Locations[teV, 1:2]
    chromo = data.table(ID = teV, chromosome_name = as.integer(str_remove(chromo$chr, "chr")), 
                                          start_position = (chromo$pos - 1), 
                                          end_position = chromo$pos)
    
    temp = matrix(nrow = nrow(chromo), ncol = 1)
    temp = as.data.table(temp)
    colnames(temp) = c("VarAnnotation")
    chromo = cbind(chromo, temp)
    
    output = rbind(output, chromo)
    
  } 
  
  teV = Biovector[str_detect(Biovector, "NR")]
  
  if(length(teV) != 0){
    # RefSeq ncRNA data
    
    refseqids = teV
    
    chromo = getBM(attributes = c("refseq_ncrna","chromosome_name", "start_position", "end_position"), 
                      filters = "refseq_ncrna",
                       values = refseqids, 
                         mart = biomart)
    
    colnames(chromo) = c("ID", "chromosome_name", "start_position", "end_position")
    
    chromo = as.data.table(chromo)
    
    temp = matrix(nrow = nrow(chromo), ncol = 1)
    temp = as.data.table(temp)
    colnames(temp) = c("VarAnnotation")
    chromo = cbind(chromo, temp)
    
    output = rbind(output, chromo)
  } 
  
  teV = Biovector[str_detect(Biovector, "NM")]
  
  if(length(teV) != 0){
    # RefSeq mRNA data
    
    refseqids = teV
    
    chromo = getBM(attributes = c("refseq_mrna","chromosome_name", "start_position", "end_position"), 
                      filters = "refseq_mrna",
                       values = refseqids, 
                         mart = biomart)
    
    colnames(chromo) = c("ID", "chromosome_name", "start_position", "end_position")
    
    chromo = as.data.table(chromo)
    
    temp = matrix(nrow = nrow(chromo), ncol = 1)
    temp = as.data.table(temp)
    colnames(temp) = c("VarAnnotation")
    chromo = cbind(chromo, temp)
    
    output = rbind(output, chromo)
  }
  
  teV = Biovector[str_detect(Biovector, "NP")]
  
  if(length(teV) != 0){
    # RefSeq peptide data
    
    refseqids = teV
    
    chromo = getBM(attributes=c("refseq_peptide","chromosome_name", "start_position", "end_position"), 
                   filters="refseq_peptide",
                   values=refseqids, 
                   mart=biomart)
    
    colnames(chromo) = c("ID", "chromosome_name", "start_position", "end_position")
    
    chromo = as.data.table(chromo)
    
    temp = matrix(nrow = nrow(chromo), ncol = 1)
    temp = as.data.table(temp)
    colnames(temp) = c("VarAnnotation")
    chromo = cbind(chromo, temp)
    
    output = rbind(output, chromo)
  }
  
  teV = Biovector[str_detect(Biovector, "ENSG")]
  
  if(length(teV) != 0){
    # RefSeq peptide data
    
    ensembl_gene_ids = teV
    
    chromo = getBM(attributes=c("ensembl_gene_id","chromosome_name", "start_position", "end_position"), 
                   filters="ensembl_gene_id",
                   values=ensembl_gene_ids, 
                   mart=biomart)
    
    colnames(chromo) = c("ID", "chromosome_name", "start_position", "end_position")
    
    chromo = as.data.table(chromo)
    
    temp = matrix(nrow = nrow(chromo), ncol = 1)
    temp = as.data.table(temp)
    colnames(temp) = c("VarAnnotation")
    chromo = cbind(chromo, temp)
    
    output = rbind(output, chromo)
  }
  
  temp = output$end_position - output$start_position
  temp = cbind(output, data.table(range = temp))
  temp = temp[order(temp$range, decreasing = TRUE), ]
  output = temp[,1:5]
  output = output[!duplicated(output$ID), ]
  
  return(output)
}