getGenomicFeatures = function(id, overlaps, features){
  library(GenomicFeatures)
  library(data.table)
  library(stringr)
  
  overlaps = overlaps[which(overlaps@to == id), ]
  features = features[overlaps@from, ]
  
  t_gene_id = unique(unlist(features$feature_by))
  
  if(length(t_gene_id) == 0){
    out = data.table(Gene_id = "intergenic", Genomic_feature = "NA")
    return(out)
  }
  
  t_gene_id = t_gene_id[!(t_gene_id == "")]
  
  t_gene_feature = character(length(t_gene_id))
  features_genes = unlist(features$feature_by)
  
  for(i in 1:length(t_gene_id)){
    
    features_gr = features[which(features_genes == t_gene_id[i]), ]
    features_gr = unique(unlist(features_gr$featuretype))
    features_gr = paste(features_gr, collapse = "")
    
    t_gene_feature[i] = ifelse(str_detect(features_gr, "threeUTR"), "UTR3",
                               ifelse(str_detect(features_gr, "fiveUTR"), "UTR5",
                                      ifelse(str_detect(features_gr, "cds"), "cds",
                                             ifelse(str_detect(features_gr, "exon"), "exon",
                                                    ifelse(str_detect(features_gr, "intron"), "intron", NA)))))
  }
  
  out = data.table(Gene_id = paste(t_gene_id, collapse = "|"), Genomic_feature = paste(unlist(t_gene_feature), collapse = "|"))
  
  return(out)
}