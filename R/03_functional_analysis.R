########## Loading libraries ########## 

source("R/libraries_enrich.R")
source("R/motif_enrich.R")
source("R/go_pathway_enrich.R")
source("R/visualization_enrich.R")

########### Inputs ########## 


#' Input parameters for TADiff part
#' 
#' @param tech Human Genome Reference used
#' 
#' @param dbs Databases used from EnrichR
#' 
#' @param type the prevously selected databases acronyms used for the names of the outputs files
#'             e.g. GO.MF for GO_Molecular_Function_2018
#' 
#' @param genes.cover gene coverage of the databases
#' 
#' @param p.adjust.method p adjustment method, the methods supported are: 
#'                        c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
#' 
#' @param criterio Enrichr result column selected as criterio: "P.value" or "Adjusted.P.value"
#' 
#' @param cut.off cut-off Enrichr enrichment (adjusted) p-value
#' 
#' @param cut.off.TF cut-off motif enrichment adjusted p-value
#' 
#' @param min.genes min number of genes in over-represented terms
#' 
#' @param system RStudio supports different fonts for different operating systems
#'  
#' @param dir_name name or filepath of the input folder
#'  
#' @param output_folder name or filepath of the output folder  
#'  
#' 

tech <- "hg19" # or "hg38"

dbs <- c("GO_Molecular_Function_2018","GO_Biological_Process_2018","KEGG_2019_Human")

type <- c("GO.MF","GO.BP","KEGG")

genes.cover <- c(11459,14433,7802)

#choose a p adjust method, the methods supported are:
#c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
p.adjust.method <- "fdr"

cut.off.TF <- 0.05

cut.off <- 0.05

criterio <- "P.value"  #"Adjusted.P.value" 

min.genes <- 3

system <- "win"

dir_name <- "test_files"

output_folder <- paste0("Outputs_test_",criterio,"_",cut.off)

#scal_test <- data.table(part = character(6),
#                        time = character(6),
#                        memory = numeric(6) )

#scal_test$part[1] <- "start"
#scal_test$time[1] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")
#scal_test$memory[1] <- mem_used()

########### Enrichment + Data Analysis ##########

dir.create(output_folder, showWarnings = FALSE)

data.type <- data.table(type =c("GO.MF","GO.BP","KEGG"),
                        cover = as.numeric(genes.cover))

files.evenDiff <- list.files(dir_name,pattern = c("evenDiff"))
files.TADiff <- list.files(dir_name,pattern = c("TADiff"))

if (!is_empty(files.evenDiff)){
   for (one.file in files.evenDiff){
      
      file.name <- str_remove(one.file,".txt")
      folder <- createFolders(paste0(output_folder,"/",file.name))
      biodata <- fread(paste0(dir_name,"/",one.file))
      
      #enrichment all
      listAll <- enrichAll(biodata,dbs, cut.off,criterio, type)
      
      #data analysis
      
      for (l in c(1:length(dbs))){
         if (nrow(listAll[[l]]) > 0){
            
            name <- names(listAll)[l]
            genes.coverage <- data.type[type == name,cover]
            result <- dataAnalysis(listAll[[l]], name,listAll$data.with.genes, genes.coverage,p.adjust.method, min.genes)
            if (!is.null(result)) {assign(paste0('list.', name, '.All'), result)}
         }
         
      }
      
      #get pathview input data
      #if (nrow(listAll$KEGG) != 0){
      
      #listPathAll <- getKEGGIds(listAll$KEGG, biodata[,Gene_id,diff])         
      #}
      
      if( exists("list.GO.MF.All") & exists("list.GO.BP.All")){
         
         #join GO Molecular Function and Biological Process outputs
         dataAll <- full_join(list.GO.MF.All$data.perTAD, list.GO.BP.All$data.perTAD, by = "TAD")
      }
      
      
      #motif enrichment
      tech <- "hg19"
      exp.parent <- 1
      report.list <- motifEnrich(biodata, folder$motifOutputsFolder,dir_name,
                                 p.adjust.method, cut.off.TF, tech, exp.parent)
      #report.list <- dget(paste0(folder$motifOutputsFolder,"/report MotifEA.txt"))
      listMotif <- motifOutputs(report.list)
      
      ########### Output Files ########## 
      
      #enrichment all
      if (exists("dataAll")){
         fwrite(dataAll, paste0(folder$goAllOutputs, "/over-represented GO terms-enrichment all.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
      }
      
      if (exists("list.GO.MF.All")){
         fwrite(list.GO.MF.All$data.perTerm, paste0(folder$goAllOutputs, "/GO MF Terms in different TADs.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
         fwrite(list.GO.MF.All$data.perTAD, paste0(folder$goAllOutputs, "/over-represented GO MF terms-enrichment all.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
      }
      
      if (exists("list.GO.BP.All")){
         
         fwrite(list.GO.BP.All$data.perTerm, paste0(folder$goAllOutputs, "/GO BP Terms in different TADs.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
         
         fwrite(list.GO.BP.All$data.perTAD, paste0(folder$goAllOutputs, "/over-represented GO BP terms-enrichment all.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
      }
      
      if (exists("list.KEGG.All")){
         
         fwrite(list.KEGG.All$data.perTAD, paste(folder$keggAllOutputs, "/over-represented KEGG Pathways-enrichment all.csv", sep = ""), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
         
         fwrite(list.KEGG.All$data.perTerm, paste0(folder$keggAllOutputs, "/KEGG Pathways in different TADs.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
         #fwrite(listPathAll$output.csv, paste0(folder$keggAllOutputs, "/Pathview input.csv"),
         #       row.names = FALSE, sep = "\t", quote = FALSE)
         
         
      }
      
     
      #motif enrichment
      fwrite(listMotif$table_perTAD, paste0(folder$motifOutputsFolder, "/over-represented TFs in each tad.csv"), 
             row.names = FALSE, sep = "\t", quote = FALSE)
      fwrite(listMotif$table_perTFs, paste0(folder$motifOutputsFolder, "/TFs in different TADs.csv"), 
             row.names = FALSE, sep = "\t", quote = FALSE)
      file.create(paste0(folder$motifOutputsFolder,"/report MotifEA.txt"), showWarnings = FALSE)
      dput(report.list, file = paste0(folder$motifOutputsFolder,"/report MotifEA.txt"))
      
      ########### Visualization ##########
      
      setGraphFonts(system)
      
      #enrich all visualization
      if (exists("list.GO.MF.All")){
         enrichrVisual(folder$goAllImages,"GO MF Terms",list.GO.MF.All$data.visual, criterio)
      }
      
      if (exists("list.GO.BP.All")){
         enrichrVisual(folder$goAllImages, "GO BP Terms",list.GO.BP.All$data.visual, criterio)
      }
      
      if (exists("list.KEGG.All")){
         enrichrVisual(folder$keggAllImages, "KEGG Pathways",list.KEGG.All$data.visual, criterio)
         #pathVisual(listPathAll$pathview.input ,folder$keggAllImages)
         
      }
      
      #motif enrichment analysis visualization
      #report.list <- dget(paste0(folder$motifOutputsFolder,"/report MotifEA.txt"))
      motifVisual(folder$motifImageOutputs, folder$motifOutputsFolder, listMotif$data.visual, report.list, criterio)
      
      rm(listMotif,listAll,list.GO.MF.All,list.GO.BP.All,list.KEGG.All, dataAll,
         listPerTAD,list.GO.MF.PerTAD,list.GO.BP.PerTAD,list.KEGG.PerTAD, dataPerTAD,
         listPathAll, listPathPerTAD)
   }
   
}

if (!is_empty(files.TADiff)){
   for (one.file in files.TADiff){
      
      file.name <- str_remove(one.file,".txt")
      folder <- createFolders(paste0(output_folder,"/",file.name))
      biodata <- fread(paste0(dir_name,"/",one.file))
      
      #enrichment all
      listAll <- enrichAll(biodata,dbs, cut.off,criterio, type)
      
      #data analysis
      
      for (l in c(1:length(dbs))){
         if (nrow(listAll[[l]]) > 0){
            
            name <- names(listAll)[l]
            genes.coverage <- data.type[type == name,cover]
            result <- dataAnalysis(listAll[[l]], name,listAll$data.with.genes, genes.coverage,p.adjust.method, min.genes)
            if (!is.null(result)) {assign(paste0('list.', name, '.PerTAD'), result)}
         }
         
      }
      
      #get pathview input data
      #if (nrow(listAll$KEGG) != 0){
      
      #listPathAll <- getKEGGIds(listAll$KEGG, biodata[,Gene_id,diff])         
      #}
      
      if( exists("list.GO.MF.All") & exists("list.GO.BP.All")){
         
         #join GO Molecular Function and Biological Process outputs
         dataAll <- full_join(list.GO.MF.All$data.perTAD, list.GO.BP.All$data.perTAD, by = "TAD")
      }
      
      #enrichment per TAD
      listPerTAD <- enrichPerTAD(biodata, dbs, cut.off,criterio, type)
         
      for (l in c(1:length(dbs))){
         if (nrow(listPerTAD[[l]]) > 0){
               
            name <- names(listPerTAD)[l]
            genes.coverage <- data.type[type == name,cover]
            result <- dataAnalysis(listPerTAD[[l]], name,listPerTAD$data.with.genes, genes.coverage,p.adjust.method, min.genes)
            if (!is.null(result)) {assign(paste0('list.', name, '.PerTAD'), result)}
         }
            
      }
         
      #get pathview input data
      #if (nrow(listAll$KEGG) != 0){
         
      #listPathPerTAD <- getKEGGIds(listPerTAD$KEGG, biodata[,Gene_id,diff])        #get pathview input data
      #}
         
      if(exists("list.GO.MF.PerTAD") & exists("list.GO.BP.PerTAD")){
            
         #join GO Molecular Function and Biological Process outputs
         dataPerTAD <- full_join(list.GO.MF.PerTAD$data.perTAD, list.GO.BP.PerTAD$data.perTAD, by = "TAD")
      }
         
      
      
      #motif enrichment
      tech <- "hg19"
      exp.parent <- 1
      report.list <- motifEnrich(biodata, folder$motifOutputsFolder,p.adjust.method, cut.off.TF, tech, exp.parent)

      #report.list <- dget(paste0(folder$motifOutputsFolder,"/report MotifEA.txt"))
      listMotif <- motifOutputs(report.list)
      
      ########### Output Files ########## 
      
      #enrichment all
      if (exists("dataAll")){
         fwrite(dataAll, paste0(folder$goAllOutputs, "/over-represented GO terms-enrichment all.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
      }
      
      if (exists("list.GO.MF.All")){
         fwrite(list.GO.MF.All$data.perTerm, paste0(folder$goAllOutputs, "/GO MF Terms in different TADs.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
         fwrite(list.GO.MF.All$data.perTAD, paste0(folder$goAllOutputs, "/over-represented GO MF terms-enrichment all.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
      }
      
      if (exists("list.GO.BP.All")){
         
         fwrite(list.GO.BP.All$data.perTerm, paste0(folder$goAllOutputs, "/GO BP Terms in different TADs.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
         
         fwrite(list.GO.BP.All$data.perTAD, paste0(folder$goAllOutputs, "/over-represented GO BP terms-enrichment all.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
      }
      
      if (exists("list.KEGG.All")){
         
         fwrite(list.KEGG.All$data.perTAD, paste(folder$keggAllOutputs, "/over-represented KEGG Pathways-enrichment all.csv", sep = ""), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
         
         fwrite(list.KEGG.All$data.perTerm, paste0(folder$keggAllOutputs, "/KEGG Pathways in different TADs.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
         
         #fwrite(listPathAll$output.csv, paste0(folder$keggAllOutputs, "/Pathview input.csv"),
         #       row.names = FALSE, sep = "\t", quote = FALSE)
         
         
      }
      
      #enrichment per TAD  
      if (exists("dataPerTAD")){
           fwrite(dataPerTAD, paste(folder$goPerOutputs, "/over-represented GO terms-enrichment per tad.csv", sep = ""), 
                row.names = FALSE, sep = "\t", quote = FALSE)
      }
         
      if (exists("list.GO.MF.PerTAD")){
         fwrite(list.GO.MF.PerTAD$data.perTerm, paste0(folder$goPerOutputs, "/GO MF Terms in different TADs.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
            
         fwrite(list.GO.MF.PerTAD$data.perTAD, paste0(folder$goPerOutputs, "/over-represented GO MF terms-enrichment per TAD.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
      
      }
         
      if (exists("list.GO.BP.PerTAD")){
            
         fwrite(list.GO.BP.PerTAD$data.perTerm, paste0(folder$goPerOutputs, "/GO BP Terms in different TADs.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
            
            
         fwrite(list.GO.BP.PerTAD$data.perTAD, paste0(folder$goPerOutputs, "/over-represented GO BP terms-enrichment per TAD.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
            
      }
         
      if (exists("list.KEGG.PerTAD")){
            
         fwrite(list.KEGG.PerTAD$data.perTAD, paste(folder$keggPerOutputs, "/over-represented KEGG Pathways-enrichment per TAD.csv", sep = ""), 
                row.names = FALSE, sep = "\t", quote = FALSE)
            
            
         fwrite(list.KEGG.PerTAD$data.perTerm, paste0(folder$keggPerOutputs, "/KEGG Pathways in different TADs.csv"), 
                row.names = FALSE, sep = "\t", quote = FALSE)
            
         #fwrite(listPathPerTAD$output.csv, paste0(folder$keggPerOutputs, "/Pathview input.csv"),
         #       row.names = FALSE, sep = "\t", quote = FALSE)
            
      }
  
      
      #motif enrichment
      fwrite(listMotif$table_perTAD, paste0(folder$motifOutputsFolder, "/over-represented TFs in each tad.csv"), 
             row.names = FALSE, sep = "\t", quote = FALSE)
      fwrite(listMotif$table_perTFs, paste0(folder$motifOutputsFolder, "/TFs in different TADs.csv"), 
             row.names = FALSE, sep = "\t", quote = FALSE)
      file.create(paste0(folder$motifOutputsFolder,"/report MotifEA.txt"), showWarnings = FALSE)
      dput(report.list, file = paste0(folder$motifOutputsFolder,"/report MotifEA.txt"))
      
      ########### Visualization ##########
      
      setGraphFonts(system)
      
      #enrich all visualization
      if (exists("list.GO.MF.All")){
         enrichrVisual(folder$goAllImages,"GO MF Terms",list.GO.MF.All$data.visual, criterio)
      }
      
      if (exists("list.GO.BP.All")){
         enrichrVisual(folder$goAllImages, "GO BP Terms",list.GO.BP.All$data.visual, criterio)
      }
      
      if (exists("list.KEGG.All")){
         enrichrVisual(folder$keggAllImages, "KEGG Pathways",list.KEGG.All$data.visual, criterio)
         #pathVisual(listPathAll$pathview.input ,folder$keggAllImages)
         
      }
      
      
      #enrich per TAD visualization 
      if (exists("list.GO.MF.PerTAD")){
         enrichrVisual(folder$goPerImages,"GO MF Terms",list.GO.MF.PerTAD$data.visual, criterio)
      }
         
      if (exists("list.GO.BP.PerTAD")){
         enrichrVisual(folder$goPerImages, "GO BP Terms",list.GO.BP.PerTAD$data.visual, criterio)
      }
         
      if (exists("list.KEGG.PerTAD")){
         enrichrVisual(folder$keggPerImages, "KEGG Pathways",list.KEGG.PerTAD$data.visual, criterio)
         #pathVisual(listPathPerTAD$pathview.input ,folder$keggPerImages)         }
         
      }
      
      
      #motif enrichment analysis visualization
      #report.list <- dget(paste0(folder$motifOutputsFolder,"/report MotifEA.txt"))
      motifVisual(folder$motifImageOutputs, folder$motifOutputsFolder, listMotif$data.visual, report.list, criterio)
      
      rm(listMotif,listAll,list.GO.MF.All,list.GO.BP.All,list.KEGG.All, dataAll,
         listPerTAD,list.GO.MF.PerTAD,list.GO.BP.PerTAD,list.KEGG.PerTAD, dataPerTAD,
         listPathAll, listPathPerTAD)
   }
   
}
