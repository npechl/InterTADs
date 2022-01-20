# Loading libraries ------------------------------------------------------------

rm(list = ls())

source("R/libraries.R")

start_tad_time <- Sys.time()

# Inputs -----------------------------------------------------------------------

#' Input parameters for TADiff part
#' 
#' @param dir_name Directory of input datasets containing feature counts and 
#' frequency tables
#' 
#' @param output_folder Folder name for printing output tables
#' 
#' @param meta meta-data file name used
#' 
#' @param names.meta meta data columns to process (names or indexes)
#' 
#' @param expr_data Parent index of expression data. 
#' If no expression is provided, place FALSE

dir_name <- "Datasets"

output_folder <- "results_bloodcancer"

meta <- "meta-data.csv"

names.meta <- c('groups')
# names.meta = c("IGHV", 
#                "gain2p25.3",
#                "del8p12",
#                "gain8q24",
#                "del9p21.3",
#                "del11q22.3",
#                "trisomy12",
#                "del13q14_any",
#                "del13q14_bi",
#                "del13q14_mono",
#                "del14q24.3",
#                "del15q15.1",
#                "del17p13",
#                "Chromothripsis",
#                "BRAF",
#                "KRAS",
#                "MYD88",
#                "NOTCH1",
#                "SF3B1",
#                "TP53",
#                "ACTN2",
#                "ATM",
#                "BIRC3",
#                "CPS1",
#                "EGR2",
#                "FLRT2",
#                "IRF2BP2",
#                "KLHL6",
#                "LRP1",
#                "MED12",
#                "MGA",
#                "MUC16",
#                "NFKBIE",
#                "PCLO",
#                "UMODL1",
#                "XPO1",
#                "ZC3H18")

expr_data <- 2



data.all <- fread(paste(output_folder, "/integrated-tad-table-methNorm.txt",
                        sep = ""),
                    sep = "\t")

data.all$ID <- paste(data.all$tad_name, data.all$ID, sep = ";")

meta <- fread(paste(dir_name, meta, sep = "/"))
who <- meta == ""
who <- apply(who, 1, sum, na.rm = TRUE)
meta <- meta[which(who == 0), ]


sample.list <- meta$newNames
data.all$Mean <- rowMeans(data.all[,..sample.list])


data.all <- data.all[which(round(data.all$Mean) > 10), ]


sign_table <- matrix(0, nrow = length(names.meta), ncol = 2)

rm(who)

for (z in 1:length(names.meta)){
  
    cat(c(names.meta[z], "\n"))
    
    analysis <- names.meta[z]
    groups <- as.character(meta[[analysis]])
    groups <- unique(groups)
    groups <- groups[which(groups != "")]
    groups <- groups[!is.na(groups)]
    
    group1 <- groups[1]
    group2 <- groups[2]
    
    meta.keep <- meta[which(meta[[analysis]] == group1 | meta[[analysis]] == 
                                group2), ]
    
    sample.list <- meta.keep$newNames
    
    df <- data.all[,..sample.list]
    
    df <- as.data.frame(df)
    row.names(df) <- data.all$ID
    
    pheno <- as.factor(meta.keep[[analysis]])
    
    phenoMat <- model.matrix(~pheno)
    colnames(phenoMat) <- sub("^pheno", "", colnames(phenoMat))
    
    fit <- lmFit(object = df, design = phenoMat)
    
    gc()
    set.seed(6)
    fit <- eBayes(fit)
    
    gc()
    top.rank <- topTable(fit, number = nrow(df), adjust.method = "fdr",
                        sort.by = "p")
    
    sign.table <- top.rank[which(top.rank$adj.P.Val <= 0.01 & 
                                    abs(top.rank$logFC) > 2), ]
    
    if (nrow(sign.table) == 0) {
      
        cat(c("No statistical significant events for:", names.meta[z], "\n"))
        
        sign_table[z,1] <- analysis 
        sign_table[z,2] <- "0" 
    }
      
    else {
      
        # annotate sign.table
        
        sign.table$ID <- row.names(sign.table)
        
        sign.table <- merge(sign.table, 
                            data.all, 
                            by.x = "ID", 
                            by.y = "ID")
        
        sign.tad.info <- sign.table %>% 
        dplyr::group_by(tad_name) %>% 
        dplyr::summarise(count = n()) 
        
        
        
        # annotate tad.info table
        
        tad.info <- data.all %>% 
        dplyr::group_by(tad_name) %>% 
        dplyr::summarize(count = n()) 
        
        sign.tad.info <- merge(sign.tad.info, 
                                tad.info, 
                                by = "tad_name")
        
        sign.tad.info$freq <- sign.tad.info$count.x /
            sign.tad.info$count.y * 100
      
        sign.tad.info$pvalue <- 1
        
        for (i in 1:nrow(sign.tad.info)) {
          
            sign.tad.info$pvalue[i] <- 1 - phyper(sign.tad.info$count.x[i], 
                                                nrow(sign.table), 
                                                nrow(df) - nrow(sign.table), 
                                                sign.tad.info$count.y[i])
          
        }  
        
        
        if(expr_data != FALSE){
          
            # get expression data
            
            expr <- data.all[which(data.all$parent == expr_data), ]
            
            df <- expr[,..sample.list]
            
            df <- as.data.frame(df)
            row.names(df) <- expr$ID
            
            # build model
            
            fit <- lmFit(object = df, design = phenoMat)
            
            gc()
            set.seed(6)
            fit <- eBayes(fit)
            
            # get top rank events
            
            gc()
            top.rank <- topTable(fit, number = nrow(df), adjust.method = "fdr")
            
            # any filtering ?
            
            # annotate sign.table (expression)
            
            sign.table.expr <- as.data.frame(top.rank)
            
            sign.table.expr$ID <- row.names(sign.table.expr)
            
            sign.table.expr <- merge(sign.table.expr, 
                                    expr, 
                                    by = "ID")
            
            sign.tad.expr <- sign.table.expr %>% 
            dplyr::group_by(tad_name) %>% 
            dplyr::summarise(mean_logFC = mean(abs(logFC)))
            
            
            
            # merge information
            
            tad.all.info <- merge(sign.tad.info, 
                                    sign.tad.expr, 
                                    by = "tad_name" )
            
            tad.all.info.f <- tad.all.info %>% 
            filter(count.x > 4) %>% 
            filter(pvalue < 0.01) %>% 
            filter(freq > 15) %>% 
            filter(mean_logFC > 2)
            
            sign_table[z,1] <- analysis 
            sign_table[z,2] <- as.character(nrow(tad.all.info.f))
            
            
            
            # create output tables
            
            full.tads <- merge(tad.all.info.f, 
                                data.all, 
                                by = "tad_name")
            
            if(nrow(full.tads) > 0) {
              
                write.table(full.tads, 
                            file = paste(output_folder, "/", analysis,
                                            "_TADiff.txt", sep = ""), 
                            col.names = TRUE, 
                            row.names = FALSE, 
                            quote = FALSE, 
                            sep = "\t")
              
            }
        }
        
        # what's the value of `sign_table` if we don't have expression data
        # which are the outputs if we don't have expression data
    } 
}

sign_table <- as.data.frame(sign_table)

write.table(sign_table, 
            file = paste(output_folder, "/Summary_TADiff.txt", sep = ""), 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")


