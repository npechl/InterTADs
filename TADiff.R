########## Loading libraries ########## 

source("libraries.R")

start_tad_time = Sys.time()

########### Inputs ###########

#' Input parameters for TADiff part
#' 
#' @param dir_name Directory of input datasets containing feature counts and frequency tables
#' 
#' @param output_folder Folder name for printing output tables
#' 
#' @param image_output_folder Folder name for print output images
#' 
#' @param tad_file BED file containing information about TADs
#' 
#' @param meta meta-data file name used
#' 
#' @param paired.data Boolean flag indicating whether input data is paired or not
#' 
#' @param FDR_criterion user defined FDR criterion

dir_name = "newDatasets"
output_folder = "onlyMethylation_newDatasets"
image_output_folder = "onlyMethylation_newDatasets_vis"

tad_file = "hglft_genome_2dab_ec1330.bed"

meta = "meta-data-file.txt"

paired.data = FALSE

FDR_criterion = 0.05

#'
#' Reading inputs
#' 


data.all = fread(paste(output_folder, "/integrated_table.csv", sep = ""))

TAD = fread(paste(dir_name, tad_file, sep = "/"), 
            header = F, 
            sep = "\t", 
            check.names = FALSE)

meta = fread(paste(dir_name, meta, sep = "/"))
who = meta == ""
who = apply(who, 1, sum, na.rm = TRUE)
meta = meta[which(who == 0), ]


#' 
#' @param groups used to compare two groups from the meta data file
#' in case there is no need for any comparison to perform, 
#' please set the variable to NULL:
#' 
#' groups = NULL
#' 
#' 

diff_col = "groups"
groups = meta[[diff_col]]
groups = as.character(groups)
groups = unique(groups)
groups = groups[which(groups != "")]
groups = groups[which(!is.na(groups))]
group1 = groups[1]
group2 = groups[2]

########### Processing ########### 

#'
#' Generate a barcode
#'  
data.all$name = paste0(data.all$ID, data.all$Gene_id)

#' Reordering
#' 
data.over = data.all[,c("name", "chromosome_name", "start_position", "end_position")]
data.over$chromosome_name = paste("chr", data.over$chromosome_name, sep = "")
data.over = data.over[,c(2,3,4,1)]

#'
#' Overlap of the TAD with events
#' 
colnames(TAD) = paste(c("chr", "start", "end", "name"))
colnames(data.over) = paste(c("chr", "start", "end", "name"))

#'
#' Make GRanges object
#' 
gr1 = with(TAD, GRanges(chr, IRanges(start = start, end = end, names = name)))
gr2 = with(data.over, GRanges(chr, IRanges(start = start, end = end, names = name)))

#'
#' Completely overlapping
#' 
type1 = findOverlaps(query = gr1, subject = gr2)

type1.df = cbind(TAD[queryHits(type1),], data.over[subjectHits(type1),])
type1.df = type1.df[,c(2,3,4,8)]
colnames(type1.df) = c("tad_start", "tad_end", "tad_name", "name.1")

#'
#' Generating output overlapping with
#' TADs table
#' 
full = base::merge(type1.df, data.all, by.x = "name.1", by.y = "name")
keep = c("chromosome_name", 
         "tad_name", 
         "tad_start", 
         "tad_end", 
         "ID", 
         "start_position", 
         "end_position", 
         "Gene_id", 
         "Gene_locus", 
         "parent",
         meta$newNames)

full = full[,..keep]

# test for NAs, remove them for the statistical analysis
# full[is.na(full)] <- 0
# full$Gene_id[which(full$Gene_id == 0)] = "NA"

########### Statistics ########### 

if(!is.null(groups)){
  
  #'
  #' Statistical differences between two groups
  #' 
  
  list1 = meta[which(meta[[diff_col]] == group1), ]$newNames
  list2 = meta[which(meta[[diff_col]] == group2), ]$newNames
  
  full.1 = full[,..list1]
  full.2 = full[,..list2]
  
  full$diff = apply(full.2, 1, mean, na.rm = TRUE) - apply(full.1, 1, mean, na.rm = TRUE)
  
  #'
  #' Exclude the events with no differences
  #' 
  full = full[which(round(abs(full$diff))>0),]
  
  #'
  #' Statistics
  #' 
  tad_sum = dplyr::group_by(full, tad_name) %>% 
    dplyr::summarise(count = n(), 
                     mean = mean(abs(diff), na.rm = TRUE),
                     IQR = IQR(diff, na.rm = TRUE))
  
  tad_sum$ttest = numeric(length = nrow(tad_sum))
  tad_sum$wilcoxon = numeric(length = nrow(tad_sum))
  tad_sum = as.data.table(tad_sum)
  
  tad_sum$ttest = NA
  tad_sum$wilcoxon = NA
  
  #'
  #' parametric t.test
  #' 
  #' non-parametric wilcoxon
  #' 
  
  if(paired.data){
    
    for (i in 1:nrow(tad_sum)){
      tad = full[which(full$name == as.character(tad_sum[i,1])), ]
      
      tad.1 = tad[,..list1]
      tad.2 = tad[,..list2]
      
      statistics2 = t.test(unlist(tad.1), unlist(tad.2), na.rm = TRUE, paired = T)
      statistics3 = wilcox.test(unlist(tad.1), unlist(tad.2), na.rm = TRUE, paired = T)
      
      tad_sum[i,5] = statistics2$p.value
      tad_sum[i,6] = statistics3$p.value
    }
    
  } else {
    
    for (i in 1:nrow(tad_sum)){
      tad = full[which(full$tad_name == as.character(tad_sum[i,1])), ]
      
      tad.1 = tad[,..list1]
      tad.2 = tad[,..list2]
      
      one.run = function(x) return(length(unique(x)) == 1)
      
      same.values.1 = sum(apply(tad.1, 1, one.run))
      same.values.2 = sum(apply(tad.2, 1, one.run))
      
      if(same.values.1 == nrow(tad.1) && same.values.2 == nrow(tad.2)){
        next
      }
      
      statistics1 = t.test(tad.1, tad.2, na.rm = TRUE)
      statistics4 = wilcox.test(unlist(tad.1), unlist(tad.2), na.rm = TRUE, correct = FALSE)
      
      tad_sum[i,5] = statistics1$p.value
      tad_sum[i,6] = statistics4$p.value
    }
    
  }
  
  tad_sum = tad_sum[which(!is.na(tad_sum$ttest)), ]
  
  if (nrow(meta) >= 30){
    
    data.b = tad_sum[,5]
    tad_sum$FDR = p.adjust(unlist(data.b), method = c("BH"), n = nrow(data.b))
    
  } else {
    
    data.b = tad_sum[,6]
    tad_sum$FDR = p.adjust(unlist(data.b), method = c("BH"), n = nrow(data.b))
    
  }
  
  #'
  #' Criterion based on mean of diff across TADs
  #' 
  criterion = as.matrix(summary(tad_sum$mean))
  criterion = criterion[5]
  
  #'
  #' Filtering
  #' 
  tad_sign = tad_sum[which(tad_sum$FDR < FDR_criterion & tad_sum$mean > criterion),]
  
  #'
  #' Merge full TAD table with significant TADs
  #' 
  full.tads = merge(tad_sign, 
                    full, 
                    by.x = "tad_name", 
                    by.y = "tad_name")
  
  genes.found = full.tads$Gene_id
  genes.found = str_split(genes.found, "\\|")
  genes.found = unlist(genes.found)
  genes.found = genes.found[genes.found != "NA"]
  
  keep = colnames(full.tads)
  keep = which(!(keep %in% c("IQR", "ttest", "wilcoxon")))
  full.tads = full.tads[,..keep]
  
  
}

########### Generating output images ########### 
dir.create(image_output_folder, showWarnings = FALSE)

#' Generating image with events in two groups
#' 
#' @param tad_to_visual TAD to visualize

tad_to_visual = c("TAD2130", "TAD854")

# tad_to_visual = c(tad_to_visual,
#                   tad_sign[which(tad_sign$FDR == min(tad_sign$FDR)), ]$tad_name)
# 
# tad_to_visual = c(tad_to_visual,
#                   tad_sign[which(tad_sign$mean == max(tad_sign$mean)), ]$tad_name)


for(i in tad_to_visual){
  tad.test = full %>% dplyr::filter(tad_name == i)
  tad.test.1 = as.matrix(unlist(tad.test[,..list1]))
  tad.test.1 = as.data.frame(tad.test.1)
  tad.test.1$status = paste(group1)
  tad.test.2 = as.matrix(unlist(tad.test[,..list2]))
  tad.test.2 = as.data.frame(tad.test.2)
  tad.test.2$status = paste(group2)
  
  tad.test.plot = rbind(tad.test.1, tad.test.2)
  
  # png(filename = paste(image_output_folder, "/", i, ".png", sep = ""), 
  #     width = 600, height = 820)
  
  # pdf(paste(image_output_folder, "/", i, ".pdf", sep = ""))
  
  gr = ggplot(tad.test.plot, aes(x = status, y = V1, fill = status)) + 
       geom_jitter(position = position_jitter(0.1), shape = 21, color = "gray61", size = 1.5) +
       theme_bw() + 
       theme(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(), 
             legend.position = "none",
             axis.text.x = element_text(size = 15),
             axis.text.y = element_text(size = 15),
             axis.title.y = element_text(size = 16)) +
       labs(y = i, x = "") +
       stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
                    geom = "crossbar", width = 0.2)
  
  saveImageHigh::save_as_pdf({print(gr)},
                             file.name = file.path(image_output_folder,
                                                   paste(i, "allValues.png", sep = "_")),
                             height = 8)
  
  # dev.off()
  
  tad.test = full %>% dplyr::filter(tad_name == i)
  tad.test = tad.test[order(abs(tad.test$diff), decreasing = TRUE), ]

  tad.test.1 = tad.test[,..list1]
  tad.test.2 = tad.test[,..list2]
  tad.test.1$mean = apply(tad.test.1, 1, mean)
  tad.test.2$mean = apply(tad.test.2, 1, mean)
  tad.test.1$status = paste(group1)
  tad.test.2$status = paste(group2)
  tad.test.1$ids = 1:nrow(tad.test.1)
  tad.test.2$ids = 1:nrow(tad.test.2)

  tad.test.1$xj = 1
  tad.test.2$xj = 2

  tad.test.1$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")
  tad.test.2$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")

  tad.test.1[1:min(30, nrow(tad.test.1)), ]$line.color = "gray38"

  tad.test.plot = rbind(tad.test.1[,c("mean", "xj", "ids", "line.color", "status")],
                        tad.test.2[,c("mean", "xj", "ids", "line.color", "status")])



  tad.test.plot$xj = jitter(tad.test.plot$xj, amount = 0.1, factor = 1)
  
  # png(filename = paste(image_output_folder, "/", i, "_topValues.png", sep = ""),
  #     width = 600, height = 820)
  
  # pdf(paste(image_output_folder, "/", i, "_topValues.pdf", sep = ""))
  
  f3 = ggplot(data = tad.test.plot, aes(y = mean, x = xj, fill = status)) +
       geom_line(aes(x = xj, group = ids), color = "lightgray") +
       geom_line(aes(x = xj, group = ids), color = tad.test.plot$line.color) +
       geom_point(aes(x = xj), size = 1.5, shape = 21, color = "gray61") +
    
       # geom_half_boxplot(data=tad.test.plot %>% filter(status == group1), 
       #                   aes(x=xj, y = mean), position = position_nudge(x = -.25),
       #                   side = "r",outlier.shape = NA, center = TRUE, 
       #                   errorbar.draw = FALSE, width = .2) +
       # 
       # geom_half_boxplot(data = tad.test.plot %>% filter(status==group2), 
       #                   aes(x=xj, y = mean), position = position_nudge(x = .15),
       #                   side = "r",outlier.shape = NA, center = TRUE, 
       #                   errorbar.draw = FALSE, width = .2) +
    
       geom_half_violin(data = tad.test.plot %>% filter(status==group1),
                        aes(x = xj, y = mean), position = position_nudge(x = -.3),
                        side = "l") +
    
       geom_half_violin(data = tad.test.plot %>% filter(status==group2),
                        aes(x = xj, y = mean), position = position_nudge(x = .3),
                        side = "r") +
    
       scale_x_continuous(breaks=c(1,2), labels=c("ss6", "ss8"), limits=c(0, 3)) +
       theme_classic() + labs(y = i, x = "") + 
       theme(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(), 
             legend.position = "none",
             axis.text.x = element_text(size = 15),
             axis.text.y = element_text(size = 15),
             axis.title.y = element_text(size = 16))
  
  saveImageHigh::save_as_pdf({print(f3)},
                             file.name = file.path(image_output_folder,
                                                   paste(i, "_topValues.png", sep = "")),
                             height = 8)
                               
  
  # dev.off()
  
  tad.test = full %>% dplyr::filter(tad_name == i)
  tad.test = tad.test[order(tad.test$diff, decreasing = TRUE), ]

  tad.test.1 = tad.test[,..list1]
  tad.test.2 = tad.test[,..list2]
  tad.test.1$mean = apply(tad.test.1, 1, mean)
  tad.test.2$mean = apply(tad.test.2, 1, mean)
  tad.test.1$status = paste(group1)
  tad.test.2$status = paste(group2)
  tad.test.1$ids = 1:nrow(tad.test.1)
  tad.test.2$ids = 1:nrow(tad.test.2)

  tad.test.1$xj = 1
  tad.test.2$xj = 2
  
  tad.test.1 = tad.test.1[order(tad.test$diff, decreasing = TRUE), ]
  tad.test.2 = tad.test.2[order(tad.test$diff, decreasing = TRUE), ]

  tad.test.1$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")
  tad.test.2$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")

  tad.test.1[1:min(30, nrow(tad.test.1)), ]$line.color = "gray38"

  tad.test.plot = rbind(tad.test.1[,c("mean", "xj", "ids", "line.color", "status")],
                        tad.test.2[,c("mean", "xj", "ids", "line.color", "status")])



  tad.test.plot$xj = jitter(tad.test.plot$xj, amount = 0.1, factor = 1)

  # png(filename = paste(image_output_folder, "/", i, "_positive_connections.png", sep = ""),
  #     width = 600, height = 820)
  
  # pdf(paste(image_output_folder, "/", i, "_positive_connections.pdf", sep = ""))

  f3 = ggplot(data=tad.test.plot, aes(y = mean, x = xj, fill = status)) +
       geom_line(aes(x=xj, group=ids), color = "lightgray") +
       geom_line(aes(x=xj, group=ids), color = tad.test.plot$line.color) +
       geom_point(aes(x=xj), size = 1.5, shape = 21, color = "gray61") +
    
       geom_half_violin(data = tad.test.plot %>% filter(status==group1),
                        aes(x = xj, y = mean), position = position_nudge(x = -.3),
                        side = "l") +
    
       geom_half_violin(data = tad.test.plot %>% filter(status==group2),
                        aes(x = xj, y = mean), position = position_nudge(x = .3),
                        side = "r") +
    
       scale_x_continuous(breaks=c(1,2), labels=c("ss6", "ss8"), limits=c(0, 3)) +
       theme_classic() + labs(y = i, x = "") +
       theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             legend.position = "none",
             axis.text.x = element_text(size = 15),
             axis.text.y = element_text(size = 15),
             axis.title.y = element_text(size = 16))
  
  saveImageHigh::save_as_pdf({print(f3)},
                             file.name = file.path(image_output_folder,
                                                   paste(i, "_topPositiveValues.png", sep = "")),
                             height = 8)

  # dev.off()
  
  tad.test = full %>% dplyr::filter(tad_name == i)
  tad.test = tad.test[order(tad.test$diff), ]

  tad.test.1 = tad.test[,..list1]
  tad.test.2 = tad.test[,..list2]
  tad.test.1$mean = apply(tad.test.1, 1, mean)
  tad.test.2$mean = apply(tad.test.2, 1, mean)
  tad.test.1$status = paste(group1)
  tad.test.2$status = paste(group2)
  tad.test.1$ids = 1:nrow(tad.test.1)
  tad.test.2$ids = 1:nrow(tad.test.2)

  tad.test.1$xj = 1
  tad.test.2$xj = 2
  
  tad.test.1 = tad.test.1[order(tad.test$diff), ]
  tad.test.2 = tad.test.2[order(tad.test$diff), ]
  
  tad.test.1$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")
  tad.test.2$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")
  
  tad.test.1[1:min(30, nrow(tad.test.1)), ]$line.color = "gray38"
  
  tad.test.plot = rbind(tad.test.1[,c("mean", "xj", "ids", "line.color", "status")],
                        tad.test.2[,c("mean", "xj", "ids", "line.color", "status")])
  
  
  
  tad.test.plot$xj = jitter(tad.test.plot$xj, amount = 0.1, factor = 1)
  
  # png(filename = paste(image_output_folder, "/", i, "_negative_connections.png", sep = ""),
  #     width = 600, height = 820)
  
  # pdf(paste(image_output_folder, "/", i, "_negative_connections.pdf", sep = ""))
  
  f3 = ggplot(data=tad.test.plot, aes(y = mean, x = xj, fill = status)) +
       geom_line(aes(x=xj, group=ids), color = "lightgray") +
       geom_line(aes(x=xj, group=ids), color = tad.test.plot$line.color) +
       geom_point(aes(x=xj), size = 1.5, shape = 21, color = "gray61") +
          
       geom_half_violin(data = tad.test.plot %>% filter(status==group1),
                        aes(x = xj, y = mean), position = position_nudge(x = -.3),
                        side = "l") +
        
       geom_half_violin(data = tad.test.plot %>% filter(status==group2),
                        aes(x = xj, y = mean), position = position_nudge(x = .3),
                        side = "r") +
        
       scale_x_continuous(breaks=c(1,2), labels=c("ss6", "ss8"), limits=c(0, 3)) +
       theme_classic() + labs(y = i, x = "") + 
       theme(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(), 
             legend.position = "none",
             axis.text.x = element_text(size = 15),
             axis.text.y = element_text(size = 15),
             axis.title.y = element_text(size = 16))
  
  saveImageHigh::save_as_pdf({print(f3)},
                             file.name = file.path(image_output_folder,
                                                   paste(i, "_topNegativeValues.png", sep = "")),
                             height = 8)
  
  # dev.off()
}

########### Clear enviroment ########### 

end_tad_time = Sys.time()

rm(list = setdiff(ls(), c("data.all", "full", "full.tads", "tad_sign", "tad_sum", 
                        "dir_name", "end_tad_time", "start_tad_time",
                        "paired_data", "groups", "FDR_criterion", "genes.found", "output_folder", "x", "res", "res2")))

tad_sign = tad_sign[,c("tad_name", "count", "mean", "FDR")]

########### Generating outputs ###########  
dir.create(output_folder, showWarnings = FALSE)

write.table(full, paste(output_folder, "/integrated_table_with_tads.csv", sep = ""),
            row.names = FALSE, sep = "\t", quote = FALSE)

if(!is.null(groups)){

  write.table(tad_sum, paste(output_folder, "/tad_statistics.csv", sep = ""),
              row.names = FALSE, sep = "\t", quote = FALSE)

  write.table(full.tads, paste(output_folder, "/integrated_table_with_sign_tads.csv", sep = ""),
              row.names = FALSE, sep = "\t", quote = FALSE)

  write.table(tad_sign, paste(output_folder, "/sign_tad_statistics.csv", sep = ""),
              row.names = FALSE, sep = "\t", quote = FALSE)

  write.table(genes.found, paste(output_folder, "/genes_found.txt", sep = ""),
              row.names = FALSE, col.names = FALSE, quote = FALSE)


}
