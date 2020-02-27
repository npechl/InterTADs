########## Loading libraries ########## 
library(tidyverse)
library(GenomicRanges)
library(ggplot2)
library(gplots)

start_tad_time = Sys.time()

########### Inputs ########## 

dir_name = "Data Integration"

# data.all = read.csv(paste(dir_name, "/integrated_table.csv", sep = ""))
data.all = fread(paste(dir_name, "/integrated_table.csv", sep = ""))

# TAD = read.csv(paste(dir_name, "/hglft_genome_2dab_ec1330.bed", sep = ""), 
#                header = F, sep = "\t", check.names = FALSE)

TAD = fread(paste(dir_name, "/hglft_genome_2dab_ec1330.bed", sep = ""), 
            header = F, sep = "\t", check.names = FALSE)

# meta = read.csv(paste(dir_name, "/meta-data.csv", sep = ""))
meta = fread(paste(dir_name, "/meta-data.csv", sep = ""))

#4. set up the groups of interest
# group1 = "ss6"
# group2 = "ss8"

groups = meta$groups
groups = as.character(groups)
groups = unique(groups)
groups = groups[which(groups != "")]

group1 = groups[1]
group2 = groups[2]

paired.data = FALSE

FDR_criterion = 0.05

########### Processing ########### 

# generate a barcode 
data.all$name <- paste0(data.all$ID, data.all$VarAnnotation, data.all$Gene_id)

# reordering
data.over <- data.all[,c("name", "chromosome_name", "start_position", "end_position")]
data.over$chromosome_name = paste("chr", data.over$chromosome_name, sep = "")
data.over = data.over[,c(2,3,4,1)]

# overlap of the TAD with events
colnames(TAD) <- paste(c("chr", "start", "end", "name"))
colnames(data.over) <- paste(c("chr", "start", "end", "name"))

# Make GRanges object
gr1 = with(TAD, GRanges(chr, IRanges(start = start, end = end, names = name)))
gr2 = with(data.over, GRanges(chr, IRanges(start = start, end = end, names = name)))

# Completely overlapping
type1 = findOverlaps(query = gr1, subject = gr2)
# type1.df = data.frame(TAD[queryHits(type1),], data.over[subjectHits(type1),])

type1.df = cbind(TAD[queryHits(type1),], data.over[subjectHits(type1),])
type1.df = type1.df[,c(2,3,4,8)]
colnames(type1.df) = c("tad_start", "tad_end", "tad_name", "name.1")
# test <- type1.df[duplicated(type1.df$name.1), ]

# merge the values
full <- merge(type1.df, data.all, by.x = "name.1", by.y = "name")
keep = c("chromosome_name", "tad_name", "tad_start", "tad_end", "ID", "start_position", "end_position", "Gene_id", "Gene_locus", meta$newNames)
full = full[,..keep]

# test for NAs, remove them for the statistical analysis
full[is.na(full)] <- 0

# average values between groups
# annotation.1 <- annotation[which(annotation$SUBSET == group1),]
# annotation.2 <- annotation[which(annotation$SUBSET == group2),]

list1 = meta[which(meta$groups == group1), ]$newNames
list2 = meta[which(meta$groups == group2), ]$newNames

# list1 <- as.character(annotation.1$RNAseq)
full.1 <- full[,..list1]

# list2 <- as.character(annotation.2$RNAseq)
full.2 <- full[,..list2]

full$diff <- apply(full.2, 1, mean, na.rm = TRUE) - apply(full.1, 1, mean, na.rm = TRUE)

# exclude the events with no differences
full <- full[which(round(abs(full$diff))>0),]

###########  TADS ########### 

# Statistics
tad_sum <- group_by(full, tad_name) %>% summarise(count = n(), 
                                              mean = mean(abs(diff), na.rm = TRUE),
                                              IQR = IQR(diff, na.rm = TRUE))

tad_sum$ttest <- numeric(length = nrow(tad_sum))
tad_sum$wilcoxon <- numeric(length = nrow(tad_sum))
# tad_sum = as.data.frame(tad_sum)
tad_sum = as.data.table(tad_sum)

# parametric t.test
# non-parametric wilcoxon

if(paired.data){
  for (i in 1:nrow(tad_sum)){
    tad <- full[which(full$name == as.character(tad_sum[i,1])), ]
    
    tad.1 <- tad[,..list1]
    tad.2 <- tad[,..list2]
    
    statistics2 <- t.test(unlist(tad.1), unlist(tad.2), na.rm = TRUE, paired = T)
    statistics3 <- wilcox.test(unlist(tad.1), unlist(tad.2), na.rm = TRUE, paired = T)
    
    tad_sum[i,5] <- statistics2$p.value
    tad_sum[i,6] <- statistics3$p.value
  }
} else {
  for (i in 1:nrow(tad_sum)){
    tad <- full[which(full$tad_name == as.character(tad_sum[i,1])), ]
    
    tad.1 <- tad[,..list1]
    tad.2 <- tad[,..list2]
    
    statistics1 <- t.test(tad.1, tad.2, na.rm = TRUE)
    statistics4 <- wilcox.test(unlist(tad.1), unlist(tad.2), na.rm = TRUE, correct = FALSE)
    
    tad_sum[i,5] <- statistics1$p.value
    tad_sum[i,6] <- statistics4$p.value
  }
}

# test for FDR - based on paired wilc test [,7] because study group <30 cases

if (nrow(meta) >= 30){
  data.b = tad_sum[,5]
  tad_sum$FDR <- p.adjust(unlist(data.b), method = c("BH"), n = nrow(data.b))
} else {
  data.b = tad_sum[,6]
  tad_sum$FDR <- p.adjust(unlist(data.b), method = c("BH"), n = nrow(data.b))
}

# criterion based on mean of diff across TADs
criterion = as.matrix(summary(tad_sum$mean))
criterion = criterion[5]

# filtering
tad_sign <- tad_sum[which(tad_sum$FDR < FDR_criterion & tad_sum$mean > criterion),]

# final file - export
full.tads <- merge(tad_sign, full, by.x = "tad_name", by.y = "tad_name")

# # plot for 1 example
# tad.test <- full %>% filter(name == "TAD1361")
# # colnames(tad.test)
# # row.names(heat.data)<- tad.test$name.1
# 
# tad.test.1 <- as.matrix(unlist(tad.test[,list1]))
# tad.test.1 <- as.data.frame(tad.test.1)
# tad.test.1$status <- paste(group1)
# tad.test.2 <- as.matrix(unlist(tad.test[,list2]))
# tad.test.2 <- as.data.frame(tad.test.2)
# tad.test.2$status <- paste(group2)
# 
# tad.test.plot <- rbind(tad.test.1, tad.test.2)
# 
# ggplot(tad.test.plot, aes(x=status, y=V1, fill=status)) + 
#   geom_jitter(position=position_jitter(0.1), shape=21, color="gray61", size=1.5) +
#   #scale_fill_manual(values=c("orange", "steelblue2","chartreuse", "yellow","red", "black", "grey")) +
#   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5)
# 
# ###heatmap
# 
# tad.test <- full %>% filter(name == "TAD2060")
# # colnames(tad.test)
# heat.data <- tad.test[,as.character(annotation$RNAseq)]
# row.names(heat.data) <- tad.test$ID
# 
# heatmap.2(as.matrix(heat.data),  
#           dendrogram = "both",  
#           notecol = "black", 
#           density = "density", 
#           trace = "none", 
#           margins = c(8, 13), 
#           col = greenred(75))

end_tad_time = Sys.time()

rm(list=setdiff(ls(), c("data.all", "full", "full.tads", "tad_sign", "tad_sum", 
                        "dir_name", "end_tad_time", "start_tad_time",
                        "paired_data", "groups", "FDR_criterion")))

########### Generating outputs ###########  

write.table(full, paste(dir_name, "/integrated_table_with_tads.csv", sep = ""), row.names = FALSE, sep = "\t")
write.table(tad_sum, paste(dir_name, "/tad_statistics.csv", sep = ""), row.names = FALSE, sep = "\t")

write.table(full.tads, paste(dir_name, "/integrated_table_with_sign_tads.csv", sep = ""), row.names = FALSE, sep = "\t")
write.table(tad_sign, paste(dir_name, "/sign_tad_statistics.csv", sep = ""), row.names = FALSE, sep = "\t")

