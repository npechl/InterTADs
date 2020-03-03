############################# Getting required libraries ###########################

library(stringr)
library(karyoploteR)
library(data.table)
library(GenomicRanges)
library(png)

############################### Inputs ##############################

dir_name = "Data Integration"

full.tads = fread(paste(dir_name, "/integrated_table_with_tads.csv", sep = ""))
tad.sign = fread(paste(dir_name, "/sign_tad_statistics.csv", sep = ""))
meta = fread(paste(dir_name, "/meta-data.csv", sep = ""))

image_folder_name = "output_visualizations"
dir.create(image_folder_name, showWarnings = FALSE)
image_folder_name = paste(image_folder_name, "/", sep = "")

groups = unique(meta$groups)

by = "FDR"
tad_to_visual = tad.sign[which(tad.sign$FDR == max(tad.sign$FDR)), ]$tad_name

# by = "diff"
# t = abs(full.tads$diff)
# tad_to_visual = full.tads[which(t == max(t)), ]$tad_name

full.vtads = full.tads[which(full.tads$tad_name %in% tad_to_visual), ]

full.vtads$chromosome_name = paste("chr", full.vtads$chromosome_name,  sep = "")

chr_to_visual = as.character(unique(full.vtads$chr))
start = min(unique(full.vtads$start))
end = max(unique(full.vtads$end))
region = toGRanges(data.frame(chr_to_visual, start, end))
tick.distance = (end - start)/5
# tick.distance = tick.distance + (tick.distance / 5.9)

####################### Reading data ##########################

columns = c("ID", "chromosome_name", "start_position", "end_position", meta$newNames, "diff")

range = full.vtads$end_position - full.vtads$start_position

rain_data = full.vtads[which(range == 1), ..columns]

box_data = full.vtads[which(range > 1), ..columns]
box_data = box_data %>% group_by(start_position, end_position) %>% filter(diff == max(diff))
box_data = as.data.table(box_data)

############################## Plotting ##############################

keep = meta$newNames[which(meta$groups == groups[1])]

rain_data$mean = apply(rain_data[,..keep], 1, mean)
box_data$mean = apply(box_data[,..keep], 1, mean)

png(filename = paste(image_folder_name, groups[1], "_", by, ".png", sep = ""), width = 900, height = 900)

kp = plotKaryotype(genome = "hg19", plot.type = 4, zoom = region, cex = 1.5)
kpDataBackground(kp, data.panel = 1)
kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = tick.distance, cex = 1.25)

if(nrow(rain_data) > 0){
  rain_color = ifelse(rain_data$mean <= 30, "green4", ifelse(rain_data$mean <= 70, "blue", "red"))
  numeric.vector = as.numeric(unlist(rain_data$mean))
  kpPoints(kp, chr = rain_data$chr,
           x = rain_data$end_position, y = numeric.vector / max(numeric.vector),
           col = rain_color, cex = 0.8)
}

if(nrow(box_data) > 0){
  kpAxis(kp, r0 = 1.55, r1 = 2, ymin = 0, ymax = max(box_data$mean))
  numeric.vector = as.numeric(unlist(box_data$mean))
  
  kpBars(kp, chr = box_data$chr,
         x0 = box_data$start_position, x1 = box_data$end_position,
         y1 = (numeric.vector / max(numeric.vector)), border = "slateblue4")
}

kpAddLabels(kp, labels = paste(tad_to_visual, collapse = ","), side = "left", srt = 90, label.margin = 0.03, cex = 1.5)

dev.off()

keep = meta$newNames[which(meta$groups == groups[2])]

rain_data$mean = apply(rain_data[,..keep], 1, mean)
box_data$mean = apply(box_data[,..keep], 1, mean)

png(filename = paste(image_folder_name, groups[2], "_", by, ".png", sep = ""), width = 900, height = 900)

kp = plotKaryotype(genome = "hg19", plot.type = 4, zoom = region, cex = 1.5)
kpDataBackground(kp, data.panel = 1)
kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = tick.distance, cex = 1.25)

if(nrow(rain_data) > 0){
  rain_color = ifelse(rain_data$mean <= 30, "green4", ifelse(rain_data$mean <= 70, "blue", "red"))
  numeric.vector = as.numeric(unlist(rain_data$mean))
  kpPoints(kp, chr = rain_data$chr,
           x = rain_data$end_position, y = numeric.vector / max(numeric.vector),
           col = rain_color, cex = 0.8)
}

if(nrow(box_data) > 0){
  kpAxis(kp, r0 = 1.55, r1 = 2, ymin = 0, ymax = max(box_data$mean))
  numeric.vector = as.numeric(unlist(box_data$mean))
  
  kpBars(kp, chr = box_data$chr,
         x0 = box_data$start_position, x1 = box_data$end_position,
         y1 = (numeric.vector / max(numeric.vector)), border = "slateblue4")
}

kpAddLabels(kp, labels = paste(tad_to_visual, collapse = ","), side = "left", srt = 90, label.margin = 0.03, cex = 1.5)

dev.off()