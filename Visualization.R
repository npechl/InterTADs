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
dir.create(image_folder_name)
image_folder_name = paste(image_folder_name, "/", sep = "")

by = "FDR"
tad_to_visual = tad.sign[which(tad.sign$FDR == max(tad.sign$FDR)), ]$tad_name

full.vtads = full.tads[which(full.tads$tad_name %in% tad_to_visual), ]

full.vtads$chromosome_name = paste("chr", full.vtads$chromosome_name,  sep = "")

chr_to_visual = as.character(unique(full.vtads$chr))
start = min(unique(full.vtads$start))
end = max(unique(full.vtads$end))
region = toGRanges(data.frame(chr_to_visual, start, end))

####################### Reading data ##########################

columns = c("ID", "chromosome_name", "start_position", "end_position", meta$newNames)

range = full.vtads$end_position - full.vtads$start_position

rain_data = full.vtads[which(range == 1), ..columns]

box_data = full.vtads[which(range > 1), ..columns]

############################## Plotting ##############################

for(i in columns[5:length(columns)]){
  png(filename = paste(image_folder_name, i, "_", by, ".png", sep = ""))

  kp = plotKaryotype(genome = "hg19", chromosomes = chr_to_visual, plot.type = 4, zoom = region)
  kpAddBaseNumbers(kp, add.units = TRUE)
  
  kpDataBackground(kp, data.panel = 1)
  
  if(nrow(rain_data) > 0){
    rain_color = ifelse(rain_data[ ,..i] <= 30, "green4", ifelse(rain_data[ ,..i] <= 70, "blue", "red"))
    numeric.vector = as.numeric(unlist(rain_data[,..i]))
    kpPoints(kp, chr = rain_data$chr, 
             x = rain_data$end_position, y = numeric.vector / max(numeric.vector), 
             col = rain_color)
  }
  
  if(nrow(box_data) > 0){
    kpAxis(kp, r0 = 1.55, r1 = 2, ymin = 0, ymax = max(box_data[ ,..i]))
    numeric.vector = as.numeric(unlist(box_data[,..i]))
    
    kpBars(kp, chr = box_data$chr, 
           x0 = box_data$start_position, x1 = box_data$end_position, 
           y1 = (numeric.vector / max(numeric.vector)), border = "slateblue4")
  }
  
  kpAddLabels(kp, labels = paste(tad_to_visual, collapse = ","), side = "left", srt = 90, label.margin = 0.03)
  
  dev.off()
}

