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

# by = "FDR"
# tad_to_visual = tad.sign[which(tad.sign$FDR == max(tad.sign$FDR)), ]$tad_name

by = "diff"
t = abs(full.tads$diff)
tad_to_visual = full.tads[which(t == max(t)), ]$tad_name

full.vtads = full.tads[which(full.tads$tad_name %in% tad_to_visual), ]

full.vtads$chromosome_name = paste("chr", full.vtads$chromosome_name,  sep = "")

chr_to_visual = as.character(unique(full.vtads$chr))
start = min(unique(full.vtads$start))
end = max(unique(full.vtads$end))
region = toGRanges(data.frame(chr_to_visual, start, end))
tick.distance = (end - start)/5
tick.distance = tick.distance + (tick.distance / 5.9)

####################### Reading data ##########################

columns = c("ID", "chromosome_name", "start_position", "end_position", meta$newNames, "diff")

range = full.vtads$end_position - full.vtads$start_position

rain_data = full.vtads[which(range == 1), ..columns]

box_data = full.vtads[which(range > 1), ..columns]
box_data = box_data %>% group_by(start_position, end_position) %>% filter(abs(diff) == max(abs(diff)))
box_data = as.data.table(box_data)

############################## Plotting ##############################

for(i in 1:length(groups)){
  keep = meta$newNames[which(meta$groups == groups[i])]
  
  rain_data$mean = apply(rain_data[,..keep], 1, mean)
  box_data$mean = apply(box_data[,..keep], 1, mean)
  
  maximum_value = 100
  minimum_value = 0
  
  png(filename = paste(image_folder_name, groups[i], "_max_", by, ".png", sep = ""), width = 1150, height = 900)
  
  kp = plotKaryotype(genome = "hg19", plot.type = 4, zoom = region, cex = 1.5)
  kpDataBackground(kp, data.panel = 1)
  kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = tick.distance, cex = 1.25)
  
  lb = c(minimum_value, (maximum_value + minimum_value)/2, maximum_value)
  lb = round(lb, digits = 2)
  lb = paste(lb, "%", sep = " ")
  kpAxis(kp, r0 = 0, r1 = 1, ymin = minimum_value, ymax = maximum_value, side = 2, labels = lb)
  
  if(nrow(rain_data) > 0){
    rain_data = rain_data[which(rain_data$mean <= maximum_value & rain_data$mean >= minimum_value), ]
    rain_color = ifelse(rain_data$mean <= 30, "green4", ifelse(rain_data$mean <= 70, "blue", "red"))
    numeric.vector = rain_data$mean
    numeric.vector = numeric.vector / maximum_value
    kpPoints(kp, chr = rain_data$chr,
             x = rain_data$end_position, y = numeric.vector,
             col = rain_color, cex = 0.8)
  }
  
  if(nrow(box_data) > 0){
    box_data = box_data[which(box_data$mean <= maximum_value & box_data$mean >= minimum_value), ]
    numeric.vector = box_data$mean
    numeric.vector = numeric.vector / maximum_value
    kpBars(kp, chr = box_data$chr,
           x0 = box_data$start_position, x1 = box_data$end_position,
           y1 = numeric.vector, border = "slateblue4")
  }
  
  kpAddLabels(kp, labels = paste(tad_to_visual, collapse = ","), side = "left", srt = 90, label.margin = 0.03, cex = 1.5)
  
  dev.off()
}

png(filename = paste(image_folder_name, "Group_differences_max_", by, "_zooming.png", sep = ""), width = 1150, height = 900)

kp = plotKaryotype(genome = "hg19", plot.type = 4, zoom = region, cex = 1.5)
kpDataBackground(kp, data.panel = 1)
kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = tick.distance, cex = 1.25)

all = c(rain_data$diff, box_data$diff)

maximum_value = max(all)
minimum_value = min(all)

lb = c(minimum_value, 0, maximum_value)
lb = round(lb, digits = 2)
lb = paste(lb, "%", sep = " ")
kpAxis(kp, r0 = 0, r1 = 1, ymin = minimum_value, ymax = maximum_value, side = 2, labels = lb)

if(nrow(rain_data) > 0){
  numeric.vector = rain_data$diff
  who = numeric.vector > 0
  numeric.vector[who] = numeric.vector[who] / maximum_value
  numeric.vector[!who] = numeric.vector[!who] / abs(minimum_value)
  numeric.vector = (numeric.vector + 1)/2

  kpPoints(kp, chr = rain_data$chr,
           x = rain_data$end_position, y = numeric.vector, cex = 0.8)
}

if(nrow(box_data) > 0){
  lb = c(min(box_data$diff), 0, max(box_data$diff))
  lb = round(lb, digits = 2)
  lb = paste(lb, "%", sep = " ")
  numeric.vector = box_data$diff
  who = numeric.vector > 0
  numeric.vector[who] = numeric.vector[who] / maximum_value
  numeric.vector[!who] = numeric.vector[!who] / abs(minimum_value)
  numeric.vector = (numeric.vector + 1)/2

  bottom = rep(0.5, length(numeric.vector))
  bottom[!who] = numeric.vector[!who]

  up = rep(0.5, length(numeric.vector))
  up[who] = numeric.vector[who]

  kpBars(kp, chr = box_data$chr,
         x0 = box_data$start_position, x1 = box_data$end_position,
         y0 = bottom, y1 = up,
         border = "slateblue4")
}

kpAddLabels(kp, labels = paste(tad_to_visual, collapse = ","), side = "left", srt = 90, label.margin = 0.03, cex = 1.5)

dev.off()

rm("t", "lb", "i", "columns", "rain_color", "numeric.vector")