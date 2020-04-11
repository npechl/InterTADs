############################# Getting required libraries ###########################

library(stringr)
library(karyoploteR)
library(data.table)
library(GenomicRanges)
library(png)
library(dplyr)

############################### Inputs ##############################

dir_name = "Datasets"
tad_folder = "output_tables"

full.tads = fread(paste(tad_folder, "/integrated_table_with_tads.csv", sep = ""))
tad.sign = fread(paste(tad_folder, "/sign_tad_statistics.csv", sep = ""))
meta = fread(paste(dir_name, "/meta-data.csv", sep = ""))
who = meta == ""
who = apply(who, 1, sum)
meta = meta[which(who == 0), ]

groups = unique(meta$groups)

image_folder_name = "output_visualizations"
dir.create(image_folder_name, showWarnings = FALSE)
image_folder_name = paste(image_folder_name, "/", sep = "")

tad_to_visual = c()

tad_to_visual = c(tad_to_visual,
                  tad.sign[which(tad.sign[,"mean"] == max(tad.sign[,"mean"])), ]$tad_name)

tad_to_visual = c(tad_to_visual,
                  tad.sign[which(tad.sign[,"FDR"] == min(tad.sign[,"FDR"])), ]$tad_name)

############################### For every sign TAD ##############################

for(j in tad_to_visual){
  full.vtads = full.tads[which(full.tads$tad_name %in% j), ]
  full.vtads$chromosome_name = paste("chr", full.vtads$chromosome_name,  sep = "")
  
  chr_to_visual = as.character(unique(full.vtads$chr))
  start = min(unique(full.vtads$start))
  end = max(unique(full.vtads$end))
  region = toGRanges(data.frame(chr_to_visual, start, end))
  tick.distance = (end - start)/6
  # tick.distance = tick.distance + (tick.distance / 5.9)
  
  ####################### Reading data ##########################
  
  columns = c("ID", "chromosome_name", "start_position", "end_position", meta$newNames, "diff")
  
  range = full.vtads$end_position - full.vtads$start_position
  
  rain_data = full.vtads[which(range == 1 & !str_detect(full.vtads$ID, ":")), ..columns]
  
  box_data = full.vtads[which(range > 1), ..columns]
  
  mut_data = full.vtads[which(range == 1 & str_detect(full.vtads$ID, ":")), ..columns]
  
  ############################## Plotting ##############################
  dir.create(paste(image_folder_name, "patients/", sep = ""), showWarnings = FALSE)
  
  for(i in meta$newNames){
    png(filename = paste(image_folder_name, "patients/", i, "_", j, ".png", sep = ""), width = 1150, height = 900)
    
    kp = plotKaryotype(genome = "hg19", plot.type = 4, zoom = region, cex = 1.5)
    kpDataBackground(kp, data.panel = 1)
    kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = tick.distance, cex = 1.25)
    lb = c("0 %", "50 %", "100 %")
    kpAxis(kp, r0 = 0, r1 = 1, ymin = 0, ymax = 100, side = 2, labels = lb)
    
    if(nrow(rain_data) > 0){
      rain_color = ifelse(rain_data[ ,..i] <= 30, "green4", ifelse(rain_data[ ,..i] <= 70, "blue", "red"))
      numeric.vector = as.numeric(unlist(rain_data[,..i]))
      numeric.vector = numeric.vector / 100
      kpPoints(kp, chr = rain_data$chromosome_name,
               x = rain_data$end_position, y = numeric.vector, cex = 0.8)
    }
    
    if(nrow(box_data) > 0){
      numeric.vector = as.numeric(unlist(box_data[,..i]))
      numeric.vector = numeric.vector / 100
      kpBars(kp, chr = box_data$chr,
             x0 = box_data$start_position, x1 = box_data$end_position,
             y1 = numeric.vector, border = "slateblue4", cex = 1.2)
    }
    
    if(nrow(mut_data) > 0){
      mut_data = mut_data[which(mut_data[,..i] > 0), ]
      numeric.vector = as.numeric(unlist(mut_data[,..i]))
      numeric.vector = numeric.vector / 100
      kpPoints(kp, chr = rain_data$chromosome_name,
               x = mut_data$end_position, y = numeric.vector,
               col = "red", cex = 1.2)
    }
    
    kpAddLabels(kp, labels = j, side = "left", srt = 90, label.margin = 0.03, cex = 1.5)
    
    dev.off()
  }
  
  dir.create(paste(image_folder_name, j, sep = ""), showWarnings = FALSE)
  
  ####################### Filtering data ##########################
  
  columns = c("ID", "chromosome_name", "start_position", "end_position", meta$newNames, "diff")
  
  range = full.vtads$end_position - full.vtads$start_position
  
  rain_data = full.vtads[which(range == 1 & !str_detect(full.vtads$ID, ":")), ..columns]
  
  box_data = full.vtads[which(range > 1), ..columns]
  box_data = box_data %>% group_by(start_position, end_position) %>% filter(abs(diff) == max(abs(diff)))
  box_data = as.data.table(box_data)
  
  mut_data = full.vtads[which(range == 1 & str_detect(full.vtads$ID, ":")), ..columns]
  
  ############################## Plotting ##############################
  
  for(i in 1:length(groups)){
    keep = meta$newNames[which(meta$groups == groups[i])]
    
    rain_data$mean = apply(rain_data[,..keep], 1, mean)
    box_data$mean = apply(box_data[,..keep], 1, mean)
    mut_data$mean = apply(mut_data[,..keep], 1, mean)
    
    maximum_value = 100
    minimum_value = 0
    
    png(filename = paste(image_folder_name, j, "/", groups[i], ".png", sep = ""), 
        width = 1150, height = 900)
    
    kp = plotKaryotype(genome = "hg19", plot.type = 4, zoom = region, cex = 1.5)
    kpDataBackground(kp, data.panel = 1)
    kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = tick.distance, cex = 1.25)
    
    lb = c(minimum_value, (maximum_value + minimum_value)/2, maximum_value)
    lb = round(lb, digits = 2)
    lb = paste(lb, "%", sep = " ")
    kpAxis(kp, r0 = 0, r1 = 1, ymin = minimum_value, ymax = maximum_value, side = 2, labels = lb)
    
    rain_data = rain_data[which(rain_data$mean <= maximum_value & rain_data$mean >= minimum_value), ]
    if(nrow(rain_data) > 0){
      rain_color = ifelse(rain_data$mean <= 30, "green4", ifelse(rain_data$mean <= 70, "blue", "red"))
      numeric.vector = rain_data$mean
      numeric.vector = numeric.vector / maximum_value
      kpPoints(kp, chr = rain_data$chromosome_name,
               x = rain_data$end_position, y = numeric.vector, cex = 0.8)
    }
    
    box_data = box_data[which(box_data$mean <= maximum_value & box_data$mean >= minimum_value), ]
    if(nrow(box_data) > 0){
      numeric.vector = box_data$mean
      numeric.vector = numeric.vector / maximum_value
      kpBars(kp, chr = box_data$chr,
             x0 = box_data$start_position, x1 = box_data$end_position,
             y1 = numeric.vector, border = "slateblue4", cex = 1.2)
    }
    
    mut_data = mut_data[which(mut_data$mean <= maximum_value & mut_data$mean >= minimum_value & mut_data$mean > 0), ]
    
    if(nrow(mut_data) > 0){
      numeric.vector = mut_data$mean
      numeric.vector = numeric.vector / 100
      kpPoints(kp, chr = mut_data$chromosome_name,
               x = mut_data$end_position, y = numeric.vector,
               col = "red", cex = 1.2)
    }
    
    kpAddLabels(kp, labels = j, side = "left", srt = 90, label.margin = 0.03, cex = 1.5)
    
    dev.off()
  }
  
  png(filename = paste(image_folder_name, j, "/", "Group_differences", ".png", sep = ""), 
      width = 1150, height = 900)
  
  kp = plotKaryotype(genome = "hg19", plot.type = 4, zoom = region, cex = 1.5)
  kpDataBackground(kp, data.panel = 1)
  kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = tick.distance, cex = 1.25)
  
  all = c(rain_data$diff, box_data$diff)
  
  maximum_value = 100
  minimum_value = -100
  
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
    
    kpPoints(kp, chr = rain_data$chromosome_name,
             x = rain_data$end_position, y = numeric.vector, cex = 0.8)
  }
  
  if(nrow(box_data) > 0){
    numeric.vector = box_data$diff
    who = numeric.vector > 0
    numeric.vector[who] = numeric.vector[who] / maximum_value
    numeric.vector[!who] = numeric.vector[!who] / abs(minimum_value)
    numeric.vector = (numeric.vector + 1)/2
    
    bottom = rep(0.5, length(numeric.vector))
    bottom[!who] = numeric.vector[!who]
    
    up = rep(0.5, length(numeric.vector))
    up[who] = numeric.vector[who]
    
    kpBars(kp, chr = box_data$chromosome_name,
           x0 = box_data$start_position, x1 = box_data$end_position,
           y0 = bottom, y1 = up,
           border = "slateblue4", cex = 1.2)
  }
  
  if(nrow(mut_data) > 0){
    numeric.vector = mut_data$diff
    who = numeric.vector > 0
    numeric.vector[who] = numeric.vector[who] / maximum_value
    numeric.vector[!who] = numeric.vector[!who] / abs(minimum_value)
    numeric.vector = (numeric.vector + 1)/2
    
    kpPoints(kp, chr = mut_data$chromosome_name,
             x = mut_data$end_position, y = numeric.vector, 
             col = "red", cex = 1.2)
  }
  
  kpAddLabels(kp, labels = j, side = "left", srt = 90, label.margin = 0.03, cex = 1.5)
  
  dev.off()
}


rm(list = ls())
