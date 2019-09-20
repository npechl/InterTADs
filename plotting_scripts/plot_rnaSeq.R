############################# Getting required libraries ###########################

library(stringr)
library(stringi)
library(karyoploteR)
library(data.table)
library(png)
library(magick)

############################### Inputs ##############################

# file to read
file_name = "biodata_pre_integration_14.csv"
pre_data = read.csv(file_name, sep = ";")

image_folder_name = "plotting_vol1"
image_folder_name = paste(image_folder_name, "/", sep = "")

# chromosomes to visualize
chr_to_visual = "chr14"

# patient to visualize
patient = colnames(pre_data)[8:length(pre_data)]
# patient = which(colnames(pre_data) %in% patient)

# scale that has been applied
scale = NULL

####################### Reading data ##########################

gene_data = pre_data[which(str_detect(pre_data$ID, "ENSG")), ]
gene_data$chromosome_name = paste("chr", gene_data$chromosome_name, sep = "")

# gene_data[,8:ncol(gene_data)] = log10(gene_data[,8:ncol(gene_data)] + 1)

############################## Plotting ##############################

for(i in patient){
  png(filename = paste(image_folder_name, i, ".png", sep = ""), width = 2494, height = 1812)
  
  kp = plotKaryotype(genome="hg19", chromosomes = chr_to_visual, plot.type = 4)
  kpAddBaseNumbers(kp, add.units = TRUE)
  
  kpDataBackground(kp, data.panel = 1)
  kpAxis(kp, ymin = 0, ymax = (max(gene_data[ ,i])))
  kpBars(kp, chr = gene_data$chromosome_name, 
         x0 = gene_data$start_position, x1 = gene_data$end_position, 
         y1 = (gene_data[ ,i] / max(gene_data[ ,i])), border = "slateblue4")
  kpAddLabels(kp, labels = "RNASeq", side = "left", srt = 90, label.margin = 0.03)
  
  dev.off()
}

rm(i)

gene_data[,8:ncol(gene_data)] = log10(gene_data[,8:ncol(gene_data)] + 1)

for(i in seq(1, length(patient), by = 2)){
  img_pre = image_read(paste(image_folder_name, patient[i], ".png", sep = ""))
  img_post = image_read(paste(image_folder_name, patient[i + 1], ".png", sep = ""))
  
  img = c(img_pre, img_post)
  img = image_append(img)
  
  image_write(img, path = paste(image_folder_name, patient[i], "_", patient[i + 1], ".png", sep = ""), format = "png")
  
  png(filename = paste(image_folder_name, patient[i], "_", patient[i + 1], "_diff.png", sep = ""), width = 2494, height = 1812)
  
  kp = plotKaryotype(genome="hg19", chromosomes = chr_to_visual, plot.type = 4)
  kpAddBaseNumbers(kp, add.units = TRUE)
  
  data_to_visualize = gene_data[ ,patient[i]] - gene_data[ ,patient[i + 1]]
  col_to_apply = ifelse(data_to_visualize > 0, "red4", "blue4")
  
  kpDataBackground(kp, data.panel = 1)
  kpAxis(kp, ymin = (min(data_to_visualize)), ymax = (max(data_to_visualize)))
  
  Y0 = ifelse(data_to_visualize > 0, 0, data_to_visualize) # / abs(min(data_to_visualize))
  Y1 = ifelse(data_to_visualize > 0, data_to_visualize, 0) # / (max(data_to_visualize))
  
  Y0 = Y0 * 0.5 / abs(min(Y0))
  Y1 = Y1 * 0.5 / max(Y1)
  
  kpRect(kp, chr = gene_data$chromosome_name,
         x0 = gene_data$start_position, x1 = gene_data$end_position,
         y0 = (Y0 + 0.5),
         y1 = (Y1 + 0.5), border = col_to_apply)
  kpAddLabels(kp, labels = "RNASeq", side = "left", srt = 90, label.margin = 0.03)
  
  dev.off()
}

rm(gene_data, pre_data, chr_to_visual, file_name, 
   i, patient, img_pre, img_post, img, kp, image_folder_name, data_to_visualize)
rm(col_to_apply, Y0, Y1, scale)