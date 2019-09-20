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

meth_data = pre_data[which(str_detect(pre_data$ID, "cg")), ]
meth_data$chromosome_name = paste("chr", meth_data$chromosome_name, sep = "")

if(!is.null(scale)){
  meth_data[,8:ncol(meth_data)] = meth_data[,8:ncol(meth_data)] / scale
}

############################## Plotting ##############################

# plot.params = getDefaultPlotParams(plot.type = 1)
# plot.params$data1height = 138
# plot.params$bottommargin = 10
# plot.params$data1inmargin = 5
# plot.params$rightmargin = 0.02
# plot.params$leftmargin = 0.04
# plot.params$ideogramheight = 10

for(i in patient){
  png(filename = paste(image_folder_name, i, ".png", sep = ""), width = 2494, height = 1812)
  
  kp = plotKaryotype(genome="hg19", chromosomes = chr_to_visual, plot.type = 4)
  kpAddBaseNumbers(kp, add.units = TRUE)
  
  kpDataBackground(kp, data.panel = 1)
  # kpAxis(kp, r0 = 0, r1 = 0.45, ymin = 0, ymax = 100)
  meth_color = character(nrow(meth_data))
  meth_color = ifelse(meth_data[ ,i] <= 0.3, "green4", ifelse(meth_data[ ,i] <= 0.7, "blue", "red"))
  kpPoints(kp, chr = meth_data$chromosome_name, x = meth_data$end_position, y = meth_data[ ,i], col = meth_color, cex = 0.7)
  kpAddLabels(kp, labels = "DNAMethylation", side = "left", srt = 90, label.margin = 0.03)
  
  dev.off()
}

rm(i)

for(i in seq(1, length(patient), by = 2)){
  img_pre = image_read(paste(image_folder_name, patient[i], ".png", sep = ""))
  img_post = image_read(paste(image_folder_name, patient[i + 1], ".png", sep = ""))
  
  img = c(img_pre, img_post)
  img = image_append(img)
  
  image_write(img, path = paste(image_folder_name, patient[i], "_", patient[i + 1], ".png", sep = ""), format = "png")
  
  png(filename = paste(image_folder_name, patient[i], "_", patient[i + 1], "_diff.png", sep = ""), width = 2494, height = 1812)
  
  kp = plotKaryotype(genome="hg19", chromosomes = chr_to_visual, plot.type = 4)
  kpAddBaseNumbers(kp, add.units = TRUE)
  
  kpDataBackground(kp, data.panel = 1)
  
  # kpAxis(kp, r0 = 1.55, r1 = 2, ymin = -5, ymax = 5, numticks = 11, side = 2)
  
  kpDataBackground(kp, data.panel = 1)
  
  data_to_visualize = meth_data[ ,patient[i]] - meth_data[ ,patient[i + 1]]
  
  neg = data_to_visualize > 0
  
  data_to_visualize[neg] = log2(31*data_to_visualize[neg] + 1) / 5
  data_to_visualize[!neg] = (-1) * log2(31*abs(data_to_visualize[!neg]) + 1) / 5
  
  rm(neg)
  
  col_to_apply = ifelse(data_to_visualize > 0, "tomato", "skyblue3")
  
  kpPoints(kp, chr = meth_data$chromosome_name, x = meth_data$end_position, y = ((data_to_visualize / 2) + 0.5), col = col_to_apply, cex = 0.7)
  kpAddLabels(kp, labels = "DNAMethylation", side = "left", srt = 90, label.margin = 0.03)
  
  dev.off()
}

rm(meth_data, pre_data, chr_to_visual, file_name, 
   i, patient, img_pre, img_post, img, kp, image_folder_name, data_to_visualize)
rm(col_to_apply, meth_color, scale)