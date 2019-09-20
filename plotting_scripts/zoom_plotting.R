## Getting required libraries

library(stringr)
library(stringi)
library(karyoploteR)
library(data.table)
library(png)
library(magick)

##################### Inputs ##################

# file to read
file_name = "biodata_pre_integration_14.csv"
pre_data = read.csv(file_name, sep = ";")

image_folder_name = "images_chr14_MOAP1"
image_folder_name = paste(image_folder_name, "/", sep = "")

# chromosomes to visualize
chr_to_visual = "chr14"

# patient to visualize
patient = colnames(pre_data)[8:length(pre_data)]
# patient = which(colnames(pre_data) %in% patient)

# specific region to zoom
start = 93648541
end = 93651273
tickDistance = round((end - start) / 7)
region = toGRanges(data.frame(chr_to_visual, start, end))
# rm(start, end)

########################## Reading data ##########################

meth_data = pre_data[which(str_detect(pre_data$ID, "cg")), ]
meth_data$chromosome_name = paste("chr", meth_data$chromosome_name, sep = "")

gene_data = pre_data[which(str_detect(pre_data$ID, "ENSG")), ]
gene_data$chromosome_name = paste("chr", gene_data$chromosome_name, sep = "")

temp = gene_data[!(gene_data$start_position > end), ]
gene_data = temp[!(temp$end_position < start), ]

rm(start, end)

rm(temp)

var_data = pre_data[is.na(pre_data$ID), ]
var_data$chromosome_name = paste("chr", var_data$chromosome_name, sep = "")

var_data[,8:ncol(var_data)] = ifelse(is.na(var_data[,8:ncol(var_data)]), 0, 1)

temp = unique(var_data[,1:4])
temp = split(temp, seq(nrow(temp)))

temp = lapply(temp, function(obj, data){
  data = data[which(data$chromosome_name == obj$chromosome_name), ]
  data = data[which(data$start_position == obj$start_position), ]
  data = data[which(data$end_position == obj$end_position), ]
  
  varAn = paste(data$VarAnnotation, collapse = ";")
  numericData = colSums(data[ ,8:ncol(data)])
  
  varAn = data.frame(VarAnnotation = varAn)
  rownames(varAn) = rownames(obj)
  numericData = t(as.data.frame(numericData))
  rownames(numericData) = rownames(obj)
  
  out = cbind(obj, varAn, numericData)
  
  return(out)
}, var_data)

var_data = rbindlist(temp)
var_data = as.data.frame(var_data)

rm(temp)

######################## Plotting #########################

plot.params = getDefaultPlotParams(plot.type = 1)
plot.params$data1height = 200
plot.params$bottommargin = 10
plot.params$data1inmargin = 5
plot.params$rightmargin = 0.02
plot.params$leftmargin = 0.04
plot.params$ideogramheight = 10
plot.params$topmargin = 0

for(i in patient){
  png(filename = paste(image_folder_name, i, ".png", sep = ""), width = 2494, height = 1812)
  
  kp = plotKaryotype(genome="hg19", plot.type = 1, plot.params = plot.params, zoom = region)
  kpAddBaseNumbers(kp, tick.dist = tickDistance, add.units = TRUE)
  
  kpDataBackground(kp, data.panel = 1, r0=0, r1=0.23)
  # kpAxis(kp, r0 = 0, r1 = 0.45, ymin = 0, ymax = 100)
  meth_color = character(nrow(meth_data))
  meth_color = ifelse(meth_data[ ,i] <= 30, "green4", ifelse(meth_data[ ,i] <= 70, "blue", "red"))
  kpPoints(kp, chr = meth_data$chromosome_name, r0 = 0, r1 = 0.23, x = meth_data$end_position, y = (meth_data[ ,i] / 100), col = meth_color, cex = 1.5)
  kpAddLabels(kp, labels = "DNAMethylation", side = "left", r0 = 0, r1 = 0.23, srt = 90, label.margin = 0.03)
  
  kpDataBackground(kp, data.panel = 1, r0=0.25, r1=0.48)
  kpAxis(kp, r0 = 0.25, r1 = 0.48, ymin = 0, ymax = max(gene_data[ ,i]))
  kpBars(kp, chr = gene_data$chromosome_name, 
         r0 = 0.25, r1 = 0.48, 
         x0 = gene_data$start_position, x1 = gene_data$end_position, 
         y1 = (gene_data[ ,i] / max(gene_data[ ,i])), border = "slateblue4", cex = 1, lwd = 3)
  kpAddLabels(kp, labels = "RNASeq", side = "left", r0 = 0.25, r1 = 0.48, srt = 90, label.margin = 0.03)
  
  kpDataBackground(kp, data.panel = 1, r0=0.5, r1=0.73)
  kpAxis(kp, r0 = 0.5, r1 = 0.73, ymin = 0, ymax = 5, numticks = 6)
  kpSegments(kp, chr = var_data$chromosome_name,
             r0 = 0.5, r1 = 0.73,
             x0 = var_data$start_position,
             x1=var_data$start_position,
             y0 = 0, y1=(var_data[ ,i] / 5), lwd = 2)
  kpAddLabels(kp, labels = "FCR_WES", side = "left", r0 = 0.5, r1 = 0.73, srt = 90, label.margin = 0.03)
  
  kpDataBackground(kp, data.panel = 1, r0=0.75, r1=1)
  kpAxis(kp, r0 = 0.75, r1 = 1, ymin = 0, ymax = max(gene_data[ ,i]))
  kpAxis(kp, r0 = 0.75, r1 = 1, ymin = 0, ymax = 5, numticks = 6, side = 2)
  kpPoints(kp, chr = meth_data$chromosome_name, r0 = 0.75, r1 = 1, x = meth_data$end_position, y = (meth_data[ ,i] / 100), col = meth_color, cex = 1.5)
  kpBars(kp, chr = gene_data$chromosome_name, 
         r0 = 0.75, r1 = 1, 
         x0 = gene_data$start_position, x1 = gene_data$end_position, 
         y1 = (gene_data[ ,i] / max(gene_data[ ,i])), border = "slateblue4", lwd = 3)
  kpSegments(kp, chr = var_data$chromosome_name,
             r0 = 0.75, r1 = 1,
             x0 = var_data$start_position,
             x1=var_data$start_position,
             y0 = 0, y1=(var_data[ ,i] / 5), lwd = 2)
  kpAddLabels(kp, labels = "Integrated table", side = "left", r0 = 0.75, r1 = 1, srt = 90, label.margin = 0.03)
  
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
  
  kp = plotKaryotype(genome="hg19", plot.type = 1, plot.params = plot.params, zoom = region)
  kpAddBaseNumbers(kp, tick.dist = tickDistance, add.units = TRUE)
  
  kpDataBackground(kp, data.panel = 1, r0=0.75, r1=1)
  
  kpAxis(kp, r0 = 0.75, r1 = 1, ymin = -5, ymax = 5, numticks = 11, side = 2)
  
  kpDataBackground(kp, data.panel = 1, r0=0, r1=0.23)
  
  data_to_visualize = meth_data[ ,patient[i]] - meth_data[ ,patient[i + 1]]
  col_to_apply = ifelse(data_to_visualize > 0, "tomato", "skyblue3")
  
  kpPoints(kp, chr = meth_data$chromosome_name, r0 = 0, r1 = 0.23, x = meth_data$end_position, y = ((data_to_visualize / 200) + 0.5), col = col_to_apply, cex = 1.5)
  kpPoints(kp, chr = meth_data$chromosome_name, r0 = 0.75, r1 = 1, x = meth_data$end_position, y = ((data_to_visualize / 200) + 0.5), col = col_to_apply, cex = 1.5)
  kpAddLabels(kp, labels = "DNAMethylation", side = "left", r0 = 0, r1 = 0.23, srt = 90, label.margin = 0.03)
  
  data_to_visualize = gene_data[ ,patient[i]] - gene_data[ ,patient[i + 1]]
  col_to_apply = ifelse(data_to_visualize > 0, "red4", "blue4")
  
  kpDataBackground(kp, data.panel = 1, r0=0.25, r1=0.48)
  kpAxis(kp, r0 = 0.25, r1 = 0.48, ymin = (min(data_to_visualize)), ymax = (max(data_to_visualize)))
  kpAxis(kp, r0 = 0.75, r1 = 1, ymin = (min(data_to_visualize)), ymax = (max(data_to_visualize)))
  
  Y0 = ifelse(data_to_visualize > 0, 0, data_to_visualize) # / abs(min(data_to_visualize))
  Y1 = ifelse(data_to_visualize > 0, data_to_visualize, 0) # / (max(data_to_visualize))
  
  Y0 = Y0 * 0.5 / abs(min(Y0))
  Y1 = Y1 * 0.5 / max(Y1)
  
  kpRect(kp, chr = gene_data$chromosome_name,
         r0 = 0.25, r1 = 0.48,
         x0 = gene_data$start_position, x1 = gene_data$end_position,
         y0 = (Y0 + 0.5),
         y1 = (Y1 + 0.5), border = col_to_apply, lwd = 2)
  kpRect(kp, chr = gene_data$chromosome_name,
         r0 = 0.75, r1 = 1,
         x0 = gene_data$start_position, x1 = gene_data$end_position,
         y0 = (Y0 + 0.5),
         y1 = (Y1 + 0.5), border = col_to_apply, lwd = 2)
  kpAddLabels(kp, labels = "RNASeq", side = "left", r0 = 0.25, r1 = 0.48, srt = 90, label.margin = 0.03)
  
  kpDataBackground(kp, data.panel = 1, r0=0.5, r1=0.73)
  kpAxis(kp, r0 = 0.5, r1 = 0.73, ymin = -5, ymax = 5, numticks = 11)
  data_to_visualize = var_data[ ,patient[i]] - var_data[ ,patient[i + 1]]
  
  Y0 = ifelse(data_to_visualize > 0, 0, data_to_visualize) * 0.5 / 5
  Y1 = ifelse(data_to_visualize > 0, data_to_visualize, 0) * 0.5 / 5
  
  col_to_apply = ifelse(data_to_visualize > 0, "red", "blue")
  
  kpSegments(kp, chr = var_data$chromosome_name,
             r0 = 0.5, r1 = 0.73,
             x0 = var_data$start_position,
             x1=var_data$start_position,
             y0 = Y0 + 0.5, y1 = Y1 + 0.5,
             col = col_to_apply, lwd = 2)
  kpAbline(kp, chr = gene_data$chromosome_name,
           h = 0.5, r0 = 0.5, r1 = 0.73)
  kpSegments(kp, chr = var_data$chromosome_name,
             r0 = 0.75, r1 = 1,
             x0 = var_data$start_position,
             x1=var_data$start_position,
             y0 = Y0 + 0.5, y1 = Y1 + 0.5,
             col = col_to_apply, lwd = 2)
  kpAddLabels(kp, labels = "FCR_WES", side = "left", r0 = 0.5, r1 = 0.73, srt = 90, label.margin = 0.03)
  kpAddLabels(kp, labels = "Integrated table", side = "left", r0 = 0.75, r1 = 1, srt = 90, label.margin = 0.03)
  dev.off()
}

rm(gene_data, meth_data, plot.params, pre_data, var_data, chr_to_visual, file_name, 
   i, patient, img_pre, img_post, img, kp, image_folder_name, data_to_visualize)
rm(col_to_apply, meth_color, Y0, Y1, region)