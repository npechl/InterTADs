# Getting required libraries ---------------------------------------------------

source("libraries.R")

# Inputs -----------------------------------------------------------------------

#'
#' Input parameters for Visualization part
#' 
#' @param dir_name Directory of input datasets 
#' 
#' @param tad_folder Folder name for reading InterTADs output tables
#' 
#' @param meta meta-data file name used
#' 
#' @param image_folder_name Directory name of output visualizations
#' 
#' @param diff_col Group column to compare
#' 
#' @param my_colors Colors for printing images
#' 
#' @param my_cex Sizes for printing images
#' 

dir_name <- "Datasets"

tad_folder <- "output_tables_test"

meta <- "meta-data.csv"

image_folder_name <- "output_visualizations_test"

diff_col <- "groups"


my_colors <- c("black", "slateblue4", "red", "goldenrod4", "darkorange3")
my_cex <- c(1, 0.8, 1.2, 1.4)


full.tads <- fread(paste(tad_folder, "/integrated_table_with_tads.csv",
                        sep = ""))
tad.sign <- fread(paste(tad_folder, "/sign_tad_statistics.csv", sep = ""))

meta <- fread(paste(dir_name, meta, sep = "/"))
who <- meta == ""
who <- apply(who, 1, sum)
meta <- meta[which(who == 0), ]

groups <- unique(meta[[diff_col]])

dir.create(image_folder_name, showWarnings = FALSE)
# image_folder_name = paste(image_folder_name, "/", sep = "")

tad_to_visual <- c("TAD2")
# tad_to_visual = c(tad_to_visual,
#                   tad.sign[which(tad.sign[,"mean"] == max(tad.sign[,"mean"])),
#                   ]$tad_name)
# 
# tad_to_visual = c(tad_to_visual,
#                   tad.sign[which(tad.sign[,"FDR"] == min(tad.sign[,"FDR"])),
#                   ]$tad_name)

# For every sign TAD ----------------------------------------------------------------------------

for(j in tad_to_visual){
    full.vtads <- full.tads[which(full.tads$tad_name %in% j), ]
    full.vtads$chromosome_name <- paste("chr", full.vtads$chromosome_name,
                                        sep = "")
    
    chr_to_visual <- as.character(unique(full.vtads$chr))
    start <- min(unique(full.vtads$start))
    end <- max(unique(full.vtads$end))
    region <- toGRanges(data.frame(chr_to_visual, start, end))
    tick.distance <- (end - start)/6
    # tick.distance = tick.distance + (tick.distance / 5.9)
    
    ## Reading data ------------------------------------------------------------
    
    columns <- c("ID", "chromosome_name", "start_position", "end_position",
                  meta$newNames, "diff")
    
    full.vtads$range <- full.vtads$end_position - full.vtads$start_position
    
    ## Plotting ----------------------------------------------------------------
  
    dir.create(paste(image_folder_name, "/patients/", sep = ""),
                showWarnings = FALSE)
    
    parents <- unique(full.vtads$parent)
    
    for(i in meta$newNames){
        # png(filename = paste(image_folder_name, "patients/", i, "_",
        #j, ".png", sep = ""), width = 1150, height = 900)
        
        kp <- expression({
        kp <- plotKaryotype(genome = "hg19", plot.type = 4, zoom = region,
                            cex = 1.5)
        kpDataBackground(kp, data.panel = 1)
        kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = tick.distance,
                        cex = 1.25)
        lb <- c("0 %", "50 %", "100 %")
        kpAxis(kp, r0 = 0, r1 = 1, ymin = 0, ymax = 100, side = 1, labels = lb)
          
        for(p in 1:length(parents)){
            vis_data <- full.vtads[which(full.vtads$parent == parents[p]), ]
            
            numeric.vector <- vis_data[[i]] 
            numeric.vector <- numeric.vector / 100
            
            if(max(vis_data$range) <= 1){
                (kp, 
                 chr = vis_data$chromosome_name,
                 x = vis_data$end_position, 
                 y = numeric.vector, 
                 col = my_colors[p],
                 cex = my_cex[p])
            } 
            else {
                kpBars(kp, 
                 chr = vis_data$chr,
                 x0 = vis_data$start_position, 
                 x1 = vis_data$end_position,
                 y1 = numeric.vector, 
                 border = my_colors[p], 
                 cex = my_cex[p])
              
            }
            
            kpAddMainTitle(kp, main = j, cex = 1.5)
            
        }
        })
        
        saveImageHigh::save_as_png({eval(kp)},
                                    file.name = file.path(image_folder_name, 
                                                            "patients", 
                                                            paste(i, "_", j,
                                                                  ".png",
                                                                  sep = "")),
                                    width = 14, height = 10)
        
        # dev.off()
        }
    
    dir.create(paste(image_folder_name, j, sep = "/"), showWarnings = FALSE)
    
    ## Filtering data ----------------------------------------------------------
    
    columns <- c("ID", "chromosome_name", "start_position", "end_position",
                meta$newNames, "diff", "parent")
    
    full.vtads$range <- full.vtads$end_position - full.vtads$start_position
  
    box_data <- full.vtads[which(full.vtads$range > 1), ..columns]
    
    
    rain_data <- full.vtads[which(full.vtads$range <= 1), ..columns]
  
    box_data <- box_data %>% dplyr::group_by(start_position, end_position) %>%
        dplyr::filter(abs(diff) == max(abs(diff)))
    box_data <- as.data.table(box_data)
  
    
    ## Plotting ----------------------------------------------------------------
    
    if(!is.null(groups)){
      
        for(i in 1:length(groups)){
            keep <- meta$newNames[which(meta[[diff_col]] == groups[i])]
            
            box_data$mean <- apply(box_data[,..keep], 1, mean, na.rm = TRUE)
            rain_data$mean <- apply(rain_data[,..keep], 1, mean, na.rm = TRUE)
            
            maximum_value <- 100
            minimum_value <- 0
            
            # png(filename = paste(image_folder_name, j, "/", groups[i], ".png",
            #sep = ""), 
            #     width = 1150, height = 900)
            
            kp <- expression({
            kp = plotKaryotype(genome = "hg19", plot.type = 4, zoom = region,
                                cex = 1.5)
            kpDataBackground(kp, data.panel = 1)
            kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = tick.distance,
                                cex = 1.25)
            
            lb <- c(minimum_value, (maximum_value + minimum_value)/2,
                    maximum_value)
            lb <- round(lb, digits = 2)
            lb <- paste(lb, "%", sep = " ")
            kpAxis(kp, r0 = 0, r1 = 1, ymin = minimum_value,
                    ymax = maximum_value, side = 1, labels = lb)
            
            rain_data <- rain_data[which(rain_data$mean <= maximum_value &
                                            rain_data$mean >= minimum_value), ]
            box_data <- box_data[which(box_data$mean <= maximum_value &
                                            box_data$mean >= minimum_value), ]
            
            for(p in 1:length(parents)){
              
                r_vis_data <- rain_data[which(rain_data$parent == parents[p]), ]
                b_vis_data <- box_data[which(box_data$parent == parents[p]), ]
                
                if(nrow(r_vis_data) > 0){
                    numeric.vector <- r_vis_data$mean
                    numeric.vector <- numeric.vector / maximum_value
                    kpPoints(kp, 
                            chr = r_vis_data$chromosome_name,
                            x = r_vis_data$end_position, 
                            y = numeric.vector, 
                            col = my_colors[p],
                            cex = my_cex[p])
                  }
                
                
                if(nrow(b_vis_data) > 0){
                    numeric.vector <- b_vis_data$mean
                    numeric.vector <- numeric.vector / maximum_value
                    kpBars(kp, 
                            chr = b_vis_data$chr,
                            x0 = b_vis_data$start_position, 
                            x1 = b_vis_data$end_position,
                            y1 = numeric.vector, 
                            border = my_colors[p], 
                            cex = my_cex[p])
                }
        
            }
            
            kpAddMainTitle(kp, main = j, cex = 1.5)
            })
            
            saveImageHigh::save_as_png({eval(kp)},
                                        file.name = file.path(image_folder_name,
                                                                j,
                                                                paste(groups[i],
                                                                ".png", 
                                                                sep = "")),
                                        width = 14, height = 10)
            # dev.off()
        }
        
        # png(filename = paste(image_folder_name, j, "/", "Group_differences",
        #".png", sep = ""), 
        #     width = 1150, height = 900)
        
        kp <- expression({
          
          
            kp <- plotKaryotype(genome = "hg19", plot.type = 4, zoom = region,
                                cex = 1.5)
            kpDataBackground(kp, data.panel = 1)
            kpAddBaseNumbers(kp, add.units = TRUE, tick.dist = tick.distance,
                                cex = 1.25)
            
            # all = c(rain_data$diff, box_data$diff)
            
            maximum_value <- 100
            minimum_value <- -100
            
            lb <- c(minimum_value, 0, maximum_value)
            lb <- round(lb, digits = 2)
            lb <- paste(lb, "%", sep = " ")
            kpAxis(kp, r0 = 0, r1 = 1, ymin = minimum_value,
                    ymax = maximum_value, side = 1, labels = lb)
            
            for(p in 1:length(parents)){
              
                r_vis_data <- rain_data[which(rain_data$parent == parents[p]), ]
                b_vis_data <- box_data[which(box_data$parent == parents[p]), ]
                
                if(nrow(r_vis_data) > 0){
                    numeric.vector <- r_vis_data$diff
                    who <- numeric.vector > 0
                    numeric.vector[who] <- numeric.vector[who] / maximum_value
                    numeric.vector[!who] <- numeric.vector[!who] /
                        abs(minimum_value)
                    numeric.vector <- (numeric.vector + 1)/2
                    
                    kpPoints(kp, 
                                chr = r_vis_data$chromosome_name,
                                x = r_vis_data$end_position, 
                                y = numeric.vector, 
                                col = my_colors[p],
                                cex = my_cex[p])
                  }
                
                if(nrow(b_vis_data) > 0){
                    numeric.vector <- b_vis_data$diff
                    who <- numeric.vector > 0
                    numeric.vector[who] <- numeric.vector[who] / maximum_value
                    numeric.vector[!who] <- numeric.vector[!who] /
                        abs(minimum_value)
                    numeric.vector <- (numeric.vector + 1)/2
                    
                    bottom <- rep(0.5, length(numeric.vector))
                    bottom[!who] <- numeric.vector[!who]
                    
                    up <- rep(0.5, length(numeric.vector))
                    up[who] <- numeric.vector[who]
                    
                    kpBars(kp, 
                            chr = b_vis_data$chromosome_name,
                            x0 = b_vis_data$start_position, 
                            x1 = b_vis_data$end_position,
                            y0 = bottom, 
                            y1 = up,
                            border = my_colors[p], 
                            cex = my_cex[p])
                }
                
            }
            
            kpAddMainTitle(kp, main = j, cex = 1.5)
            
        })
        
        saveImageHigh::save_as_png({eval(kp)},
                                    file.name = file.path(image_folder_name, j,
                                        "Group_differences.png"),
                                        width = 14, height = 10)
        
        # dev.off()
        
    }

}


rm(list = ls())
