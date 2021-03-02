# TADiff exclude

#' ########### Generating output images ########### 
#' dir.create(image_output_folder, showWarnings = FALSE)
#' 
#' #' Generating image with events in two groups
#' #' 
#' #' @param tad_to_visual TAD to visualize
#' 
#' tad_to_visual = c("TAD2130", "TAD854")
#' 
#' # tad_to_visual = c(tad_to_visual,
#' #                   tad_sign[which(tad_sign$FDR == min(tad_sign$FDR)), ]$tad_name)
#' # 
#' # tad_to_visual = c(tad_to_visual,
#' #                   tad_sign[which(tad_sign$mean == max(tad_sign$mean)), ]$tad_name)
#' 
#' 
#' for(i in tad_to_visual){
#'   tad.test = full %>% dplyr::filter(tad_name == i)
#'   tad.test.1 = as.matrix(unlist(tad.test[,..list1]))
#'   tad.test.1 = as.data.frame(tad.test.1)
#'   tad.test.1$status = paste(group1)
#'   tad.test.2 = as.matrix(unlist(tad.test[,..list2]))
#'   tad.test.2 = as.data.frame(tad.test.2)
#'   tad.test.2$status = paste(group2)
#'   
#'   tad.test.plot = rbind(tad.test.1, tad.test.2)
#'   
#'   # png(filename = paste(image_output_folder, "/", i, ".png", sep = ""), 
#'   #     width = 600, height = 820)
#'   
#'   # pdf(paste(image_output_folder, "/", i, ".pdf", sep = ""))
#'   
#'   gr = ggplot(tad.test.plot, aes(x = status, y = V1, fill = status)) + 
#'        geom_jitter(position = position_jitter(0.1), shape = 21, color = "gray61", size = 1.5) +
#'        theme_bw() + 
#'        theme(panel.grid.major = element_blank(), 
#'              panel.grid.minor = element_blank(), 
#'              legend.position = "none",
#'              axis.text.x = element_text(size = 15),
#'              axis.text.y = element_text(size = 15),
#'              axis.title.y = element_text(size = 16)) +
#'        labs(y = i, x = "") +
#'        stat_summary(fun = mean, fun.min = mean, fun.max = mean, 
#'                     geom = "crossbar", width = 0.2)
#'   
#'   saveImageHigh::save_as_pdf({print(gr)},
#'                              file.name = file.path(image_output_folder,
#'                                                    paste(i, "allValues.png", sep = "_")),
#'                              height = 8)
#'   
#'   # dev.off()
#'   
#'   tad.test = full %>% dplyr::filter(tad_name == i)
#'   tad.test = tad.test[order(abs(tad.test$diff), decreasing = TRUE), ]
#' 
#'   tad.test.1 = tad.test[,..list1]
#'   tad.test.2 = tad.test[,..list2]
#'   tad.test.1$mean = apply(tad.test.1, 1, mean)
#'   tad.test.2$mean = apply(tad.test.2, 1, mean)
#'   tad.test.1$status = paste(group1)
#'   tad.test.2$status = paste(group2)
#'   tad.test.1$ids = 1:nrow(tad.test.1)
#'   tad.test.2$ids = 1:nrow(tad.test.2)
#' 
#'   tad.test.1$xj = 1
#'   tad.test.2$xj = 2
#' 
#'   tad.test.1$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")
#'   tad.test.2$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")
#' 
#'   tad.test.1[1:min(30, nrow(tad.test.1)), ]$line.color = "gray38"
#' 
#'   tad.test.plot = rbind(tad.test.1[,c("mean", "xj", "ids", "line.color", "status")],
#'                         tad.test.2[,c("mean", "xj", "ids", "line.color", "status")])
#' 
#' 
#' 
#'   tad.test.plot$xj = jitter(tad.test.plot$xj, amount = 0.1, factor = 1)
#'   
#'   # png(filename = paste(image_output_folder, "/", i, "_topValues.png", sep = ""),
#'   #     width = 600, height = 820)
#'   
#'   # pdf(paste(image_output_folder, "/", i, "_topValues.pdf", sep = ""))
#'   
#'   f3 = ggplot(data = tad.test.plot, aes(y = mean, x = xj, fill = status)) +
#'        geom_line(aes(x = xj, group = ids), color = "lightgray") +
#'        geom_line(aes(x = xj, group = ids), color = tad.test.plot$line.color) +
#'        geom_point(aes(x = xj), size = 1.5, shape = 21, color = "gray61") +
#'     
#'        # geom_half_boxplot(data=tad.test.plot %>% filter(status == group1), 
#'        #                   aes(x=xj, y = mean), position = position_nudge(x = -.25),
#'        #                   side = "r",outlier.shape = NA, center = TRUE, 
#'        #                   errorbar.draw = FALSE, width = .2) +
#'        # 
#'        # geom_half_boxplot(data = tad.test.plot %>% filter(status==group2), 
#'        #                   aes(x=xj, y = mean), position = position_nudge(x = .15),
#'        #                   side = "r",outlier.shape = NA, center = TRUE, 
#'        #                   errorbar.draw = FALSE, width = .2) +
#'     
#'        geom_half_violin(data = tad.test.plot %>% filter(status==group1),
#'                         aes(x = xj, y = mean), position = position_nudge(x = -.3),
#'                         side = "l") +
#'     
#'        geom_half_violin(data = tad.test.plot %>% filter(status==group2),
#'                         aes(x = xj, y = mean), position = position_nudge(x = .3),
#'                         side = "r") +
#'     
#'        scale_x_continuous(breaks=c(1,2), labels=c("ss6", "ss8"), limits=c(0, 3)) +
#'        theme_classic() + labs(y = i, x = "") + 
#'        theme(panel.grid.major = element_blank(), 
#'              panel.grid.minor = element_blank(), 
#'              legend.position = "none",
#'              axis.text.x = element_text(size = 15),
#'              axis.text.y = element_text(size = 15),
#'              axis.title.y = element_text(size = 16))
#'   
#'   saveImageHigh::save_as_pdf({print(f3)},
#'                              file.name = file.path(image_output_folder,
#'                                                    paste(i, "_topValues.png", sep = "")),
#'                              height = 8)
#'                                
#'   
#'   # dev.off()
#'   
#'   tad.test = full %>% dplyr::filter(tad_name == i)
#'   tad.test = tad.test[order(tad.test$diff, decreasing = TRUE), ]
#' 
#'   tad.test.1 = tad.test[,..list1]
#'   tad.test.2 = tad.test[,..list2]
#'   tad.test.1$mean = apply(tad.test.1, 1, mean)
#'   tad.test.2$mean = apply(tad.test.2, 1, mean)
#'   tad.test.1$status = paste(group1)
#'   tad.test.2$status = paste(group2)
#'   tad.test.1$ids = 1:nrow(tad.test.1)
#'   tad.test.2$ids = 1:nrow(tad.test.2)
#' 
#'   tad.test.1$xj = 1
#'   tad.test.2$xj = 2
#'   
#'   tad.test.1 = tad.test.1[order(tad.test$diff, decreasing = TRUE), ]
#'   tad.test.2 = tad.test.2[order(tad.test$diff, decreasing = TRUE), ]
#' 
#'   tad.test.1$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")
#'   tad.test.2$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")
#' 
#'   tad.test.1[1:min(30, nrow(tad.test.1)), ]$line.color = "gray38"
#' 
#'   tad.test.plot = rbind(tad.test.1[,c("mean", "xj", "ids", "line.color", "status")],
#'                         tad.test.2[,c("mean", "xj", "ids", "line.color", "status")])
#' 
#' 
#' 
#'   tad.test.plot$xj = jitter(tad.test.plot$xj, amount = 0.1, factor = 1)
#' 
#'   # png(filename = paste(image_output_folder, "/", i, "_positive_connections.png", sep = ""),
#'   #     width = 600, height = 820)
#'   
#'   # pdf(paste(image_output_folder, "/", i, "_positive_connections.pdf", sep = ""))
#' 
#'   f3 = ggplot(data=tad.test.plot, aes(y = mean, x = xj, fill = status)) +
#'        geom_line(aes(x=xj, group=ids), color = "lightgray") +
#'        geom_line(aes(x=xj, group=ids), color = tad.test.plot$line.color) +
#'        geom_point(aes(x=xj), size = 1.5, shape = 21, color = "gray61") +
#'     
#'        geom_half_violin(data = tad.test.plot %>% filter(status==group1),
#'                         aes(x = xj, y = mean), position = position_nudge(x = -.3),
#'                         side = "l") +
#'     
#'        geom_half_violin(data = tad.test.plot %>% filter(status==group2),
#'                         aes(x = xj, y = mean), position = position_nudge(x = .3),
#'                         side = "r") +
#'     
#'        scale_x_continuous(breaks=c(1,2), labels=c("ss6", "ss8"), limits=c(0, 3)) +
#'        theme_classic() + labs(y = i, x = "") +
#'        theme(panel.grid.major = element_blank(),
#'              panel.grid.minor = element_blank(),
#'              legend.position = "none",
#'              axis.text.x = element_text(size = 15),
#'              axis.text.y = element_text(size = 15),
#'              axis.title.y = element_text(size = 16))
#'   
#'   saveImageHigh::save_as_pdf({print(f3)},
#'                              file.name = file.path(image_output_folder,
#'                                                    paste(i, "_topPositiveValues.png", sep = "")),
#'                              height = 8)
#' 
#'   # dev.off()
#'   
#'   tad.test = full %>% dplyr::filter(tad_name == i)
#'   tad.test = tad.test[order(tad.test$diff), ]
#' 
#'   tad.test.1 = tad.test[,..list1]
#'   tad.test.2 = tad.test[,..list2]
#'   tad.test.1$mean = apply(tad.test.1, 1, mean)
#'   tad.test.2$mean = apply(tad.test.2, 1, mean)
#'   tad.test.1$status = paste(group1)
#'   tad.test.2$status = paste(group2)
#'   tad.test.1$ids = 1:nrow(tad.test.1)
#'   tad.test.2$ids = 1:nrow(tad.test.2)
#' 
#'   tad.test.1$xj = 1
#'   tad.test.2$xj = 2
#'   
#'   tad.test.1 = tad.test.1[order(tad.test$diff), ]
#'   tad.test.2 = tad.test.2[order(tad.test$diff), ]
#'   
#'   tad.test.1$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")
#'   tad.test.2$line.color = rgb(255, 255, 255, max = 255, alpha = 0, names = "white")
#'   
#'   tad.test.1[1:min(30, nrow(tad.test.1)), ]$line.color = "gray38"
#'   
#'   tad.test.plot = rbind(tad.test.1[,c("mean", "xj", "ids", "line.color", "status")],
#'                         tad.test.2[,c("mean", "xj", "ids", "line.color", "status")])
#'   
#'   
#'   
#'   tad.test.plot$xj = jitter(tad.test.plot$xj, amount = 0.1, factor = 1)
#'   
#'   # png(filename = paste(image_output_folder, "/", i, "_negative_connections.png", sep = ""),
#'   #     width = 600, height = 820)
#'   
#'   # pdf(paste(image_output_folder, "/", i, "_negative_connections.pdf", sep = ""))
#'   
#'   f3 = ggplot(data=tad.test.plot, aes(y = mean, x = xj, fill = status)) +
#'        geom_line(aes(x=xj, group=ids), color = "lightgray") +
#'        geom_line(aes(x=xj, group=ids), color = tad.test.plot$line.color) +
#'        geom_point(aes(x=xj), size = 1.5, shape = 21, color = "gray61") +
#'           
#'        geom_half_violin(data = tad.test.plot %>% filter(status==group1),
#'                         aes(x = xj, y = mean), position = position_nudge(x = -.3),
#'                         side = "l") +
#'         
#'        geom_half_violin(data = tad.test.plot %>% filter(status==group2),
#'                         aes(x = xj, y = mean), position = position_nudge(x = .3),
#'                         side = "r") +
#'         
#'        scale_x_continuous(breaks=c(1,2), labels=c("ss6", "ss8"), limits=c(0, 3)) +
#'        theme_classic() + labs(y = i, x = "") + 
#'        theme(panel.grid.major = element_blank(), 
#'              panel.grid.minor = element_blank(), 
#'              legend.position = "none",
#'              axis.text.x = element_text(size = 15),
#'              axis.text.y = element_text(size = 15),
#'              axis.title.y = element_text(size = 16))
#'   
#'   saveImageHigh::save_as_pdf({print(f3)},
#'                              file.name = file.path(image_output_folder,
#'                                                    paste(i, "_topNegativeValues.png", sep = "")),
#'                              height = 8)
#'   
#'   # dev.off()
#' }
#' 
#' ########### Clear enviroment ########### 
#' 
#' end_tad_time = Sys.time()
#' 
#' rm(list = setdiff(ls(), c("data.all", "full", "full.tads", "tad_sign", "tad_sum", 
#'                         "dir_name", "end_tad_time", "start_tad_time",
#'                         "paired_data", "groups", "FDR_criterion", "genes.found", "output_folder", "x", "res", "res2")))
#' 
#' tad_sign = tad_sign[,c("tad_name", "count", "mean", "FDR")]
#' 
#' ########### Generating outputs ###########  
#' dir.create(output_folder, showWarnings = FALSE)
#' 
#' write.table(full, paste(output_folder, "/integrated_table_with_tads.csv", sep = ""),
#'             row.names = FALSE, sep = "\t", quote = FALSE)
#' 
#' if(!is.null(groups)){
#' 
#'   write.table(tad_sum, paste(output_folder, "/tad_statistics.csv", sep = ""),
#'               row.names = FALSE, sep = "\t", quote = FALSE)
#' 
#'   write.table(full.tads, paste(output_folder, "/integrated_table_with_sign_tads.csv", sep = ""),
#'               row.names = FALSE, sep = "\t", quote = FALSE)
#' 
#'   write.table(tad_sign, paste(output_folder, "/sign_tad_statistics.csv", sep = ""),
#'               row.names = FALSE, sep = "\t", quote = FALSE)
#' 
#'   write.table(genes.found, paste(output_folder, "/genes_found.txt", sep = ""),
#'               row.names = FALSE, col.names = FALSE, quote = FALSE)
#' 
#' 
#' }