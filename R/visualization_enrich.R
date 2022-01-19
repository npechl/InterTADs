#' This file contains the functions used in the "Visualization"
#' section of the `03_functional_analysis.R` script.


#' RStudio doesn't support all graph fonts in Windows OS.
#' This function sets the appropriate fonts for the graphs. 
set_graph_fonts <- function(device_system) {
  
  if (device_system == "win") {
    
    font_import() 
    loadfonts(device = "win")
    windowsFonts(Times = windowsFont("TT Times New Roman"))
    theme_set(theme_bw(base_size = 12, base_family = 'Times New Roman') + 
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()))
  }
}


#' This function is being called by the `03_functional_nalysis.R` script
#' 
#' It creates graphs of the Motif EA.
motif_visual <- function(image_output_folder,
                         motif_output_folder, 
                         data_visual, 
                         report_list, 
                         criterio) {
  
  data_density <- per_tad_plots(report_list, image_output_folder)
  
  per_tfs_plots(report_list, data_visual, image_output_folder)
  
  density_plot(data_density, image_output_folder)
  
  rm(data_density)
  
}



#' This function is being called by the `03_functional_analysis.R` 
#' script.
#' 
#' It creates graphs of the EA using Enrichr
enrichr_visual <- function(folder,
                           type,
                           data_visual, 
                           criterio) {
  
  # Set or/and create image output folder
  if (str_detect(type, "GO")) {
    
    string <- paste0(folder, "/", type)
    dir.create(string, showWarnings = FALSE)
    
  } else {
    
    string <- folder
  }
  
  compare_density_plot(string, type, data_visual)
  
  histogram_plot(string, type, data_visual)
  
  network_plot(string, type, data_visual) 
  
  group_plots(string, type, data_visual)
  
}


#' This fuction is called from the `enrichr_visual` function.
#' It creates two density plots in the same graph.
#' 
#' @details
#' The first density plot is made with the P values of the data 
#' after the analysis.The second density plot is made with the 
#' Adjusted P values.
compare_density_plot <- function(string,
                                 type,
                                 data_visual) {
  
  data_plot <- data.table(p_value = c(data_visual$p_value, data_visual$p_adjust),
                          Type = c(rep_len("P value", nrow(data_visual)),
                                  rep_len("Adjusted P value", nrow(data_visual))))
 
  p <- ggplot(data = data_plot, aes(x = p_value, group = Type, fill = Type)) +
    geom_density(adjust = 1.5, alpha = .4) +
    ggtitle(paste0("P values of " , type)) +
    xlab("P value") +
    ylab("Density") +
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 22),
          axis.title.x = element_text(size = 18, vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size = 18, vjust = 1.5, hjust = 0.5),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16)
    )

  save_image(print(p), 
             file.name = paste0(string, "/", "Density plot-P values of ", type, ".png"),
             res = 200)
  
  rm(data_plot, p)
}


#' This function is called from the `enrichr_visual` function.
#' It creates two histogram graphs.
#' 
#' @description 
#' The first one shows the Terms that were found in the most TADs 
#' in descending order according to the number of TADs they were 
#' found into.
#' 
#' The second one shows the same Terms, but adds the information of 
#' the mean Adjusted P value of each Term to be found in these TADs. 
#' The Terms are arranged in ascending order, according to the mean 
#' Adjusted P value.The bar color corresponds to the number of TADs 
#' each Term was found into.
histogram_plot <- function(string, 
                           type,
                           data_visual) {
  
  terms <- dplyr::count(data_visual, Term)
  
  data_plot <- merge(terms, data_visual) %>%
    dplyr::select(Term, n, p_adjust) %>%
    group_by(Term, n) %>%
    summarise(Term, n, p_adjust = mean(p_adjust)) %>%
    unique() 
  
  data_plot <- setorder(data_plot, -n) 
  
  if (nrow(data_plot) > 19)  data_plot <- data_plot[1:20, ]
  
  p <- data_plot %>%
    ggplot(aes(x = as.factor(reorder(Term, n)), y = n)) +
    geom_histogram(fill = "#66ccff", color = "#e9ecef", alpha = 0.9, stat = "identity") +
    labs(title = paste0("Top ", type, " in different TADs")) +
    xlab(type) +
    ylab("Number of TADs") +
    coord_flip() +
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.title.x = element_text(size = 20, vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size = 20, vjust = 0.5, hjust = 0.5)
    )
 
  save_image(print(p), 
             file.name = paste0(string, "/",type," in different TADs", ".png"),
             height = (0.5 * nrow(data_plot) + 2), 
             width = 22,
             res = 100)
  
  rm(p)
  
  #plot 2.5
  colnames(data_plot) <- str_replace(colnames(data_plot), "n", "number_of_tads")
 
  p <- data_plot %>%
    ggplot(aes(x = as.factor(reorder(Term,  -p_adjust)), 
               y = p_adjust, 
               fill = number_of_tads)) +
    geom_histogram(color = "#e9ecef", stat = "identity") +
    labs(title = paste0("Top ", type, " in different TADs")) +
    xlab(type) +
    ylab("Mean p_adjust") +
    coord_flip() +
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 28),
          plot.subtitle = element_text(size = 22),
          axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.title.x = element_text(size = 28, vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size = 28, vjust = 0.5, hjust = 0.5),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 14)
    ) +
    guides(fill = guide_legend(title = "number of TADs"))
  
  save_image(print(p), 
             file.name = paste0(string, "/Top ",type, ".png"),
             height = (0.5 * nrow(data_plot)) + 2,
             width = 22,
             res = 100)
  
  rm(p, data_plot, terms)
  
}


#' This function is called from the `enrichr_visual` function.
#' It creates a network graph of the Terms (max 10), 
#' that were found in the most TADs.
#' 
#' @details 
#' Each Term is a group.
#' Each node is unique and corresponds to a TAD.
#' The nodes are connected with lines. 
#' Each different line color corresponds to a group.
network_plot <- function(string, 
                         type,
                         data_visual) {
  
    # Choose only 10 largest groups
    go_count <- dplyr::count(data_visual, Term)
    go_count <- setorder(go_count, -n)
    go_count <- go_count[1:10, ]
    data_plot <- merge(data_visual, go_count, by = "Term")
    
    if (nrow(data_plot) == 10)  return(NULL)
    
    nodes <- data_plot %>%
      dplyr::select(TAD, Term)
    
    gos <- group_by(nodes, Term)
    groups <- group_split(gos)
    
    edges <- data.frame(from = character(), 
                        to = character() , 
                        group = character())
    
    for (i in c(1:length(groups))) {
      
      temp <- groups[[i]]
     
      for (j in c(1:nrow(temp))) {
        from_to <- data.frame(from = temp$TAD[j], 
                              to =temp$TAD[1:nrow(temp)],
                              group = temp$Term[1]) 
        from_to <- from_to[-j, ]
        edges <- dplyr::bind_rows(edges, from_to)
      }
    }
    
    edges <- unique(edges)
    
    if (nrow(edges) == 0)  return(NULL)
    
    g <- graph_from_data_frame(edges, directed = FALSE)
    
   
    network <- ggraph(g) +
      geom_edge_link(aes(color = group), alpha = 1, linemitre = 10, n = 1000) +     
      geom_node_point(size = 5, shape = 21, stroke = 1,
                      fill = 'white', color = 'black') +
      geom_node_text(aes(label = name), repel = TRUE, size = 5 , fontface = "bold") +
      theme_void() +
      ggtitle(label = paste0("Top ", type, " Network Graph")) +
      theme(
        plot.title = element_text(size = 20, hjust = 0.5, vjust = 0.5, face = "bold"),
        legend.text = element_text(size = 14),
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size = 14)
      )
    
    save_image(print(network), 
               file.name = paste0(string, "/Top ",type," network graph", ".png"),
               height = 7, width = 17, res = 100)

    rm(g, network, edges, from_to, temp, nodes, groups, gos,
       data_plot, go_count)
}


#' This function is called by the `enrichrVisual` function.
#' It creates two folders with histogram graphs.
#'
#' @description 
#' The first folder is called "Adjusted P values per TADs".
#' Each histogram shows the Terms found in each TAD in ascending order,
#' according to their Adjusted P value.
#'
#' The second folder is called "Adjusted P values per Term".
#' Each histogram shows the TADs, each Term was found into,
#' in ascending order, according to their Adjusted P value.
group_plots <- function(string,
                        type, 
                        data_visual) {
  
  if (str_detect(type,"GO")) {
    
    data_visual$ID <- paste0("GO ", data_visual$ID)
    data_plot <- data_visual
    
  } else if (str_detect(type,"KEGG")) {
    
    data_plot <- data_visual
    data_plot$ID <- data_plot$Term
  }
  
  tad_grouping <- data_plot %>%
    group_by(Term)
  new_grouping <- group_split(tad_grouping)
  
  dir.create(paste0(string, "/Adjusted P values per ", type), showWarnings = FALSE)
  
  iterations <- c(1:length(new_grouping))
  
  for (i in iterations) {
    
    temp <- new_grouping[[i]]
    if (str_detect(type,"GO", negate = TRUE))   temp$Term <- rep_len("", nrow(temp))
    
    temp$ID <- str_replace(temp$ID, "/", "-")
    
    p <- temp %>%
      ggplot(aes(x = as.factor(reorder(TAD, -p_adjust)), y = p_adjust)) +
      geom_histogram(fill="#66ccff", color="#e9ecef", stat = "identity", width = 0.8) +
      labs(title = paste0("Adjusted P values of ", temp$ID[1], " in different TADs"), 
           subtitle = temp$Term[1]) +
      ylab("Adjusted P value")+
      xlab("TAD number") +
      theme_bw() + 
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11),
            axis.title.x = element_text(size = 13, hjust = 0.5, vjust = 0.5, face = "bold"),
            axis.title.y = element_text(size = 13, hjust = 0.5, vjust = 0.5, face = "bold"),
            axis.text.x = element_text(size = 13, color = "black"),
            axis.text.y = element_text(size = 13, color = "black")
      ) +
      coord_flip() 

    save_image(print(p), 
               file.name = paste0(string, "/Adjusted P values per ", type, 
                                  "/Adjusted P values of ",temp$ID[1], ".png"),
               height = (0.5 * nrow(temp) + 2),
               width = 7, res = 100)
  }
  
 
  tad_grouping <- data_visual %>%
    group_by(TAD)
  new_grouping <- group_split(tad_grouping)
  
  dir.create(paste0(string, "/Adjusted P values per TAD"), showWarnings = FALSE)
  
  iterations <- c(1:length(new_grouping))
  
  for (i in iterations) {
    
    temp <- new_grouping[[i]]
    if (nrow(temp) > 30) {
      text = paste0("Showing top 30 Terms, for the complete list consult the ",
                    type, " in different TADs.csv")
      temp <- temp[1:30, ]
    } else {
      text = ""
    }
   
    p <- temp %>%
      ggplot(aes(x = as.factor(reorder(Term, -p_adjust)), y = p_adjust)) +
      geom_histogram(fill="#66ccff", color="#e9ecef", stat = "identity", width = 0.8) +
      labs(title = paste0("Adjusted P values of ", type, " in ", temp$TAD[1]), 
           subtitle = text) +
      ylab("Adjusted P value") +
      xlab(type) +
      theme_bw() + 
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 11),
            axis.title.x = element_text(size = 13, hjust = 0.5, vjust = 0.5,face = "bold"),
            axis.title.y = element_text(size = 13, hjust = 0.5, vjust = 0.5, face = "bold"),
            axis.text.x = element_text(size = 13, color = "black"),
            axis.text.y = element_text(size = 13, color = "black")
      ) +
      coord_flip() 
    
    save_image(print(p), 
               file.name =paste0(string, "/Adjusted P values per TAD/",
                                 "Adjusted P values of ",temp$TAD[1], ".png"), 
               height = (0.5 * nrow(temp) + 2), 
               width = 13, res = 100)
    
  }
  
  rm(data_plot, tad_grouping, new_grouping, temp, p)
  
}




# #This function is called from the "enrichrVisual.R" script
# #It creates a folder with graphs of KEGG pathways' maps using the Pathview tool
# #The folder name is "Pathview"
# pathVisual <- function(pathview.input , output_folder){
#   
#   new.dir <- paste0(output_folder,"/Pathview")
#   dir.create(new.dir, showWarnings = FALSE)
#   
#   pathview.input <- group_by(pathview.input,ID)
#   
#   pathview.input$diff <- as.numeric(pathview.input$diff)
#   
#   groups <- group_split(pathview.input)
#   
#   loops <- c(1:length(groups))
#   for (l in loops){
#     
#     path <- groups[[l]]
#     
#     path <- path %>% remove_rownames %>% column_to_rownames(var = "Genes")
#     
#     if (nrow(path)>1){
#       
#       dir <- paste0(new.dir,"/",path$ID[1])
#       dir.create(dir, showWarnings = FALSE)
#       
#       current.folder <- paste0(getwd(),"/",path$ID[1],".pathview.png")
#       new.folder <- paste0(dir,"/",path$ID[1],".pathview.png")
#       
#       p.input <- as.matrix(path[,3])
#       rownames(p.input) <- rownames(path)
#       
#       p <- pathview(gene.data  = p.input,
#                     pathway.id =  path$ID[1],
#                     species    = "hsa",
#                     gene.idtype = "SYMBOL",
#                     same.layer = F,   #two-layer graph (node colors and labels are added + official gene symbols)
#                     kegg.dir = dir ,
#                     kegg.native = T,
#                     na.col = "white",keys.align = "y") 
#       
#       
#       file.copy(current.folder,new.folder,overwrite = T)
#       
#       file.remove(current.folder)
#       
#     }
#   }
# }


#' This function is called by the `motif_visual` function.
#' 
#' @description
#' It creates a folder with graphs of the Motif EA for each TAD.
#' The folder is called "Plots per TADs".
#' For each TAD it creates a folder with two graphs.

#' The first one shows the summary of the Motif EA for this TAD, 
#' it is called PWMEnrich_Image.
#' It is created with the `PWMEnrich::plot()` function.

#' The second one is called Histogram Image.
#' It shows the TFs found in the TAD, arranged in ascending order,
#' according to the mean P value.

#' It also returns the data.table data_density, with the P values 
#' of the Motif EA of all TADs.
per_tad_plots <- function(report_list,
                          image_output_folder) {
  
  dir.create(paste0(image_output_folder, "/Plots per TADs"), showWarnings = FALSE)
 
  data_plot <- data.table(P.value = numeric(),
                          Adjusted.P.value = numeric())
  
  for (i in c(1:length(report_list))) {
    
    tad_dir <- paste0(image_output_folder, "/Plots per TADs/", report_list[[i]]$tad[1])
    dir.create(tad_dir,showWarnings = F)
    
    # plot 1
    png(filename = paste0(tad_dir, "/PWMEnrich_Image.png"),
        width = 1100, height = 1200)
    
    temp <- report_list[[i]]
    temp <- temp[, c(1:4, 7, 6)]
    
    if (nrow(temp) > 10)  temp <- temp[1:10, ]
    
    temp <- new("MotifEnrichmentReport", 
                d = temp, pwms = PWMLogn.hg19.MotifDb.Hsap@pwms[temp$id])
    
    p <- plot(temp, fontsize = 20, id.fontsize = 20)
    
    dev.off()

    tabl <- report_list[[i]]
    
    data_plot <- rbind(data_plot, data.table(P.value = tabl$p.value,
                                             Adjusted.P.value = tabl$adjusted.p.value))
    
    # Exclude uncharacterized motifs
    tabl <- tabl[str_detect(tabl$id, "UW.Motif.", negate = TRUE),]
    
    tabl <- tabl %>%
      dplyr::select(target, tad, adjusted.p.value) %>%
      group_by(target, tad) %>%
      summarise(target, tad, 
                adjusted.p.value = mean(adjusted.p.value)) %>%
      as.data.table() %>%
      unique()
  
    tabl <- setorder(tabl, adjusted.p.value)
    
    if (nrow(tabl) > 30){
      
      text <- "Showing top 30 Terms, for the complete list consult the over-represented TFs in each tad.csv"
      tabl <- tabl[1:30, ]
    
    } else {
      
      text <- "" 
    }
    
    # plot 2
    p <- tabl %>%
      ggplot(aes(x = as.factor(reorder(target, -adjusted.p.value)), 
                 y = adjusted.p.value)) +
      geom_histogram(fill="#66ccff", color="#e9ecef", stat = "identity", width = 0.8) +
      labs(title = paste0("Adjusted P values of TFs in ", tabl$tad[1]), 
           subtitle = text) +
      ylab("Mean Adjusted P value") +
      xlab("Transcription Factors") +
      theme_bw() + 
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size = 28, face = "bold"),
            plot.subtitle = element_text(size = 15),
            axis.title.x = element_text(size = 24,hjust = 0.5, vjust = 0.5, face = "bold"),
            axis.title.y = element_text(size = 24,hjust = 0.5, vjust = 0.5, face = "bold"),
            axis.text.x = element_text(size = 18, color = "black"),
            axis.text.y = element_text(size = 18, color = "black")
      ) +
      coord_flip() 
    
    save_image(print(p), 
               file.name = paste0(tad_dir, "/Histogram Image.png"),
               height = 0.5 * nrow(tabl) + 1,
               width = 15, res = 100)
  }
  
  rm(temp, p, tabl, tad_dir, text)
  
  return(data_plot)
  
}


#' This function is called by the `motif_visual` function.
#' 
#' @description 
#' It creates a folder with graphs of the Motif EA for each TF found.
#' The folder is called "Transcription Factors".
#' For each TF it creates a folder with two graphs.

#' The first one shows the binding motifs of the TF 
#' that were found in the analysis and is called Motifs.
#' It is created with the `ggseglogo::ggseqlogo()`` function.

#' The second one is called Barplot.
#' It shows the TADs that the TF was found into, arranged in
#' ascending order, according to the mean P value.

#' It also creates a histogram of the TFs found in most TADs, 
#' called "Top Transcription Factors".
#' The TFs are arranged in descending order, according to the 
#' number of TADs they were found into.

#' Lastly it creates a summary graph of the Top 10 TFs 
#' and their Top binding Motifs.
#' It is called Top 10 Transcription Factors.
per_tfs_plots <- function(report_list, 
                          data_visual, 
                          image_output_folder) {
  
  # Get PFMs for TFs
  PFMs <- list()
  names_list <- c()
  
  for (l in c(1:length(PWMLogn.hg19.MotifDb.Hsap@pwms))) {
    
    names_list <- c(names_list, PWMLogn.hg19.MotifDb.Hsap@pwms[[l]]@id)
    PFMs[[l]] <- PWMLogn.hg19.MotifDb.Hsap@pwms[[l]]@pfm 
    
  }
  names(PFMs) <- names_list
  
  rm(names_list)
  
  dir.create(paste0(image_output_folder, "/Transcription Factors"),
             showWarnings = FALSE)
  
  data_visual <- data_visual %>%
    separate_rows(adjusted_p_value, sep = "\\|")
  
  data_visual$adjusted_p_value <- as.numeric(data_visual$adjusted_p_value)
  data_visual <- data_visual %>%
    group_by(tfs, tad, top_motif)%>%
    summarise(tfs, tad, adjusted_p_value = mean(adjusted_p_value), motifs, top_motif) %>%
    unique()
  
  data_visual <- data_visual[which(data_visual$tfs != "NA"), ]
  data_plot <- group_by(data_visual, tfs)
  groups_tfs <- group_split(data_plot) 
 
  for (i in c(1:length(groups_tfs))) {
    
    temp <- groups_tfs[[i]]
    
    folder_name <- temp$tfs
    charac <- c(1:round(nchar(folder_name[1]) / 2))
    
    folder_name <- str_replace_all(folder_name,"/","-")
    folder_name <- str_replace_all(folder_name," / ","-")
    folder_name <- str_replace_all(folder_name,":","-") 
    
    string <- paste0(image_output_folder, "/Transcription Factors/", folder_name[1])
    dir.create(string, showWarnings = FALSE)
    
    # plot 1
    temp$adjusted_p_value <- as.numeric(temp$adjusted_p_value)
    
    temp <- temp %>%
      group_by(tad) %>%
      summarise(tfs, tad, 
                adjusted_p_value = mean(adjusted_p_value)) %>%
      unique()
    
    temp <- setorder(temp, adjusted_p_value)
    
    if (nrow(temp) > 29) {
      
      text <- "Showing top 30 Terms, for the complete list consult: TFs in different TADs.csv"
      temp <- temp[1:30, ]
    
    } else {
      
      text <- ""
    }

    p <- temp %>%
      ggplot(aes(x = as.factor(reorder(tad, -adjusted_p_value)),
                 y = adjusted_p_value)) +
      geom_histogram(fill="#66ccff", color="#e9ecef", stat = "identity", binwidth = 0.8) +
      labs(title = paste0("Adjusted P values of ", folder_name[1], " in different TADs"), 
           subtitle = text) +
      ylab("Mean Adjusted P value") +
      xlab("TAD number") +
      theme_bw() + 
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 13),
            axis.title.x = element_text(size = 15, hjust = 0.5, vjust = 0.5,face = "bold"),
            axis.title.y = element_text(size = 15, hjust = 0.5, vjust = 0.5, face = "bold"),
            axis.text.x = element_text(size = 14, color = "black"),
            axis.text.y = element_text(size = 14, color = "black")
      ) +
      coord_flip() 

    save_image(print(p), 
               file.name =paste0(string, "/Barplot.png"),
               height = 0.2 * nrow(temp) + 3,
               width = 13, res = 100)
    
    
    # plot 2
    temp <- groups_tfs[[i]]
    motifs <- c(str_split(temp$motifs[1], pattern = "\\|", simplify = TRUE))
    matrices_motifs <- list()
    
    loops <- c(1:length(motifs))
    for (l in loops) {
      
      matrices_motifs[[l]] <- PFMs[[motifs[l]]]
    }
    
    names(matrices_motifs) <- motifs
    
    p <- ggseqlogo(matrices_motifs, 
                   method = 'bits',
                   font = "roboto_slab_bold")
    
    save_image(print(p), 
               file.name = paste0(string, "/Motifs.png"),
               width = 14.5, height = 6.1, res = 100)
    
  }
  
  rm(p, data_plot, motifs, matrices_motifs, text, temp, groups_tfs)
  
  # plot all TFs
  
  data_visual$adjusted_p_value <- as.numeric(data_visual$adjusted_p_value)
  
  data_visual <- data_visual %>%
    group_by(tfs, tad, motifs, top_motif) %>%
    summarise(tfs, tad, motifs, 
              adjusted_p_value = mean(adjusted_p_value),
              top_motif) %>%
    as.data.table() %>%
    unique()
  
  terms <- dplyr::count(data_visual, tfs)
  
  data_plot <- merge(terms, data_visual) %>%
    dplyr::select(tfs, n, adjusted_p_value) %>%
    group_by(tfs, n) %>%
    summarise(tfs, n, adjusted_p_value = mean(adjusted_p_value)) %>%
    unique() 
  
  data_plot <- setorder(data_plot, -n)
  
  if (nrow(data_plot) > 29)  data_plot <- data_plot[1:30, ]

  colnames(data_plot) <- str_replace(colnames(data_plot), "n", "number_of_tads")
  
  p <- data_plot %>%
    ggplot(aes(x = as.factor(reorder(tfs,  number_of_tads)),
               y = adjusted_p_value, 
               fill = number_of_tads)) +
    geom_histogram(color = "#e9ecef", stat = "identity", width = 0.8) +
    labs(title = paste0("Top Transcription Factors in different TADs"), 
         subtitle = "Bar color corresponds to number of TADs each Term was found") +
    xlab("Transcription Factor") +
    ylab("Mean Adjusted P value") +
    coord_flip() +
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 24),
          plot.subtitle = element_text(size = 16),
          axis.text.x = element_text(size = 22, color = "black"),
          axis.text.y = element_text(size = 16, color = "black"),
          axis.title.x = element_text(size = 24, vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size = 24, vjust = 0.5, hjust = 0.5),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)
    ) +
    guides(fill = guide_legend(title = "number of TADs"))
  
  save_image(print(p), 
             file.name = paste0(image_output_folder, "/Top Transcription Factors.png"),
             height = 8, width = 14, res = 100)
  
  rm(data_plot, p)
  
  # plot TFs and motifs
  
  data_plot <- merge(terms, data_visual) %>%
    dplyr::select(tfs, n, adjusted_p_value, top_motif) %>%
    group_by(tfs, n, top_motif) %>%
    summarise(tfs, n, top_motif, 
              adjusted_p_value = mean(adjusted_p_value)) %>%
    unique() 
  
  data_plot <- setorder(data_plot, -n)
  
  data
  
  if (nrow(data_plot) > 9)  data_plot <- data_plot[1:10, ]
  
  colnames(data_plot) <- str_replace(colnames(data_plot), "n", "number_of_tads")
  
  data_plot <- setorder(data_plot, -adjusted_p_value)
  
  hist_plot <- data_plot %>%
    ggplot(aes(x = as.factor(reorder(tfs, -adjusted_p_value)),
               y = adjusted_p_value)) +
    geom_histogram(color = "#e9ecef",stat = "identity", width = 0.5) +
    labs(title = "Top 10 Transcription Factors in different TADs and their top motifs") +  
    xlab("Transcription Factor") +
    ylab("Mean Adjusted P value") +
    coord_flip() +
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 28),
          plot.subtitle = element_text(size = 18),
          axis.text.x = element_text(size = 26),
          axis.text.y = element_text(size = 28, color = "black"),
          axis.title.x = element_text(size = 26, vjust = 0.5,hjust = 0.5),
          axis.title.y = element_text(size = 26, vjust = 0.5,hjust = 0.5)
    ) 
 
  # plot8
 
  pfm_motifs <- list()
  
  data_plot <- setorder(data_plot, adjusted_p_value)

  for (l in c(1:10)) {
    pfm_motifs[[l]] <- PFMs[[data_plot$top_motif[l]]]
  }
  names(pfm_motifs) <- data_plot$top_motif
  
  logos <- ggseqlogo(pfm_motifs, method = 'bits', nrow = 10, ncol = 1) 
  
  empty <- plot.new()
  
  # labels_plot <- ggplot() +
  #   theme_void() +
  #   geom_text(aes(0, 0.1, label = data_plot$top_motif[10]), size = 5) +
  #   geom_text(aes(0, 0.2, label = data_plot$top_motif[9]), size = 5) +
  #   geom_text(aes(0, 0.3, label = data_plot$top_motif[8]), size = 5) +
  #   geom_text(aes(0, 0.4, label = data_plot$top_motif[7]), size = 5) +
  #   geom_text(aes(0, 0.5, label = data_plot$top_motif[6]), size = 5) +
  #   geom_text(aes(0, 0.6, label = data_plot$top_motif[5]), size = 5) +
  #   geom_text(aes(0, 0.7, label = data_plot$top_motif[4]), size = 5) +
  #   geom_text(aes(0, 0.8, label = data_plot$top_motif[3]), size = 5) +
  #   geom_text(aes(0, 0.9, label = data_plot$top_motif[2]), size = 5) +
  #   geom_text(aes(0, 1, label = data_plot$top_motif[1]), size = 5) +
  #   xlab(NULL)
  
  figure <- ggarrange(hist_plot,
                      ggarrange(empty, logos, empty, 
                                ncol = 1, nrow = 3, heights = c(0.4, 10, 0.5)),
                      ncol = 2, nrow = 1, 
                      widths = c(5, 2))
  
  save_image(print(figure), 
             file.name =paste0(image_output_folder, "/Top 10 Transcription Factors.png"),
             height = 15, width = 15, res = 100)
  
  rm(hist_plot, labels_plot, logos, empty, figure,
     data_plot, pfm_motifs, PFMs)
  
}


#' This function is called by the `motif_visual` function.
#' It creates a density plot of the P values of the Motif EA.
density_plot <- function(data_density, 
                         image_output_folder) {
  
  data_plot <- data.table(P.value = c(data_density$P.value, data_density$Adjusted.P.value),
                          Type = c(rep_len("P value", nrow(data_density)),
                                   rep_len("Adjusted P value", nrow(data_density))))
  
  p <- ggplot(data = data_plot, aes(x = P.value, group = Type, fill = Type)) +
    geom_density(adjust = 2, alpha = .8) +
    ggtitle("P values of Motifs") +
    xlab("P value") +
    ylab("Density") +
    theme_bw() + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 16),
          axis.title.x = element_text(size = 15, vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size = 15, vjust = 1.5, hjust = 0.5),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15)
    )
  
  save_image(print(p), 
             file.name =paste0(image_output_folder, "/Density plot-P values of Motifs.png"),
             res = 100)
  
  rm(data_plot, p)
}
