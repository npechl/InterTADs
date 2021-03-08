#This file contains the functions used in the "Visualization" section of the "enrichmentAnalysis.R" script


#RStudio doesn't support all graph fonts in Windows OS
#This function sets the appropriate fonts for the graphs 
setGraphFonts <- function(deviceSystem){
  
  if (deviceSystem == "win"){
    #font_import() 
    loadfonts(device = "win")
    windowsFonts(Times=windowsFont("TT Times New Roman"))
    theme_set(theme_bw(base_size=12, base_family = 'Times New Roman')+ 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()))
  }
}


#This fuction is called from the "enrichrVisual" function
#It creates two density plots in the same graph
#The first density plot is made with the P values of the data after the analysis
#The second density plot is made with the Adjusted P values
compareDensityPlot <- function(string,type,data.visual){
  
  data.plot1 <- data.table(P.value = c(data.visual$P.value, data.visual$P.adjust),
                           Type = c(rep_len("P value",nrow(data.visual)),rep_len("Adjusted P value",nrow(data.visual))),
                           stringsAsFactors = F)
 
  p1 <- ggplot(data=data.plot1, aes(x=P.value, group=Type, fill=Type)) +
    geom_density(adjust=1.5, alpha=.4) +
    ggtitle(paste0("P values of " ,type)) +
    xlab("P value")+
    ylab("Density")+
    theme_ipsum()+
    theme(
      plot.title = element_text(size=28),
      axis.title.x = element_text(size = 24, vjust = 0.5,hjust = 0.5),
      axis.title.y = element_text(size = 24, vjust = 1.5,hjust = 0.5),
      legend.text = element_text(size=24),
      legend.title = element_text(size=24)
    )

  save_image(print(p1), file.name =paste(string, "/","Density plot-P values of ",type, ".png", sep =""))
  
}


#This function is called from the "enrichrVisual" function
#It creates two histogram graphs

#The first one shows the Terms that were found in the most TADs in descending order 
#according to the number of TADs they were found into.

#The second one shows the same Terms, but adds the information of the mean Adjusted P value 
#of each Term to be found in these TADs. 
#The Terms are arranged in ascending order, according to the mean Adjusted P value.
#The bar color corresponds to the number of TADs each Term was found into.
histogramPlot <- function(string,type,data.visual){
  
  terms <- dplyr::count(data.visual,Term)
  
  data.plot2 <- merge(terms, data.visual) %>%
    dplyr::select(Term, n, P.adjust) %>%
    group_by(Term,n) %>%
    summarise(Term, n, P.adjust = mean(P.adjust)) %>%
    unique() 
  
  data.plot2 <- setorder(data.plot2,-n) 
  
  if (nrow(data.plot2)>19){
    data.plot2 <- data.plot2[1:20,]
  }
  
 
  p2 <- data.plot2 %>%
    ggplot(aes(x = as.factor(reorder(Term, n)),y = n)) +
    geom_histogram(fill="#66ccff", color="#e9ecef", alpha=0.9,stat = "identity") +
    labs(title = paste0("Top ",type," in different TADs"))+
    xlab(type)+
    ylab("Number of TADs")+
    coord_flip() +
    theme_ipsum() +
    theme(
      plot.title = element_text(size=15),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_text(size = 13, vjust = 0.5,hjust = 0.5),
      axis.title.y = element_text(size = 13, vjust = 0.5,hjust = 0.5)
    )
 
  save_image(print(p2), file.name = paste(string, "/",type," in different TADs", ".png", sep = ""), height = (0.5*nrow(data.plot2)+2), width = 13)#, height = 7, width = 13)
  
  
  #plot 2.5
  colnames(data.plot2) <- str_replace(colnames(data.plot2),"n","number.of.TADs")
 
  p2.5 <- data.plot2 %>%
    ggplot(aes(x = as.factor(reorder(Term,  -P.adjust)),y = P.adjust, fill = number.of.TADs)) +
    geom_histogram(color = "#e9ecef",stat = "identity") +
    labs(title = paste0("Top ",type," in different TADs"), 
         subtitle = "Bar color corresponds to number of TADs each Term was found")+
    xlab(type)+
    ylab("Mean P.adjust")+
    coord_flip() +
    theme_ipsum() +
    theme(
      plot.title = element_text(size=28),
      plot.subtitle = element_text(size=22),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 22, colour = "black"),
      axis.title.x = element_text(size = 28, vjust = 0.5,hjust = 0.5),
      axis.title.y = element_text(size = 28, vjust = 0.5,hjust = 0.5),
      legend.title = element_text(size =24),
      legend.text = element_text(size =24)
    )
  
  save_image(print(p2.5), file.name = paste(string, "/Top ",type, ".png", sep = ""),height = (0.5*nrow(data.plot2))+2, width = 18)#, height = 7, width = 13)
  
}


#This function is called from the "enrichrVisual" function
#It creates a network graph of the Terms(max 10), that were found in the most TADs
#Each Term is a group
#Each node is unique and corresponds to a TAD
#The nodes are connected with lines 
#Each different line color corresponds to a group
networkPlot <- function(string,type,data.visual){
  
    
    #choose only 10 largest groups
    go_count <- dplyr::count(data.visual,Term)
    go_count <- setorder(go_count, -n)
    go_count <- go_count[1:10,]
    data.plot5 <- merge(data.visual,go_count, by = "Term")
    
    if (nrow(data.plot5)==10){return(NULL)}
    
    nodes <- data.plot5%>%
      dplyr::select(TAD,Term)
    
    gos <- nodes %>% 
      group_by(Term)
    groups <- group_split(gos)
    
    iterations <- c(1:length(groups))
    
    edges <- data.frame(from = character(), 
                        to = character() , 
                        group = character(), 
                        stringsAsFactors = FALSE)
    
    for (i in iterations){
      temp <- groups[[i]]
      iter <- c(1:nrow(temp))
      for (j in iter){
        from_to <- data.frame(from = temp$TAD[j], to =temp$TAD[1:nrow(temp)] , group = temp$Term[1]) 
        from_to <- from_to[-j,]
        edges <- dplyr::bind_rows(edges, from_to)
      }
    }
    
    edges <- unique(edges)
    if (nrow(edges) == 0) {return(NULL)}
    g <- graph_from_data_frame(edges, directed = FALSE)
    
   
    network <- ggraph(g) +
      geom_edge_link(aes(color = group), alpha = 1, linemitre = 10, n = 1000) +     
      geom_node_point(size = 5, shape = 21, stroke = 1,
                      fill = 'white', color = 'black') +
      geom_node_text(aes(label = name), repel = TRUE, size = 7 , fontface = "bold") +
      theme_void()+
      ggtitle(label =paste0("Top ",type," Network Graph"))+
      theme(
        plot.title = element_text(size=24, hjust = 0.5, vjust = 0.5, face = "bold") ,
        legend.text = element_text(size = 18),
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size = 18)
      )
    
    save_image(print(network), file.name =paste(string, "/Top ",type," network graph", ".png", sep = ""), height = 7, width = 17)

}


#This function is called by the "enrichrVisual" function
#It creates two folders with histogram graphs

#The first folder is called "Adjusted P values per TADs"
#Each histogram shows the Terms found in each TAD in ascending order, according to their Adjusted P value

#The second folder is called "Adjusted P values per Term"
#Each histogram shows the TADs, each Term was found into,
#in ascending order, according to their Adjusted P value
groupPlots <- function(string,type, data.visual){
  
  #group data according to Term 
  if (str_detect(type,"GO")){
    data.visual$ID <- paste0("GO ",data.visual$ID)
    data.plot4 <- data.visual
  }else if (str_detect(type,"KEGG")){
    data.plot4 <- data.visual
    data.plot4$ID <- data.plot4$Term
  }
  
  TAD_grouping <- data.plot4 %>%
    group_by(Term)
  new_grouping <- group_split(TAD_grouping)
  dir.create(paste(string, "/Adjusted P values per ",type, sep = ""), showWarnings = FALSE)
  iterations <- c(1:length(new_grouping))
  
  for (i in iterations) {
    
    temp <- new_grouping[[i]]
    if(str_detect(type,"GO", negate = T)){
      temp$Term <- rep_len("",nrow(temp))
    }
    temp$ID <- str_replace(temp$ID,"/","-")
    
    p4 <- temp %>%
      ggplot(aes(x = as.factor(reorder(TAD, -P.adjust)),y = P.adjust)) +
      geom_histogram(fill="#66ccff", color="#e9ecef", alpha=0.9, stat = "identity", binwidth = 0.5) +
      labs(title = paste0("Adjusted P values of ",temp$ID[1]," in different TADs"), 
           subtitle = temp$Term[1])+
      ylab("Adjusted P value")+
      xlab("TAD number")+
      theme_ipsum() +
      theme(
        plot.title = element_text(size=14, face = "bold"),
        plot.subtitle = element_text(size=11),
        axis.title.x = element_text(size = 13,hjust = 0.5, vjust = 0.5,face = "bold"),
        axis.title.y = element_text(size = 13,hjust = 0.5, vjust = 0.5, face = "bold")
        
      )+coord_flip() 

    save_image(print(p4), file.name =paste(string, "/Adjusted P values per ",type, "/Adjusted P values of ",temp$ID[1], ".png", sep = ""), height = (0.5*nrow(temp)+2), width = 13)
  }
  
  
  #group data according to TAD 
  TAD_grouping <- data.visual %>%
    group_by(TAD)
  new_grouping <- group_split(TAD_grouping)
  dir.create(paste(string, "/Adjusted P values per TAD", sep = ""), showWarnings = FALSE)
  iterations <- c(1:length(new_grouping))
  
  for (i in iterations) {
    
    temp <- new_grouping[[i]]
    if (nrow(temp)>30){
      text = paste0("Showing top 30 Terms, for the complete list consult the ",type, " in different TADs.csv")
      temp <- temp[1:30,]
    }else{
      text = ""
    }
   
    p5 <- temp %>%
      ggplot(aes(x = as.factor(reorder(Term, -P.adjust)),y = P.adjust)) +
      geom_histogram(fill="#66ccff", color="#e9ecef", alpha=0.9, stat = "identity") +
      labs(title = paste0("Adjusted P values of ",type," in ",temp$TAD[1]), 
           subtitle = text)+
      ylab("Adjusted P value")+
      xlab(type)+
      theme_ipsum() +
      theme(
        plot.title = element_text(size=14, face = "bold"),
        plot.subtitle = element_text(size=11),
        axis.title.x = element_text(size = 13,hjust = 0.5, vjust = 0.5,face = "bold"),
        axis.title.y = element_text(size = 13,hjust = 0.5, vjust = 0.5, face = "bold")
        
      )+coord_flip() 
    
    save_image(print(p5), file.name =paste(string, "/Adjusted P values per TAD/","Adjusted P values of ",temp$TAD[1], ".png", sep = ""), height = (0.5*nrow(temp)+2), width = 13)
    
  }
}


#This function is being called by the "enrichmentAnalysis.R" script
#It creates graphs of the EA using Enrichr
enrichrVisual <- function(folder, type, data.visual, criterio){
  
  #set or/and create image output folder
  if (str_detect(type,"GO")){
    string <- paste(folder,"/", type, sep = "")
    dir.create(string, showWarnings = FALSE) 
  }else{
   string <- folder
  }
  
  #P values and Adjusted P values Density Plots  
  compareDensityPlot(string,type,data.visual)
  
  #basic histogram Plots of the Terms found in most TADs
  histogramPlot(string,type,data.visual)

  networkPlot(string,type,data.visual) 
  
  groupPlots(string,type, data.visual)

}


#This function is called from the "enrichrVisual.R" script
#It creates a folder with graphs of KEGG pathways' maps using the Pathview tool
#The folder name is "Pathview"
pathVisual <- function(pathview.input , output_folder){
  
  new.dir <- paste0(output_folder,"/Pathview")
  dir.create(new.dir, showWarnings = FALSE)
  
  pathview.input <- group_by(pathview.input,ID)
  
  pathview.input$diff <- as.numeric(pathview.input$diff)
  
  groups <- group_split(pathview.input)
  
  loops <- c(1:length(groups))
  for (l in loops){
    
    path <- groups[[l]]
    
    path <- path %>% remove_rownames %>% column_to_rownames(var = "Genes")
    
    if (nrow(path)>1){
      
      dir <- paste0(new.dir,"/",path$ID[1])
      dir.create(dir, showWarnings = FALSE)
      
      current.folder <- paste0(getwd(),"/",path$ID[1],".pathview.png")
      new.folder <- paste0(dir,"/",path$ID[1],".pathview.png")
      
      p.input <- as.matrix(path[,3])
      rownames(p.input) <- rownames(path)
      
      p <- pathview(gene.data  = p.input,
                    pathway.id =  path$ID[1],
                    species    = "hsa",
                    gene.idtype = "SYMBOL",
                    same.layer = F,   #two-layer graph (node colors and labels are added + official gene symbols)
                    kegg.dir = dir ,
                    kegg.native = T,
                    na.col = "white",keys.align = "y") 
      
      
      file.copy(current.folder,new.folder,overwrite = T)
      
      file.remove(current.folder)
      
    }
  }
}


#This function is called by the "motifVisual" function
#It creates a folder with graphs of the Motif EA for each TAD
#The folder is called "Plots per TADs"
#For each TAD it creates a folder with two graphs

#The first one shows the summary of the Motif EA for this TAD, it is called PWMEnrich_Image
#It is created with the PWMEnrich::plot() function

#The second one is called Histogram Image
#It shows the TFs found in the TAD, arranged in ascending order, according to the mean P value

#It also returns the data.table data.plot5, with the P values of the Motif EA of all TADs
perTADPlots <- function(report.list,image_output_folder){
  
  #plot1+2
  dir.create(paste(image_output_folder, "/Plots per TADs", sep = ""), showWarnings = FALSE)
  iterations <- c(1:length(report.list))
  
  
  data.plot5 <- data.table(P.value = numeric(),
                           Adjusted.P.value = numeric())
  for (i in iterations){
    
    tad.dir = paste0(image_output_folder, "/Plots per TADs/",report.list[[i]]@d$tad[1])
    dir.create(tad.dir,showWarnings = F)
    
    #plot 1
    png(filename = paste(tad.dir, "/PWMEnrich_Image.png", sep = ""),width = 1100, height = 1200)
    
    temp <-report.list[[i]]
    temp@d <- temp@d[,c(1:4,7,6)]
    if (nrow(temp@d)>10){
      temp <- temp[1:10]
    }
    p <- plot(temp, fontsize=20, id.fontsize=20)
    
    dev.off()

    tabl <- report.list[[i]]@d
    
    tableFor5 <- tabl %>%
      dplyr::select(p.value,adjusted.p.value)
    colnames(tableFor5) <- c("P.value", "Adjusted.P.value")
    
    data.plot5 <- rbind(data.plot5,tableFor5)
    
    #exclude uncharacterized motifs
    tabl <- tabl[str_detect(tabl$id,"UW.Motif.",negate = TRUE),]
    tabl <- tabl %>%
      dplyr::select(target,tad,adjusted.p.value) %>%
      group_by(target,tad) %>%
      summarise(target,tad, adjusted.p.value = mean(adjusted.p.value))%>%
      as.data.table() %>%
      unique()
    #tabl = tabl %>%
     #   dplyr::select(target,tad,adjusted.p.value) %>%
      #tabl[,.(adjusted.p.value = mean(adjusted.p.value)), by = list(target),]
    tabl <- setorder(tabl, adjusted.p.value)
    
    if (nrow(tabl)>30){
      text = "Showing top 30 Terms, for the complete list consult the over-represented TFs in each tad.csv"
      tabl <- tabl[1:30,]
    }else{
      text = "" 
    }
    
    #plot 2
    p2 <- tabl %>%
      ggplot(aes(x = as.factor(reorder(target, -adjusted.p.value)),y = adjusted.p.value)) +
      geom_histogram(fill="#66ccff", color="#e9ecef", alpha=0.9, stat = "identity") +
      labs(title = paste0("Adjusted P values of TFs in ",tabl$tad[1]), 
           subtitle = text)+
      ylab("Mean Adjusted P value")+
      xlab("Transcription Factors")+
      theme_ipsum() +
      theme(
        plot.title = element_text(size=14, face = "bold"),
        plot.subtitle = element_text(size=11),
        axis.title.x = element_text(size = 11.5,hjust = 0.5, vjust = 0.5,face = "bold"),
        axis.title.y = element_text(size = 13,hjust = 0.5, vjust = 0.5, face = "bold")
        
      )+coord_flip() 
    
    save_image(print(p2), file.name =paste(tad.dir, "/Histogram Image.png", sep = ""), height = (0.5*nrow(tabl)+2), width = 13)
  }
  
  return(data.plot5)
  
}


#This function is called by the "motifVisual" function
#It creates a folder with graphs of the Motif EA for each TF found
#The folder is called "Transcription Factors"
#For each TF it creates a folder with two graphs

#The first one shows the binding motifs of the TF that were found in the analysis and is called Motifs
#It is created with the ggseglogo::ggseqlogo() function

#The second one is called Barplot
#It shows the TADs that the TF was found into, arranged in ascending order, according to the mean P value

#It also creates a histogram of the TFs found in most TADs, called "Top Transcription Factors"
#The TFs are arranged in descending order, according to the number of TADs they were found into

#Lastly it creates a summary graph of the Top 10 TFs and their Top binding Motifs
#It is called Top 10 Transcription Factors
perTFsPlots <- function(report.list, data.plot3.4, image_output_folder){
  
  #plot 3+4
  #get PWMs for TFs
  PFMs <- list()
  names.list <- c()
  
  for (l in c(1:length(PWMLogn.hg19.MotifDb.Hsap@pwms))){
    names.list <- c(names.list,PWMLogn.hg19.MotifDb.Hsap@pwms[[l]]@id)
    PFMs[[l]] <- PWMLogn.hg19.MotifDb.Hsap@pwms[[l]]@pfm 
    
  }
  names(PFMs) <- names.list
  
  dir.create(paste(image_output_folder, "/Transcription Factors", sep = ""), showWarnings = FALSE)
  data.plot3.4 <- data.plot3.4 %>%
    separate_rows(Adjusted.P.value, sep = "\\|")
  data.plot3.4$Adjusted.P.value <- as.numeric(data.plot3.4$Adjusted.P.value)
  data.plot3.4 <- data.plot3.4 %>%
    group_by(TFs,TAD,top.motif)%>%
    summarise(TFs,TAD, Adjusted.P.value = mean(Adjusted.P.value),Motifs, top.motif) %>%
    unique()
  data.plot3.4 <- data.plot3.4[which(data.plot3.4$TFs != "NA"),]
  data.plot3 <- data.plot3.4 %>%
    group_by(TFs)
  groups.TFs <- group_split(data.plot3) 
  iterations <- c(1:length(groups.TFs))
  
  for (i in iterations) {
    
    temp <- groups.TFs[[i]]
    temp$folder.name <- temp$TFs
    charac <- c(1:round(nchar(temp$folder.name[1])/2))
    for(c in charac){
      temp$folder.name <- str_replace(temp$folder.name,"/","-")
      temp$folder.name <- str_replace(temp$folder.name," / ","-")
      temp$folder.name <- str_replace(temp$folder.name,":","-") 
    }
    string <- paste(image_output_folder, "/Transcription Factors/",temp$folder.name[1], sep = "")
    dir.create(string, showWarnings = F)
    
    #plot3
    temp$Adjusted.P.value <- as.numeric(temp$Adjusted.P.value)
    temp <- setorder(temp, Adjusted.P.value)
    
    if (nrow(temp)>29){
      text = "Showing top 30 Terms, for the complete list consult: TFs in different TADs.csv"
      temp <- temp[1:30,]
    }else{
      text = ""
    }

    p3 <- temp %>%
      ggplot(aes(x = as.factor(reorder(TAD,-Adjusted.P.value)),y = Adjusted.P.value)) +
      geom_histogram(fill="#66ccff", color="#e9ecef", alpha=0.9, stat = "identity", binwidth = 0.5) +
      labs(title = paste0("Adjusted P values of ",temp$folder.name[1]," in different TADs"), 
           subtitle = text)+
      ylab("Mean Adjusted P value")+
      xlab("TAD number")+
      theme_ipsum() +
      theme(
        plot.title = element_text(size=16, face = "bold"),
        plot.subtitle = element_text(size=13),
        axis.title.x = element_text(size = 15,hjust = 0.5, vjust = 0.5,face = "bold"),
        axis.title.y = element_text(size = 15,hjust = 0.5, vjust = 0.5, face = "bold")
        
      )+coord_flip() 

    save_image(print(p3), file.name =paste(string, "/Barplot.png", sep = ""), height = (0.2*nrow(temp)+3), width = 13)
    
    
    #plot4
    motifs <- c(str_split(temp$Motifs[1], pattern = "\\|",simplify = T))
    matrices.motifs <- list()
    
    loops <- c(1:length(motifs))
    for (l in loops){
      matrices.motifs[[l]] <- PFMs[[motifs[l]]]
    }
    names(matrices.motifs) <- motifs
    
    p4 <- ggseqlogo(matrices.motifs, method = 'bits', font = "roboto_slab_bold")
    
    save_image(print(p4), file.name =paste(string, "/Motifs.png", sep = ""), width = 14.5, height = 6.1)
    
  }
  
  
  #plot6
  
  data.plot3.4$Adjusted.P.value <- as.numeric(data.plot3.4$Adjusted.P.value)
  data.plot3.4 <- data.plot3.4 %>%
    group_by(TFs,TAD,Motifs, top.motif) %>%
    summarise(TFs, TAD, Motifs, Adjusted.P.value = mean(Adjusted.P.value),top.motif) %>%
    as.data.table() %>%
    unique()
  
  terms <- dplyr::count(data.plot3.4,TFs)
  
  data.plot6 <- merge(terms, data.plot3.4) %>%
    dplyr::select(TFs, n, Adjusted.P.value) %>%
    group_by(TFs,n) %>%
    summarise(TFs, n, Adjusted.P.value = mean(Adjusted.P.value)) %>%
    unique() 
  
  data.plot6 <- setorder(data.plot6,-n)
  
  if (nrow(data.plot6)>29){
    data.plot6 <- data.plot6[1:30,]
  }
  
  
  colnames(data.plot6) <- str_replace(colnames(data.plot6),"n","number.of.TADs")
  
  p6.5 <- data.plot6 %>%
    ggplot(aes(x = as.factor(reorder(TFs,  number.of.TADs)),y = Adjusted.P.value, fill = number.of.TADs)) +
    geom_histogram(color = "#e9ecef",stat = "identity") +
    labs(title = paste0("Top Transcription Factors in different TADs"), 
         subtitle = "Bar color corresponds to number of TADs each Term was found")+
    xlab("Transcription Factor")+
    ylab("Mean Adjusted P value")+
    coord_flip() +
    theme_ipsum() +
    theme(
      plot.title = element_text(size=30),
      plot.subtitle = element_text(size=24),
      axis.text.x = element_text(size = 28),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 28, vjust = 0.5,hjust = 0.5),
      axis.title.y = element_text(size = 28, vjust = 0.5,hjust = 0.5),
      legend.title = element_text(size =20),
      legend.text = element_text(size =20)
    )
  
  save_image(print(p6.5), file.name =paste(image_output_folder, "/Top Transcription Factors.png", sep = ""), height = 8, width = 14)
  
  #plot7
  
  data.plot7 <- merge(terms, data.plot3.4) %>%
    dplyr::select(TFs, n,Adjusted.P.value,top.motif) %>%
    group_by(TFs,n, top.motif) %>%
    summarise(TFs, n,top.motif, Adjusted.P.value = mean(Adjusted.P.value)) %>%
    unique() 
  
  data.plot7 <- setorder(data.plot7,-n)
  
  if (nrow(data.plot7)>9){
    data.plot7 <- data.plot7[1:10,]
  }
  
  
  colnames(data.plot7) <- str_replace(colnames(data.plot7),"n","number.of.TADs")
  
  data.plot7 <- setorder(data.plot7, -Adjusted.P.value)
  p7 <- data.plot7 %>%
    ggplot(aes(x = as.factor(reorder(TFs, -Adjusted.P.value)),y = Adjusted.P.value, fill = number.of.TADs)) +
    geom_histogram(color = "#e9ecef",stat = "identity", width = 0.5) +
    labs(title = paste0("Top 10 Transcription Factors in different TADs and their top motifs"))+ #, 
    #    subtitle = "Bar color corresponds to number of TADs each Term was found")+
    xlab("Transcription Factor")+
    ylab("Mean Adjusted P value")+
    coord_flip() +
    theme_ipsum() +
    theme(
      plot.title = element_text(size=28),
      plot.subtitle = element_text(size=18),
      axis.text.x = element_text(size = 26),
      axis.text.y = element_text(size = 28, color = "black"),
      axis.title.x = element_text(size = 26, vjust = 0.5,hjust = 0.5),
      axis.title.y = element_text(size = 26, vjust = 0.5,hjust = 0.5),
      legend.title = element_text(size =20),
      legend.text = element_text(size =20)
    )
 
  #plot8
 
  pfm.motifs <- list()
  
  data.plot7 <-setorder(data.plot7,Adjusted.P.value)
  
  loops <- c(1:10)
  for (l in loops){
    pfm.motifs[[l]] <- PFMs[[data.plot7$top.motif[l]]]
  }
 
  
  p8 <- ggseqlogo(pfm.motifs, method = 'bits', nrow = 10, ncol = 1) 
  
  p9 <- plot.new()
  
  p10 <- ggplot() +
    theme_void() +
    geom_text(aes(0,0.1,label= data.plot7$top.motif[10]), size = 5) +
    geom_text(aes(0,0.2,label= data.plot7$top.motif[9]), size = 5) +
    geom_text(aes(0,0.3,label= data.plot7$top.motif[8]), size = 5) +
    geom_text(aes(0,0.4,label= data.plot7$top.motif[7]), size = 5) +
    geom_text(aes(0,0.5,label= data.plot7$top.motif[6]), size = 5) +
    geom_text(aes(0,0.6,label= data.plot7$top.motif[5]), size = 5) +
    geom_text(aes(0,0.7,label= data.plot7$top.motif[4]), size = 5) +
    geom_text(aes(0,0.8,label= data.plot7$top.motif[3]), size = 5) +
    geom_text(aes(0,0.9,label= data.plot7$top.motif[2]), size = 5) +
    geom_text(aes(0,1,label= data.plot7$top.motif[1]), size = 5) +
    xlab(NULL)
  
  
  
  figure <- ggarrange(p7, ggarrange(p9,p10,p9,ncol = 1, nrow = 3, heights = c(0.4,10, 0.5)),
                      ggarrange(p9,p8,p9,ncol = 1, nrow = 3, heights = c(0.4,10, 0.5)),
                      ncol = 3, nrow =1, widths = c(4,1.3,1.3), heights = c(1.2,0.5, 0.5) )
  
  save_image(print(figure), file.name =paste(image_output_folder, "/Top 10 Transcription Factors.png", sep = ""), height = 15, width = 15)
  
}


#This function is called by the "motifVisual" function
#It creates a density plot of the P values of the Motif EA
densityPlot <- function(data.plot5, image_output_folder){
  
  data.plot <- data.table(P.value = c(data.plot5$P.value, data.plot5$Adjusted.P.value),
                           Type = c(rep_len("P value",nrow(data.plot5)),rep_len("Adjusted P value",nrow(data.plot5))),
                           stringsAsFactors = F)
  
  p5 <- ggplot(data=data.plot, aes(x=P.value, group=Type, fill=Type)) +
    geom_density(adjust=2, alpha=.8) +
    ggtitle("P values of Motifs") +
    xlab("P value")+
    ylab("Density")+
    theme_ipsum()+
    theme(
      plot.title = element_text(size=16),
      axis.title.x = element_text(size = 15, vjust = 0.5,hjust = 0.5),
      axis.title.y = element_text(size = 15, vjust = 1.5,hjust = 0.5),
      legend.text = element_text(size=14),
      legend.title = element_text(size=15)
    )
  
  save_image(print(p5), file.name =paste(image_output_folder, "/Density plot-P values of Motifs.png", sep = ""))
  
  
}


#This function is being called by the "enrichmentAnalysis.R" script
#It creates graphs of the Motif EA using
motifVisual <- function(image_output_folder,motif_output_folder, data.plot3.4, report.list, criterio){
  
  data.plot5 <- perTADPlots(report.list,image_output_folder)
  
  perTFsPlots(report.list, data.plot3.4, image_output_folder )
  
  densityPlot(data.plot5, image_output_folder)
  
}

