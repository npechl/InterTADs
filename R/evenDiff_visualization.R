########## Loading libraries ########## 

# source("R/libraries.R")

library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(ggplot2)


rm(list = ls())

########### Inputs ###########

#' Input parameters for TADiff part
#' 
#' @param dir_name Directory of input datasets containing feature counts and frequency tables
#' 
#' @param output_folder Folder name for printing output tables
#' 
#' @param meta meta-data file name used
#' 
#' @param names.meta meta data columns to process (names or indexes)
#' 
#' @param expr_data Parent index of expression data. If no expression is provided, place FALSE

dir_name = "Datasets_bloodcancer"

output_folder = "results_bloodcancer"

meta = "metaData_groups.csv"

names.meta = c("IGHV", 
               "gain2p25.3", 
               "del8p12", 
               "gain8q24",
               "del9p21.3",
               "del11q22.3",
               "trisomy12",
               "del13q14_any",
               "del13q14_bi",
               "del13q14_mono",
               "del14q24.3",
               "del15q15.1",
               "del17p13",
               "Chromothripsis",
               "BRAF",
               "KRAS",
               "MYD88",
               "NOTCH1",
               "SF3B1",
               "TP53",
               "ACTN2",
               "ATM",
               "BIRC3",
               "CPS1",
               "EGR2",
               "FLRT2",
               "IRF2BP2",
               "KLHL6",
               "LRP1",
               "MED12",
               "MGA",
               "MUC16",
               "NFKBIE",
               "PCLO",
               "UMODL1",
               "XPO1",	
               "ZC3H18")

expr_data = 1

data.all = fread("results_bloodcancer/TP53_evenDiff.txt",
                 sep = "\t")

data.all$source = "methylation"
data.all[which(data.all$parent == expr_data), ]$source = "expression"

meta = fread(paste(dir_name, meta, sep = "/"))
who = meta == ""
who = apply(who, 1, sum, na.rm = TRUE)
meta = meta[which(who == 0), ]

########### Inputs ########## 

column_ha = HeatmapAnnotation(Oakes = meta$ConsClust, 
                              TTT = anno_points(meta$T5), 
                              Mut_Status= meta$IGHV,
                              #IGHV= meta$`IGHV Uppsala gene usage`,
                              #SHM = anno_points(meta$`IGHV Uppsala % SHM`),
                              col = list(Oakes = c("HP" = "darkorange2", "IP" = "lightgreen", "LP" = "lightblue3"),
                                         Mut_Status = c("M" = "gray", "U" = "black")),
                              na_col = "white")


visualize_matrix = as.matrix(data.all[ ,meta$newNames, with = FALSE])

Heatmap(name = "TP53_evenDiff",
        visualize_matrix, 
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete", 
        
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        
        row_split = data.all$source,
        
        top_annotation = column_ha, 
        column_km = 3)


rm(visualize_matrix, column_ha)


# count TAD events -------------

nevents_tad = data.all %>% group_by(tad_name) %>% count()
# nevents_tad = nevents_tad[order(nevents_tad$n, decreasing = TRUE), ]
# nevents_tad = nevents_tad[which(nevents_tad$n > 1), ]

nevents_tad = merge(data.all, nevents_tad, by = "tad_name")
nevents_tad = nevents_tad[,c(1:16, ncol(nevents_tad)), with = FALSE]

nevents_tad = nevents_tad[order(nevents_tad$n, decreasing = TRUE), ]

nevents_tad_annot = nevents_tad[which(nevents_tad$adj.P.Val <= 0.001), ]
nevents_tad_annot$event_id = str_split(nevents_tad_annot$ID, ";", simplify = TRUE)[,2]

library(ggrepel)

ggplot(nevents_tad, aes(x = tad_name, y = -log10(adj.P.Val))) +
  
  # Show all points
  geom_point(aes(color = tad_name), size = 1.3, alpha = 0.5) +
  
  
  geom_label_repel(data = nevents_tad_annot, 
                  aes(x = factor(tad_name, levels = unique(tad_name)), 
                      y = -log10(adj.P.Val), 
                      label = event_id,
                      fill = as.character(parent))) +
  
  scale_color_manual(values = rep(c("red", "blue"), length(unique(nevents_tad$tad_name)))) +
  
  # custom X axis:
  # scale_x_continuous( label = axisdf$tad_name, breaks= axisdf$center ) +
  # scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )

################################

data.hist = t(data.all[, meta$newNames, with = FALSE])

fit = hclust(dist(data.hist, method = "euclidean"), 
             method = "complete")

plot(fit)
cut = cutree(fit, k = 3)

cut.info = as.data.frame(cut)
cut.info$newNames = row.names(cut.info)
meta.all = merge(meta, cut.info, by.x = "newNames", by.y = "patient_id")

#write.table(meta.all, paste("metaData_groups.csv"), 
#row.names = FALSE, sep = "\t", quote = FALSE)

################################################

library(dendextend)
library(dplyr)

dend_obj = as.dendrogram(fit)

col_dend = color_branches(dend_obj, k = 3)
plot(col_dend)

data.hist_cl <- mutate(as.data.frame(data.hist), cluster = cut)
#count(as.numeric(data.hist_cl),cluster)


##################################################
####statistics 3x3

# names<- names.meta[c(12:47)]

analysis_table = matrix(0, nrow = length(meta$newNames), ncol = 2)

for (i in 1:length(meta$newNames)) {
  
  meta.test = meta.all %>% select(meta$newNames[i], cut)
  
  inter = meta$newNames[i]
  
  meta.test = as.data.frame(meta.test)
  meta.test = meta.test[!is.na(meta.test[,1]),]
  
  colnames(meta.test) = c("factor", "cut")
  
  
  test_data = meta.test %>% 
    
    group_by(cut) %>% 
    summarise(count = n(),
              with = length(which(factor == inter))) %>% 
    mutate(non = count - with)
  
  test_data = as.matrix(test_data[,c(3,4)])
  
  x = fisher.test(test_data, alternative = "two.sided")
  
  
  analysis_table[i, 1] = paste(names[i])
  analysis_table[i, 2] = paste(x$p.value)
  
}

#M vs U
meta.test<- meta.all %>% select(IGHV, cut)
inter<- "U"
colnames(meta.test)
meta.test<- as.data.frame(meta.test)
meta.test<-meta.test[!is.na(meta.test[,1]),]
colnames(meta.test)<- c("factor","cut")


test_data<-meta.test %>% 
  group_by(cut) %>% 
  summarise(count = n(),
            with = length(which(factor == inter))) %>% 
  mutate(non=count-with)

test_data<- as.matrix(test_data[,c(3,4)])

x2<-fisher.test(test_data,
                alternative="two.sided")


#Int VS other
colnames(meta.all)
meta.test<- meta.all %>% select(ConsClust, cut)
inter<- "IP"
colnames(meta.test)
meta.test<- as.data.frame(meta.test)
meta.test<-meta.test[!is.na(meta.test[,1]),]
colnames(meta.test)<- c("factor","cut")


test_data<-meta.test %>% 
  group_by(cut) %>% 
  summarise(count = n(),
            with = length(which(factor == inter))) %>% 
  mutate(non=count-with)

test_data<- as.matrix(test_data[,c(3,4)])

x1<-fisher.test(test_data,
                alternative="two.sided")


x.all<- data.frame(c("IGHV","Oakes (IP vs others)"), c(x2$p.value,x1$p.value))
analysis_table<- as.data.frame(analysis_table)
names(x.all)<- colnames(analysis_table)
analysis_table_all<- rbind(x.all,analysis_table)

#plot p value

write.table(analysis_table_all,file="analysis_table_allcomparison_accordingmetadata.txt", col.names=T, row.names=F, quote=F, sep="\t")

analysis_table_all<-as.data.frame(analysis_table_all)
gos <- analysis_table_all[order(-as.numeric(analysis_table_all$V2)), ]  # sort
gos$marker <- factor(gos$V1, levels=gos$V1)
head(gos)
# Diverging Barcharts
ggplot(gos, aes(x=marker, y=as.numeric(V2) , label=as.numeric(V2))) + 
  geom_bar(stat='identity', width=.4,position="dodge") +
  theme_bw() + geom_hline(yintercept = 0.049, linetype="dotted", 
                          color = "black", size=0.9) +  coord_flip()


####1 vs others Fisher 2X2
names<- names.meta[c(12:47)]

analysis_table<-matrix(0, nrow = 36, ncol = 2)

for (i in 1: length(names))
{
  meta.test<- meta.all %>% select(names[i], cut)
  inter<- names[i]
  colnames(meta.test)
  meta.test<- as.data.frame(meta.test)
  meta.test<-meta.test[!is.na(meta.test[,1]),]
  colnames(meta.test)<- c("factor","cut")
  
  
  test_data<-meta.test %>% 
    group_by(cut) %>% 
    summarise(count = n(),
              with = length(which(factor == inter))) %>% 
    mutate(non=count-with)
  
  test_data<- as.matrix(test_data[,c(3,4)])
  t.data1<- test_data[1,]
  t.data2<- as.data.frame(test_data[c(2,3),])
  sum.with<-sum(t.data2$with)
  sum.non<-sum(t.data2$non)
  t.data2 <- data.frame(
    with = sum.with,
    non = sum.non
  )
  
  t.data<- rbind(t.data1, t.data2)
  
  
  x<-fisher.test(t.data,
                 alternative="two.sided")
  
  
  analysis_table[i,1]<-paste(names[i])
  analysis_table[i,2]<-paste(x$p.value)
}

#M vs U
meta.test<- meta.all %>% select(IGHV, cut)
inter<- "U"
colnames(meta.test)
meta.test<- as.data.frame(meta.test)
meta.test<-meta.test[!is.na(meta.test[,1]),]
colnames(meta.test)<- c("factor","cut")


test_data<-meta.test %>% 
  group_by(cut) %>% 
  summarise(count = n(),
            with = length(which(factor == inter))) %>% 
  mutate(non=count-with)

test_data<- as.matrix(test_data[,c(3,4)])

x<-fisher.test(test_data,
               alternative="two.sided")




##############################
colnames(meta.all)
colnames(meta.all[50])<- paste("SHM")
groups_data<-meta.all %>% 
  group_by(cut) %>% 
  summarise(count = n(), 
            mean_T5 = mean(T5, na.rm = TRUE),
            mean_SHM = mean(50, na.rm = TRUE),
            Oakes_HP= length(which(ConsClust == "HP")),
            Oakes_IP= length(which(ConsClust == "IP")),
            Oakes_LP= length(which(ConsClust == "LP")))

meta.all$IC50beforeTreatment
nrow(filter(meta.all, IGHV == "U"))
#########################################



library("survival")
data.k<- meta.all

data.k$SurvObj <- with(data.k, Surv(T6))
fit.s<-survfit(SurvObj ~ cut, data = data.k)

survdiff(SurvObj ~ cut, data.k,  rho=0)
plot(fit.s, col= c("red","green","blue"), 
     xlab="Survival", ylab="Cum survival", lwd=3)
text(0.2,0.1,"p= 0.07")
