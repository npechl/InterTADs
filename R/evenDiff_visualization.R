# Loading libraries -----------------

library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel)


rm(list = ls())

# Inputs ----------------------

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

names.meta = c("ConsClust",
               "IC50beforeTreatment",
               "treatedAfter",
               "died",
               "IGHV", 
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

file_to_read = "TP53_evenDiff.txt"

data.all = fread(paste(output_folder, "/", file_to_read, sep = ""),
                 sep = "\t")

data.all$source = "methylation"
data.all[which(data.all$parent == expr_data), ]$source = "expression"

meta = fread(paste(dir_name, meta, sep = "/"))
who = meta == ""
who = apply(who, 1, sum, na.rm = TRUE)
meta = meta[which(who == 0), ]

# Complex heatmap ----------------------

column_ha = HeatmapAnnotation(Oakes = meta$ConsClust, 
                              TTT = anno_points(meta$T5), 
                              Mut_Status= meta$IGHV,
                              #IGHV= meta$`IGHV Uppsala gene usage`,
                              #SHM = anno_points(meta$`IGHV Uppsala % SHM`),
                              col = list(Oakes = c("HP" = "darkorange2", "IP" = "lightgreen", "LP" = "lightblue3"),
                                         Mut_Status = c("M" = "gray", "U" = "black")),
                              na_col = "white")


visualize_matrix = as.matrix(data.all[ ,meta$newNames, with = FALSE])

Heatmap(name = str_remove(file_to_read, ".txt"),
        visualize_matrix, 
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete", 
        
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        
        row_split = data.all$source,
        
        show_column_names = FALSE,
        
        top_annotation = column_ha, 
        column_km = 3)


rm(visualize_matrix, column_ha)


# Count TAD events ---------------------

nevents_tad = data.all %>% group_by(tad_name) %>% count()
# nevents_tad = nevents_tad[order(nevents_tad$n, decreasing = TRUE), ]
# nevents_tad = nevents_tad[which(nevents_tad$n > 1), ]

nevents_tad = merge(data.all, nevents_tad, by = "tad_name")
nevents_tad = nevents_tad[,c(1:16, ncol(nevents_tad)), with = FALSE]

nevents_tad = nevents_tad[order(nevents_tad$n, decreasing = TRUE), ]

nevents_tad_annot = nevents_tad[which(nevents_tad$adj.P.Val <= 0.001), ]
nevents_tad_annot$event_id = str_split(nevents_tad_annot$ID, ";", simplify = TRUE)[,2]


ggplot(nevents_tad, aes(x = tad_name, y = -log10(adj.P.Val))) +
  
  # Show all points
  geom_point(aes(color = tad_name), size = 1.3, alpha = 0.5) +
  
  
  geom_label_repel(data = nevents_tad_annot, 
                   aes(x = factor(tad_name, levels = unique(tad_name)), 
                       y = -log10(adj.P.Val), 
                       label = event_id,
                       fill = as.character(parent)),
                   min.segment.length = 0.3) +
  
  scale_color_manual(values = rep(c("red", "blue"), length(unique(nevents_tad$tad_name)))) +
  
  theme_bw() +
  theme( 
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )

# Hierarchical clustering ---------------- 

data.hist = t(data.all[, meta$newNames, with = FALSE])

fit = hclust(dist(data.hist, method = "euclidean"), 
             method = "complete")

cut = cutree(fit, k = 3)

cut.info = as.data.frame(cut)
cut.info$newNames = row.names(cut.info)

colnames(cut.info) = c("hclust_cut", "newNames")
meta.all = merge(meta, cut.info, by = "newNames")

# Statistical significance ------------------------

analysis_table = matrix(0, nrow = length(names.meta), ncol = 2)

analysis_table = as.data.table(analysis_table)

colnames(analysis_table) = c("Property", "P.value")

analysis_table$Property = as.character(analysis_table$Property)
analysis_table$P.value = as.numeric(analysis_table$P.value)

for (i in 1:length(names.meta)) {
  
  meta.test = meta.all %>% select(names.meta[i], hclust_cut)
  
  inter = names.meta[i]
  
  meta.test = as.data.frame(meta.test)
  meta.test = meta.test[!is.na(meta.test[[1]]),]
  
  colnames(meta.test) = c("factor", "cut")
  
  test_data = meta.test %>% 
    group_by(cut) %>% count(factor)
  
  test_data = test_data %>% tidyr::spread(factor, n)
  
  test_data[is.na(test_data)] = 0
  
  test_data = as.matrix(test_data[,2:ncol(test_data)])
  
  x = fisher.test(test_data, alternative = "two.sided")
  
  
  analysis_table[i,]$Property = inter
  analysis_table[i,]$P.value = x$p.value
  
}

write.table(analysis_table,
            file = paste(output_folder, "/", 
                         str_replace(file_to_read, ".txt", "_metadata_significance.txt"), 
                         sep = ""), 
            row.names = FALSE, 
            quote = FALSE, 
            sep = "\t")

gos = analysis_table[order(analysis_table$P.value, decreasing = TRUE), ]
gos$marker = factor(gos$Property, levels = gos$Property)

gos$col = "skip"
gos[which(gos$P.value <= 0.05), ]$col = "sign"

# Diverging Barcharts
gr = ggplot(gos, aes(x = marker, y = P.value , label = P.value, fill = col)) + 
  geom_bar(stat = 'identity', width = 0.4, position = "dodge") +
  scale_fill_manual(values = c("skip" = "gray75", "sign" = "gray30")) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_hline(yintercept = 0.049, 
             linetype = "dotted", 
             color = "black", size = 0.9) +  
  labs(y = "P value", x = "", title = "Meta data significance between groups") +
  coord_flip()
  


####1 vs others Fisher 2X2
# names<- names.meta[c(12:47)]

analysis_table = matrix(0, nrow = length(names.meta), ncol = 2)

analysis_table = as.data.table(analysis_table)

colnames(analysis_table) = c("Property", "P.value")

analysis_table$Property = as.character(analysis_table$Property)
analysis_table$P.value = as.numeric(analysis_table$P.value)

for (i in 1:length(names.meta)) {
  
  meta.test = meta.all %>% select(names.meta[i], hclust_cut)
  inter = names.meta[i]
  
  meta.test = as.data.frame(meta.test)
  meta.test = meta.test[!is.na(meta.test[[1]]),]
  colnames(meta.test) = c("factor", "cut")
  
  test_data = meta.test %>% 
    group_by(cut) %>% count(factor)
  
  test_data = test_data %>% tidyr::spread(factor, n)
  
  test_data[is.na(test_data)] = 0
  
  test_data = as.matrix(test_data[,2:ncol(test_data)])
  
  t.data1 = test_data[1,]
  t.data2 = test_data[c(2:3),]
  
  t.data2 = colSums(t.data2)
  
  # sum.with = sum(t.data2$with)
  # sum.non = sum(t.data2$non)
  # 
  # 
  # t.data2 <- data.frame(
  #   with = sum.with,
  #   non = sum.non
  # )
  
  t.data = rbind(t.data1, t.data2)
  
  
  x = fisher.test(t.data, alternative = "two.sided")
  
  
  analysis_table[i,1] = inter
  analysis_table[i,2] = x$p.value
}


gos = analysis_table[order(analysis_table$P.value, decreasing = TRUE), ]
gos$marker = factor(gos$Property, levels = gos$Property)

gos$col = "skip"
gos[which(gos$P.value <= 0.05), ]$col = "sign"

# Diverging Barcharts
gr = ggplot(gos, aes(x = marker, y = P.value , label = P.value, fill = col)) + 
  geom_bar(stat = 'identity', width = 0.4, position = "dodge") +
  scale_fill_manual(values = c("skip" = "gray75", "sign" = "gray30")) +
  theme_bw() + 
  theme(legend.position = "none") +
  geom_hline(yintercept = 0.049, 
             linetype = "dotted", 
             color = "black", size = 0.9) +  
  labs(y = "P value", x = "", title = "Meta data significance between 1st group and others") +
  coord_flip()

# Clean eniroment --------------------

rm(analysis_table, gos, 
   cut.info, data.all, data.hist, 
   fit, gr, meta.test, nevents_tad, 
   nevents_tad_annot, t.data, test_data, x)

# Survival time ----------------------

library("survminer")
library(survival)

fit = survfit(Surv(data.k$T6, data.k$died) ~ hclust_cut, data = data.k)


ggsurvplot(fit,
                data = data.k, 
                size = 1, 
                # palette = c("#1B9E77", "#D95F02", "#D95F02"),   # custom color palettes
                conf.int = FALSE,                     # Add confidence interval
                pval = TRUE,                         # Add p-value
                risk.table = TRUE,                   # Add risk table
                risk.table.col = "strata",           # Risk table color by groups
                # legend.labs = c("Cluster", "Female"),               # Change legend labels
                risk.table.height = 0.2,            # Useful to change when you have multiple groups
                ggtheme = theme_bw()                 # Change ggplot2 theme
)

# plot(fit.s, col= c("red","green","blue"), 
#      xlab="Survival", ylab="Cum survival", lwd=3)
# text(0.2,0.1,"p= 0.07")

# ------------------------------

# library(dendextend)
# library(dplyr)

# dend_obj = as.dendrogram(fit)
# 
# col_dend = color_branches(dend_obj, k = 3)
# plot(col_dend)
# 
# data.hist_cl = mutate(as.data.frame(data.hist), cluster = cut)

# count(as.numeric(data.hist_cl),cluster)

# -----------------------------

# colnames(meta.all)

# colnames(meta.all[50])<- paste("SHM")
# groups_data<-meta.all %>% 
#   group_by(cut) %>% 
#   summarise(count = n(), 
#             mean_T5 = mean(T5, na.rm = TRUE),
#             mean_SHM = mean(50, na.rm = TRUE),
#             Oakes_HP= length(which(ConsClust == "HP")),
#             Oakes_IP= length(which(ConsClust == "IP")),
#             Oakes_LP= length(which(ConsClust == "LP")))
# 
# meta.all$IC50beforeTreatment
# nrow(filter(meta.all, IGHV == "U"))
