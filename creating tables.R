library(stringr)
library(data.table)
library(dplyr)
library(formattable)
library(tidyr)
library(png)

temp = read.csv("all_pre_data/all_pre_data_with_features.csv", sep = ";")
temp = as.data.table(temp)

meth_data = temp[which(str_detect(temp$ID, "cg")), ]

trans_data = temp[which(str_detect(temp$ID, "ENSG")), ]

var_data = temp[is.na(temp$ID), ]

rm(temp)

meth_data = meth_data[sample(nrow(meth_data), size = 5), ]
trans_data = trans_data[sample(nrow(trans_data), size = 3), ]
var_data = var_data[sample(nrow(var_data), size = 2), ]

temp = rbind(meth_data, trans_data, var_data)
# temp = temp[sample(nrow(temp)), ]

rm(meth_data, trans_data, var_data)


