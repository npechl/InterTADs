library(data.table)
library(microbenchmark)


x = runif(2029)

res  = microbenchmark::microbenchmark(sort(x, method = "quick"), 
                                      source("Data_Integration.R"))
                                      
res2 = microbenchmark::microbenchmark(sort(x, method = "quick"), 
                                      source("TADiff.R"))                                     