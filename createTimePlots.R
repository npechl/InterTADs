library(data.table)
library(gtools)

# timePerChr = fread("timePerChr.csv")
rowsPerChr = fread("rowsPerChr.csv")

# basic = cbind(rowsPerChr[1:22, ]$numOfRows, timePerChr[1:22, ])
# basic = rbind(basic, cbind(rowsPerChr[264:265, ]$numOfRows, timePerChr[264:265, ]))
# basic$chromosome_name = paste("chr", basic$chromosome_name, sep = "")
# 
# plot(basic$elapsed, type = "b", col = "blue", xlab = "Chromosomes", ylab = "time (sec)", xaxt='n', lwd = 2)
# lines(basic$user.self, type = "b", col = "red", xaxt='n', ann=FALSE, lwd = 2)
# lines(basic$sys.self, type = "b", col = "black", xaxt='n', ann=FALSE, lwd = 2)
# 
# axis(1, at = 1:24, labels = basic$chromosome_name)
# axis(3, at = 1:24, labels = basic$V1)
# # text(seq(1, 24, by=1), par("usr")[3] - 0.2, labels = basic$V1, pos = 1, xpd = TRUE)
# grid(nx = 25)
# 
# legend(18, 1000, legend = c("elapsed", "user.self", "sys.self"), col = c("blue", "red", "black"), 
#        lwd = 2, cex=0.8, text.font = 2, text.width = 1, x.intersp = 0.25)

timePerChr = fread("times_for_getting_features.csv")
timePerChr = timePerChr[mixedorder(timePerChr$chromosome_name), ]
# rowsPerChr = fread("rowsPerChr.csv")

basic = cbind(rowsPerChr[1:22, ]$numOfRows, timePerChr[1:22, ])
basic = rbind(basic, cbind(rowsPerChr[264:265, ]$numOfRows, timePerChr[264:265, ]))
basic$chromosome_name = paste("chr", basic$chromosome_name, sep = "")

plot(basic$elapsed, type = "b", col = "blue", xlab = "Chromosomes", ylab = "time (sec)", xaxt='n', lwd = 2)
lines(basic$user.self, type = "b", col = "red", xaxt='n', ann=FALSE, lwd = 2)
lines(basic$sys.self, type = "b", col = "black", xaxt='n', ann=FALSE, lwd = 2)

axis(1, at = 1:24, labels = basic$chromosome_name)
axis(3, at = 1:24, labels = basic$V1)
# text(seq(1, 24, by=1), par("usr")[3] - 0.2, labels = basic$V1, pos = 1, xpd = TRUE)
grid(nx = 25)

legend(18, 1000, legend = c("elapsed", "user.self", "sys.self"), col = c("blue", "red", "black"),
       lwd = 2, cex=0.8, text.font = 2, text.width = 1, x.intersp = 0.25)