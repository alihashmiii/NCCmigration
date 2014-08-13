setwd("~/Dropbox/DPhil/DysonModel/all_vers2/integrateAndSwitch/")
library(convert)
library(betr)
ISdata = read.csv("All relative to Baseline group log10 results 071414 renamed NA 2 replicates.csv")
ISmatrix <- data.matrix(ISdata[,2:97])
rownames(ISmatrix) <- ISdata[,1]
samples <- rownames(ISmatrix)
timePoints <- c(-120,-120,-120,-120,-120,0,0,0,0,0,0,0,16,16,2,2,30,30,45,45,4,4,60,60,8,8,90,90,106,106,92,92,120,120,135,135,94,94,150,150,98,98,180,180)
ISnoVEGF = ISmatrix[13:28,];
ISyesVEGF = ISmatrix[29:44,];
ISset <- ExpressionSet(assayData=t(ISmatrix[13:44,]))
prob <- betr(eset=ISset, cond=as.factor(c(rep("noVEGF",2*8),rep("yesVEGF",2*8))), timepoint=timePoints[13:44], replicate=rep(c(1:2),16), alpha=0.05)
