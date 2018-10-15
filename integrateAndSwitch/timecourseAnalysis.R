setwd("~/Dropbox/DPhil/DysonModel/all_vers2/integrateAndSwitch/")
library(timecourse)
ISdata = read.csv("All relative to Baseline group log10 results 071414 renamed NA 2 replicates.csv")
ISmatrix <- data.matrix(ISdata[,2:97])
rownames(ISmatrix) <- ISdata[,1]
samples <- rownames(ISmatrix)
timePoints <- c(-120,-120,-120,-120,-120,0,0,0,0,0,0,0,16,16,2,2,30,30,45,45,4,4,60,60,8,8,90,90,106,106,92,92,120,120,135,135,94,94,150,150,98,98,180,180)
ISnoVEGF = ISmatrix[13:28,];
ISyesVEGF = ISmatrix[29:44,];
tc1 = mb.long(t(ISnoVEGF), times=8, reps=rep(2,96), rep.grp=rep(c(1:2),8), time.grp=timePoints[13:28])