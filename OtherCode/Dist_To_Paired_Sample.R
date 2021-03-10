rm(list=ls())

setwd("/data/Users/kmccauley/PROSE_NEW/AnalysisData")
wuf <- "BetaDiv_qiime/weighted_unifrac_dm.txt"
open.map <- read.table("FINAL_MAP_FILE.txt", header=T, sep="\t",comment="")

dm <- read.table(wuf,header=T,sep="\t",row.names=1)

betadiv.dist <- NULL
obs.count <- table(open.map$studyid)
tt <- nrow(obs.count)
for(i in 1:tt)  {
id.nm <- names(obs.count[i])
samp.names <- as.character(open.map$X.SampleID[open.map$studyid == id.nm])
first.samp <- as.character(samp.names[1])
#If somebody has MORE than one sample (which is expected), do the following
if(sum(samp.names != first.samp)>1) {
newdm <- dm[c(samp.names[1]),c(samp.names[2])]
sampid <- names(newdm)
beta.dist.1 <- as.numeric(newdm)
newdm <- cbind(id.nm,beta.dist.1)
}
#For when there is only one sample
else {
sampid <- samp.names
newdm <- cbind(id.nm,0)
}
betadiv.dist <- rbind(newdm,betadiv.dist)
}

