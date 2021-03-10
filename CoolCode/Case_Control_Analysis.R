#Case-control data analysis

```{r}
rm(list=ls())
set.seed(123)
#Read in the mapping file
map <- read.table("/data/Users/kmccauley/PROSE_NEW/AnalysisData/FINAL_MAP_FILE.txt", header=TRUE, sep="\t",comment="",check.names=FALSE)

#Create a group of exacerbation samples
exac.samp.dat <- map[map$exac_sample == 1 & map$exac == 1,]
#There's somebody that has an exacerbation sample, but isn't considered to be exacerbation-prone. What does this mean?
#How many individuals and samples do I have? Are the numbers different?
nrow(exac.samp.dat)
exac.samp.dat$studyid <- factor(exac.samp.dat$studyid)
length(table(exac.samp.dat$studyid))

non.exac.samp.dat <- map[map$exac_sample == 0 & map$exac == 0,]

#Try Kei's suggestion of a for-loop
#get a list of exacerbation samples
exac.sampids <- exac.samp.dat[,1]
matched.data <- NULL
for(i in 1:length(exac.sampids)) {
dim(non.exac.samp.dat)
getid <- exac.sampids[i]
samp.to.match <- exac.samp.dat[exac.samp.dat[,1] == getid,]
crit1 <- samp.to.match[,"trt_group"]
crit2 <- samp.to.match[,"group"]
crit3 <- samp.to.match[,"viral_colldy"] + seq(-7:7)
potential.dat <- non.exac.samp.dat[non.exac.samp.dat$trt_group == crit1 & non.exac.samp.dat$group == crit2 & non.exac.samp.dat$viral_colldy %in% crit3,]
get.match <- potential.dat[sample(nrow(potential.dat),1),]
keep.match <- rbind(samp.to.match,get.match)
keep.match$matchid <- i
#drop the sample from the potential pool of selection
non.exac.samp.dat <- non.exac.samp.dat[non.exac.samp.dat$studyid != get.match$studyid,]
dim(non.exac.samp.dat)
matched.data <- rbind(keep.match,matched.data)
}

#test to see if it really worked!
table(matched.data$matchid)
```
