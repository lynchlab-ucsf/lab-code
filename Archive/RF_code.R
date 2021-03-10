rm(list=ls())
#setwd("/data/Users/kmccauley/PROSE_DATA")
#source("PROSE_merge.R")
setwd("/data/Users/kmccauley/PROSE/")
mapfile <- "CLEAN_MAP_FILE.txt"
map <- read.table(mapfile,sep="\t", comment = "",header=T,check.names=F)
#Comment="" : ensures that the pound sign isn't read as a comment
#check.names=F : Keeps the hashtag in the SampleID name

#This is me cleaning up my OTU table, making the columns the OTU, and the rows the SampleID
otutable <- "/data/Users/kmccauley/PROSE/RawData/PROSE_OTU_Rarefy_2000_tax.txt"
raw.otutab <- read.table(otutable, sep="\t",header=T)
#Remove the first (OTUID) and last (taxonomy) column
mat.otutab <- data.matrix(raw.otutab[,2:(ncol(raw.otutab)-1)])
#Now rotate (transform function)
t.otutab <- t(mat.otutab)
#Make the OTUID variable the column names
colnames(t.otutab)<- as.character(raw.otutab$OTUID)
#Assign rownames to a "SampleID" variable (currently just row names, and I want to be able to merge on SampleID)
SampleID <-rownames(t.otutab)
#Somthing I was trying to do to make "taxonomy" shorter (don't run!)
#cut.tax <- strsplit(as.character(raw.otutab$taxonomy),split="; ")
#raw.otutab$family <- as.character(lapply(cut.tax,"[",5))
#raw.otutab$genus <- as.character(lapply(cut.tax,"[",6))
#names(t.otutab) <- paste0(raw.otutab$family,raw.otutab$genus)

#Finally, paste SampleID and the transformed OTU table togethermat.otutab <- data.matrix(raw.otutab[,2:(ncol(raw.otutab)-1)])
final.otutab <- data.frame(SampleID, t.otutab)



#Code to rename the variables to their taxa names, not OTU# (OK to ignore this section and go down to the next section
taxanames <- strsplit(as.character(raw.otutab$taxonomy),"; ")
mat <- t(sapply(taxanames,
function(x,m) c(x,rep(NA,m-length(x))),
max(rapply(taxanames,length))))
newnames <- substr(mat,4,35)
final.names <- paste(newnames[,4],newnames[,5],newnames[,6],newnames[,7])
final.names[final.names == "NA NA NA NA"] <- "Unassigned"
final.names <- c("SampleID",final.names)

names(final.otutab) <- final.names
otu.reftab <- raw.otutab[,c("OTUID","taxonomy")]


#source("PROSE_otu_data_sandbox.R")
#final data frame from the script is "dom.tax.final"
#map$cold_split[map$ncoldspostrandomization < 5] <- 0
#map$cold_split[map$ncoldspostrandomization >=5] <- 1
#I tried running a RF with "nonexac" as my variable to see if predictors of non-exacerbation were any different from predictors of exacerbation
#Got the same exact findings
df.for.rf <- merge(final.otutab,map[,c("RV.comb","#SampleID")],by.x=c("SampleID"),by.y=c("#SampleID"))
#Remove the sampleID (because you don't want your model to consider that variable
df.for.rf2 <- df.for.rf[c(-1)]
#Remove your response variable (here, RV.comb)
redo.rf.x <- df.for.rf2[,-ncol(df.for.rf2)]
#In order for RF to work, your variable needs to be numeric. RV.comb was character, so this is how I fix that problem.
redo.rf.y <- as.numeric(as.factor(df.for.rf2$RV.comb))
#Initialize randomForest
library(randomForest)

#rf.test <- randomForest(as.factor(RV.comb) ~ .,data=df.for.rf2,ntree=1000,na.action=na.omit,proximity=TRUE,importance=TRUE)

#An attempt in parallel (it kinda works!)
library('doMC')
#Tell R the number of cores to use
coreno <- 25
totreps <- 1000

registerDoMC(coreno)
#The first value is for the number of trees you want R to generate. So here, I am saying I want 1000 total trees, and 25 cores to be run, so 40 trees need to be run within each core to make it work out evenly (I now have variables for this; see above)
#In order to make my RF work, I needed to change the complexity parameter
rfinput <- foreach(ntree=rep(totreps/coreno,coreno),.combine=combine,.packages='randomForest') %dopar%
randomForest(redo.rf.x,redo.rf.y,ntree=ntree,importance=TRUE)

#rf.imp <- data.frame(importance(rfinput))
#rf.imp[order(rf.imp$X.IncMSE),]

setwd("/data/Users/kmccauley/PROSE/RandomForest/")
png("RANDOMFOREST_RV.png",width=1400,height=700)
varImpPlot(rfinput,sort=TRUE,cex=1,main="OTUs Most-Predictive of Rhinovirus")
dev.off()

########################################
### CART (not RF) for visualization ###
#######################################

#At some point, I was getting a "root only" error, but I was able to change the "cp" parameter
#I now get pretty complex trees, which is both exciting and bad (results may be due to error/random chance)
#"cp" gets smaller, more nodes in the tree
#(From R documentation)cp = "complexity parameter": any split that does not decrease the overall lack of fit by a factor of "cp" is not attempted. Helps prune off splits that are obviously not worthwhile
#Given that I needed to reduce the cp drastically, maybe the original "trunk only" model is a better indicator of predictiveness
setwd("/data/Users/kmccauley/PROSE/RandomForest/")
fit2 <- rpart(factor(RV.comb) ~ .,data=df.for.rf2,cp=0.003)
png("rpart_RV.png",width=800,height=800)
plot(fit2)
text(fit2,cex=1.4)
dev.off()
#Trying a slightly larger CP number (there must be some optimal cp here...)
#Interesting.... would optimal be 0.01/4 (groups)?
fit2 <- rpart(factor(RV.comb) ~ .,data=df.for.rf2,cp=0.0027)
pdf("rpart_RV_0027.pdf")
plot(fit2)
text(fit2,cex=0.5)
dev.off()



#previous RandomForest
rf.results <- readRDS("RandomForest_predALLOTU_outEXAC.rds")
rf.imp <- data.frame(importance(rf.results))
rf.imp[order(rf.imp$MeanDecreaseAccuracy),]
rf.imp$OTUname <- rownames(rf.imp)

library(rpart)
test.cart <- rpart(as.factor(ncoldspostrandomization) ~ ., data=df.for.rf2, method="class")
#to look up OTUs chosen:
otu.reftab[otu.reftab$OTUID == "_418","taxonomy"]

