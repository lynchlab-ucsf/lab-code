# Script to determine OTUs that predict clinicial outcomes
## Katie McCauley
## Last Modified: August 25, 2017

#Read in OTU table:
workdir <- "/Volumes/data/Users/kmccauley/URECA_Nasal_Dust/Analysis/Subsets/" #Working directory
otutable <- "nasal_yr1.txt" #OTU table as a text file
mapping_file <- "Nasal_MapFile_yr1.txt" #Mapping file
response_var <- "site" #Variable you want to predict
classification <- TRUE # Do you want to use classification or regression

setwd(workdir)
if(readLines(otutable,n=1)=="# Constructed from biom file") {
  otu <- read.table(otutable, header=TRUE, sep="\t", skip=1, comment="", check.names=FALSE,row.names=1)
} else {
  otu <- read.table(otutable, header=TRUE, sep="\t", comment="", check.names=FALSE,row.names=1)
}

#Need to do some background work to clean up the taxonomies so they aren't long
otunames <- data.frame(OTUname=row.names(otu),tax=otu$taxonomy)
taxanames <- strsplit(as.character(otunames$tax),"; ")
cleantax <- lapply(taxanames, function(x) {
  x <- x[-7]
  x <- x[which(!x %in% c("g__","f__","o__","c__","p__","k__"))]
  x[x != "Unassigned"] <- substring(x[x != "Unassigned"], 4)
  x[length(x)]
})
genus_name <- sapply(cleantax, paste, collapse=" ")
otunames$highest_tax <- genus_name
otunames$otulabs <- gsub("_","~",otunames$OTUname)
otunames$fig_names <- paste0(otunames$highest_tax, " (", otunames$otulabs,")")
otu$taxonomy <- NULL

map <- read.table(mapping_file, header=TRUE, check.names=F, comment="", sep="\t",row.names=1)

predictors <- t(otu)
dim(predictors)

#Sorts the data so that the data.frame command later is ok
map <- map[rownames(predictors),]
response <- map[,response_var]

rf.df <- data.frame(response,predictors)
#Confirms categorical variable
if(classification == TRUE) {
  rf.df$response <- factor(rf.df$response)
} else {
  rf.df$response <- as.numeric(rf.df$response)
}
pacman::p_load(randomForest)
set.seed(123)
predoutcome <- randomForest(response ~ ., data=rf.df, ntree=5000, importance=TRUE)
print(predoutcome)

#Now I need to take the varImp data and make a pretty figure from it...
imp.df <- as.data.frame(importance(predoutcome))

imp.df <- merge(imp.df, otunames, by.x=0, by.y="OTUname")
imp.df <- imp.df[order(imp.df$MeanDecreaseAccuracy),]

first.30 <- imp.df[(nrow(imp.df)-29):nrow(imp.df),]
first.30$fig_names <- factor(first.30$fig_names,levels=first.30$fig_names)

labs1 <- sapply(strsplit(as.character(first.30$fig_names), " "),
               function(x) {
                 parse(text=paste0("italic('", x[1], "')~", x[2]))
               })
library(ggplot2)
a <- ggplot(first.30, aes(y=fig_names, x=MeanDecreaseAccuracy)) + theme_light() + geom_point() + ylab(" ") + ggtitle("Top 30 Predictive OTUs") + scale_y_discrete(labels=labs1)
quartz()
a