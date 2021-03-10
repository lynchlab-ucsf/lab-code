# Script to determine OTUs that predict clinicial outcomes
## Katie McCauley
## Last Modified: February 18, 2021
## Edited to utilize phyloseq objects

workdir<- "/data/Users/kmccauley/URECA_Nasal_Dust/" ## where do you want your output files to go
myphy <- readRDS("/data/Users/kmccauley/URECA_Nasal_Dust/KMcCauley02/anly/GeneticAnalyses/ureca_phyloseq_updt.rds") # where is your phyloseq object
response_var <- "pheno_y7" #Variable you want to predict
classification <- TRUE # Do you want to use classification or regression

setwd(workdir)


nonmiss <- subset_samples(myphy, !is.na(sample_data(myphy)[,response_var]))

predictors <- t(otu_table(nonmiss))
dim(predictors)

response <- data.frame(sample_data(myphy)[,response_var])

rf.df <- data.frame(response,predictors)

#Confirms categorical variable
if(classification == TRUE) {
  rf.df[,response_var] <- as.factor(rf.df[,response_var])
} else {
  rf.df[,response_var] <- as.numeric(rf.df[,response_var])
}

pacman::p_load(randomForest)
set.seed(123)
myformula <- paste0(response_var, " ~ .")
predoutcome <- randomForest(as.formula(myformula), data=rf.df, ntree=1000, importance=TRUE)
print(predoutcome)
plot(predoutcome)
#Now I need to take the varImp data and make a pretty figure from it...
imp.df <- as.data.frame(importance(predoutcome))

imp.df <- merge(imp.df, tax_table(myphy)@.Data, by=0)
imp.df <- imp.df[order(imp.df$MeanDecreaseAccuracy),]

first.30 <- imp.df[(nrow(imp.df)-29):nrow(imp.df),]
first.30$fig_names <- paste0(first.30$Genus, " (", first.30$Row.names, ")")
first.30$fig_names <- factor(first.30$fig_names,levels=first.30$fig_names)

labs1 <- sapply(strsplit(as.character(first.30$Genus), " "),
               function(x) {
                 parse(text=paste0("italic('", x[1], "')~", x[2]))
               })
library(ggplot2)
a <- ggplot(first.30, aes(y=fig_names, x=MeanDecreaseAccuracy)) + theme_light() + geom_point() + ylab(" ") + ggtitle("Top 30 Predictive OTUs") + scale_y_discrete(labels=labs1)
ggsave("RandomForestResults.pdf", plot=a, device=cairo_pdf)
