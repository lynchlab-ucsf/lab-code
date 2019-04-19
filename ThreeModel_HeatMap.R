# Three-Model Heat-Map

```{r}
library(extrafont)

# Function to bring in and format all Diff-Abundance Data
make_plots <- function(csvfile,outcome_val) {
read.data <- read.csv(csvfile)
sig.data <- read.data[read.data$qval.best < 0.1,]
sig.data$outcome <- outcome_val
sig.data
}
tip.order <- read.table("/data/Users/kmccauley/PROSE_NEW/AnalysisData/Tip_Order.txt",sep="\t",header=F)

fin <- rbind(make_plots("/data/Users/kmccauley/PROSE_NEW/AnalysisData/exac_ThreeModel_subset.csv","Exac (Overall)"),
make_plots("/data/Users/kmccauley/PROSE_NEW/AnalysisData/RV_ThreeModel_subset.csv","RV (Overall)"),
make_plots("/data/Users/kmccauley/PROSE_NEW/AnalysisData/HRV_A_ThreeModel_subset.csv","HRV-A (Overall)"),
make_plots("/data/Users/kmccauley/PROSE_NEW/AnalysisData/HRV_B_ThreeModel_subset.csv","HRV-B (Overall)"),
make_plots("/data/Users/kmccauley/PROSE_NEW/AnalysisData/HRV_C_ThreeModel_subset.csv","HRV-C (Overall)"))

fin$outcome <- factor(fin$outcome, levels=c("Exac (Overall)","RV (Overall)","HRV-A (Overall)","HRV-B (Overall)","HRV-C (Overall)","Exac (Placebo)","RV (Placebo)","HRV-A (Placebo)","HRV-B (Placebo)","HRV-C (Placebo)","RV (ICS)","HRV-A (ICS)","HRV-B (ICS)","HRV-C (ICS)","Exac (Xolair)","RV (Xolair)","HRV-A (Xolair)","HRV-B (Xolair)","HRV-C (Xolair)"))
taxanames <- strsplit(as.character(fin$taxonomy),";")
mat <- t(sapply(taxanames,
function(x,m) c(x,rep(NA,m-length(x))),
max(rapply(taxanames,length))))
newnames <- as.data.frame(cbind(substr(mat[,1],4,35),substr(mat[,2:7],5,35)))
names(newnames) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
newnames$bactnames <- as.character(newnames$Genus)
newnames$bactnames[newnames$bactnames == "" | is.na(newnames$bactnames)] <- as.character(newnames$Family[newnames$bactnames == "" | is.na(newnames$bactnames)])
newnames$bactnames[newnames$bactnames == "" | is.na(newnames$bactnames)] <- as.character(newnames$Order[newnames$bactnames == "" | is.na(newnames$bactnames)])
bactnames <- newnames$bactnames
bactnames[bactnames == "Planococcaceae"] <- "Staphylococcaceae"
fin <- cbind(fin,bactnames)
fin$OTUname2 <- gsub("_","~",fin$OTUname)
fin$bactnames2 <- paste0(bactnames," (",fin$OTUname2,")")
fin$finalnames <- factor(fin$bactnames2,levels=names(sort(table(fin$bactnames2))))

fin$mean.diff.bin[fin$best.coef < 1] <- 1
fin$mean.diff.bin[fin$best.coef > 1] <- 0
fin$mean.diff.bin <- factor(fin$mean.diff.bin,labels=c("Enriched","Depleted"))
# Change the direction of the "weighted mean difference"
fin$wgt_mean_diff <- -fin$wgt_mean_diff

dim(fin)

fin$OTUname_sorted <- factor(fin$OTUname,levels=tip.order$V1)
fin <- fin[order(fin$OTUname_sorted),]

#Drop obs sig in fewer than 3 analyses
#obs.to.drop <- table(fin$OTUname)[table(fin$OTUname) < mean(table(fin$OTUname))]
fin1 <- fin[1:(nrow(fin)/2),]
fin2 <- fin[(nrow(fin)/2)+1:nrow(fin),]
fin1$finalnames <- factor(fin1$finalnames)
fin1$finalnames <- factor(fin1$finalnames,unique(as.character(fin1$finalnames)))
fin2$finalnames <- factor(fin2$finalnames)
fin2$finalnames <- factor(fin2$finalnames,unique(as.character(fin2$finalnames)))
```

### Heat Map
```{r fig.height=9,fig.width=4,dpi=500}
library(ggplot2)
library(plotrix)
#reorder_size <- function(x) {
#  factor(x, levels = tip.order$V1)
#}

labs1 <- sapply(strsplit(as.character(unique(fin$finalnames)), " "),
  function(x) {
    parse(text = paste0("italic('", x[1], "')~", x[2]))
})

labs2 <- sapply(strsplit(as.character(unique(fin2$finalnames)), " "),
  function(x) {
    parse(text = paste0("italic('", x[1], "')~", x[2]))
})

#The names weren't lining up, so I thought it might have something to do with the inherent underlying order of the factors
fin$finalnames <- factor(fin$finalnames,levels=unique(as.character(fin$finalnames)))
#png("HeatMapFig.png",width=500,height=2500)
p <- ggplot(fin, aes(outcome, finalnames,fill=mean.diff.bin)) + geom_tile() + theme(axis.text.x = element_text(angle = 60, hjust = 1), text=element_text(family="Avenir",size=5),legend.title=element_blank()) + xlab(" ") + ylab(" ") + scale_y_discrete(labels=labs1) + coord_fixed(ratio=1)
p
dev.off()
q <- ggplot(fin2, aes(outcome, finalnames,fill=mean.diff.bin)) + geom_tile() + theme(axis.text.x = element_text(angle = 60, hjust = 1),legend.title=element_blank()) + xlab(" ") + ylab(" ") #+ scale_y_discrete(labels=labs2)
q

```
### Make Table
```{r}
#keep only a certain set of variables
table <- subset(fin, select=c("OTUname","wgt_mean_diff","best.mod","best.pval","qval.best","mean.diff.bin","outcome","taxonomy"))
#Thinking that I should separate out into my four groups and then merge back together somehow... Maybe consider changing variable names instead of doing suffixes, though.
#Also need to figure out how to highlight cells a certain way, though
#Also, check the weighted mean difference value to make sure that it's right (or the enriched/depleted value)
exac <- subset(table, outcome=="Exac (Overall)")
rv <- subset(table, outcome=="RV (Overall)")
hrva <- subset(table, outcome=="HRV-A (Overall)")
hrvb <- subset(table, outcome=="HRV-B (Overall)")
hrvc <- subset(table, outcome=="HRV-C (Overall)")
make.table1 <- merge(exac,rv,all=TRUE,by="OTUname",suffixes=c(".a",".b"))
make.table2 <- merge(make.table1,hrva, all=TRUE,by="OTUname",suffixes=c(".c",".d"))
make.table3 <- merge(make.table2,hrvb, all=TRUE,by="OTUname",suffixes=c(".e",".f"))
make.table4 <- merge(make.table3,hrvc, all=TRUE,by="OTUname",suffixes=c(".g",".h"))
write.table(make.table4,"/data/Users/kmccauley/PROSE_NEW/PublicationTables/OTU_DE_Table.txt",sep="\t",quote=F,row.names=FALSE)
```

