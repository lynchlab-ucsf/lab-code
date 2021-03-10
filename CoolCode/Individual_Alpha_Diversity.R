#This code is intended to create plots of individual alpha diversity over time
## It can include additional data that may be relevant to your analysis
### Here, I included information on which samples were associated with an exacerbation and/or rhinovirus infection

# My alpha diversity variables are included in my mapping file

# The sample-specific dominant taxa is indicated by the color, and the size of the circle is indicative of the proportion of dominance (available from Kei's script)

#Set your working directory
setwd("/data/Users/kmccauley/PROSE_NEW/AnalysisData")
#Define the name of the mapping file (which already includes the alpha diversity scores)
mapfile <- "FINAL_MAP_FILE.txt"
#Bring in the mapping file
map <- read.table(mapfile,sep="\t", comment = "",header=T)

#I want to only include subjects with more than one sample, and I also want to have a list of unique ID numbers in your dataset
idnames <- table(map$studyid)
#Save study ids where there is more than one sample available
SID <- row.names(idnames[idnames >1])

#Here, I'm taking the dominant taxa names and making them look less "clunky"
taxanames <- strsplit(as.character(map$domgen_bysamp),";")
tnames <- t(matrix(unlist(taxanames),nrow=6))
tnames <- substr(tnames,4,30) #start at 4 because first 3 characters are of the format "k__"
map2 <- cbind(map,tnames)
library(RColorBrewer) #Rcolorbrewer allows me to "pick" random colors for each dominant taxa (ie, I have more than 8 dominant taxa, and R's color options typically max out at 8
#Making my own color palette
set.seed(128) #I used "randomness" to pick my color pallete, and changed the seed if I didn't like the color combination the random numbers generated
#In the list of colors, about 25% of them are various forms of grey, so if I don't remove the gray, then 25% of the colors, by radom chance, will be grey; I also don't want white
graylist <- grep("gray", colors(), value = TRUE)
greylist <- grep("grey", colors(),value=TRUE)
comblist <- c("white",graylist,greylist)
colorlist <- colors()
colorlist <- colorlist[!colorlist %in% comblist]
col <- sample(colorlist, 32, replace = FALSE) #sampling without replacement; 32 is the number of unique dominant taxa
#Initialize the color pallete we've created, named 'col'
palette(col)

#If you plan to use a text variable, make sure it's a character variable and not a factor variable
as.character(map2$hrv_type)->map2$hrv_type

map2$domtax.print <- factor(paste(map2[,"5"],map2[,"6"]))
levels(map2$domtax.print)[levels(map2$domtax.print) %in% c("er er"," ")] <- "Unknown"

setwd("/data/Users/kmccauley/PROSE_NEW/Figures")
map2 <- map2[order(map2$studyid, map2$viral_colldy),]
all <- unique(as.character(map2$studyid))
#placebo.exac <- unique(as.character(map2$studyid[map2$group == "Placebo" & map2$exac == 1]))
#placebo.nonexac <- unique(as.character(map2$studyid[map2$group == "Placebo" & map2$exac == 0]))
#xolair.exac <- unique(as.character(map2$studyid[map2$group == "Xolair" & map2$exac == 1]))
#xolair.nonexac <- unique(as.character(map2$studyid[map2$group == "Xolair" & map2$exac == 0]))
#ics.exac <- unique(as.character(map2$studyid[map2$group == "ICS" & map2$exac == 1]))
#ics.nonexac <- unique(as.character(map2$studyid[map2$group == "ICS" & map2$exac == 0]))

pdf("all_participants.pdf",paper="letter", pointsize=9,width=7.5,height=10)#previously the width was 880
par(mfrow=c(5,4),mar=c(1.1,1.1,1.1,1.1),xpd=TRUE) #previously the margins were c(5.1,4.1,4.1,20.1)
for(i in 1:length(all)) {
idno <- all[i]
ind.data <- map2[map2$studyid %in% idno,]
RV.plot.var <- as.factor(as.character(ind.data$hrv_type[!is.na(ind.data$hrv_type)]))
RV.plot.time <- ind.data$time.diff[!is.na(ind.data$viral_hrv_species)]
plot(ind.data$viral_colldy,ind.data$PD_whole_tree,col=ind.data$domtax.print,ylim=c(0,max(map$PD_whole_tree)),xlim=c(0,200),xlab="Days Since First Sample Collected",ylab="Faith's Diversity",cex=c(ind.data$domgen_bysamp_pct*5),pch=16)
text(ind.data$viral_colldy,ind.data$PD_whole_tree,pos=3,offset=1.5,labels=ind.data$hrv_type,cex=1.3)
points(ind.data$viral_colldy[ind.data$exacid>0],ind.data$PD_whole_tree[ind.data$exacid>0]-1.25,pch=17,cex=1.5)
#legend("topright",inset=c(-0.4,0),legend=as.factor(levels(map2$domtax.print)),col=1:length(levels(map2$domtax.print)),pch=16,cex=0.8)
legend("topleft",legend=unique(as.character(ind.data$domtax.print)),col=as.numeric(unique(ind.data$domtax.print)),pch=16,cex=0.9)
lines(ind.data$PD_whole_tree ~ ind.data$viral_colldy)
title(unique(ind.data$studyid))
}
dev.off()


