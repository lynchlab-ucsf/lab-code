#In R: create sampleID lists for each participant
setwd("/data/Users/kmccauley/PROSE_NEW/AnalysisData")
mapfile <- "/data/Users/kmccauley/PROSE_NEW/AnalysisData/FINAL_MAP_FILE.txt"
rawdata <- read.table(mapfile, header=T, comment="",sep="\t",check.names=F)
#get list of all unique ids (generate a frequency table of IDs and pull the row names)
unique.idlist <- row.names(table(rawdata$studyid))
adonis.output <- NULL
for(i in 1:length(unique.idlist))   {
   id <- unique.idlist[i]
#extract Sample IDs for the selected ID from the loop (and make data frame)
   SID <- as.data.frame(rawdata[,"#SampleID"][rawdata$studyid %in% id])
   ind.data <- as.data.frame(rawdata[rawdata$studyid %in% id,])
#   names(ind.data)[names(ind.data) == "viral_colldy"] <- "viral_colldy"
#   names(ind.data)[names(ind.data) == "time.count"] <- "time_count"
#   dir.create(paste0("/data/Users/kmccauley/PROSE_NEW/AnalysisData/IndPCOA/",id,"/"),recursive=TRUE)
   write.table(SID,paste0("/data/Users/kmccauley/PROSE_NEW/AnalysisData/IndPCOA/",id,"/samplist.txt"),col.names=F,quote=F,row.names=F)
   write.table(ind.data, paste0("/data/Users/kmccauley/PROSE_NEW/AnalysisData/IndPCOA/",id,"/map_subset.txt"),row.names=F,quote=F,sep="\t")
}

#extract samples by participant and generate weighted unifrac PCs
cd /data/Users/kmccauley/PROSE_NEW/AnalysisData/IndPCOA
for f in ./*; do
   if [ -d "$f" ]; then
      echo $f
      cd "$f"
      filter_samples_from_otu_table.py -i ../../PROSE_OTUtable_June.biom -o otu_subset.biom -m ../../FINAL_MAP_FILE.txt --sample_id_fp samplist.txt
      parallel_beta_diversity.py -i otu_subset.biom -o indPCOA -O 10 -m weighted_unifrac -t ../../../../PROSE_OTUPICK/filtered_otutable_nochimera_core.tre
      principal_coordinates.py -i indPCOA/weighted_unifrac_otu_subset.txt -o indPCOA/weighted_unifrac_pc.txt
      cd ..
   fi
done


library(scatterplot3d)
#Read in the main map file
map <- read.table("/data/Users/kmccauley/PROSE_NEW/AnalysisData/FINAL_MAP_FILE.txt",sep="\t",header=T,comment="")

#get a vector with the dominant taxa
taxanames <- strsplit(as.character(map$domgen_bysamp),";")
#put the list into a matrix and transform it
tnames <- t(matrix(unlist(taxanames),nrow=6))
#get rid of the "x__" at the beginning of the names
tnames <- substr(tnames,4,30)
#attach these 6 new vectors to the map file
map2 <- cbind(map,tnames)
#Create a "printable" vector (ie, one that looks good, and not like code)
map2$domtax.print <- factor(paste(map2[,"5"],map2[,"6"]))
#We have a couple of blank names (ie only know something higher than Family, or Unassigned)
levels(map2$domtax.print)[levels(map2$domtax.print) %in% c("er er"," ")] <- "Unknown"
#Take out people with only one sample
tab.ids <- table(map2$studyid)
ids <- row.names(tab.ids[tab.ids>1])
map2 <- map2[map2$studyid %in% ids,]

#format RV variable (currently viral_hrv_species)
#as.character(map2$viral_hrv_species)->map2$hrv.clean
#map2$hrv.clean <- gsub("ENTEROVIRUS", "E", map2$hrv.clean)
#map2$hrv.clean[is.na(map2$hrv.clean)] <- " "
#map2$hrv.clean <- as.factor(map2$hrv.clean)

library(RColorBrewer)
#Creating a set of 32 colors by random chance
#First, set the seed number
seedno <- 128
set.seed(seedno)
#Remove obviously bad colors (ie, greys and white)
graylist <- grep("gray", colors(), value = TRUE)
greylist <- grep("grey", colors(),value=TRUE)
comblist <- c("white",graylist,greylist)
colorlist <- colors()
colorlist <- colorlist[!colorlist %in% comblist]
#Randomly select from new list
col <- sample(colorlist, 32, replace = FALSE)
palette(col)

#Create vectors of study ids belonging to each group of interest; I wanted to create stratified documents
all <- unique(as.character(map2$studyid))
#placebo.exac <- unique(as.character(map2$studyid[map2$group == "Placebo" & map2$exac == 1]))
#placebo.nonexac <- unique(as.character(map2$studyid[map2$group == "Placebo" & map2$exac == 0]))
#xolair.exac <- unique(as.character(map2$studyid[map2$group == "Xolair" & map2$exac == 1]))
#xolair.nonexac <- unique(as.character(map2$studyid[map2$group == "Xolair" & map2$exac == 0]))
#ics.exac <- unique(as.character(map2$studyid[map2$group == "ICS" & map2$exac == 1]))
#ics.nonexac <- unique(as.character(map2$studyid[map2$group == "ICS" & map2$exac == 0]))

#Make into a list for the loop
#totstrat <- list(placebo.exac,placebo.nonexac,xolair.exac,xolair.nonexac,ics.exac,ics.nonexac)
#Give names to the sections to create the document names
#names(totstrat) <- c("placebo.exac","placebo.nonexac","xolair.exac","xolair.nonexac","ics.exac","ics.nonexac")

#This is the loop. Go through the entire thing six times
#for(j in 1:length(all))   {
#Set the working directory
setwd("/data/Users/kmccauley/PROSE_NEW/Figures")
#Initialize the PDF document, naming it based on the name of the IDs. Give a 1" margin
pdf("IndBetaDiv_ALL.pdf",paper="letter",pointsize=9,width=7.5,height=10)
par(mfrow=c(4,3),mar=c(1.1,1.1,1.1,1.1))
#unique.idlist <- totstrat[[j]]
for(i in 1:length(all))   {
   #assign the current working ID to a variable to use later
   id <- all[i]
#extract Sample IDs for the selected ID from the loop (and make data frame)
   setwd(paste0("/data/Users/kmccauley/PROSE_NEW/AnalysisData/IndPCOA/",id))
   indmap <- map2[map2$studyid %in% id,]
   #I need to know how many samples each ID has so that I can restrict the PC table (see nrows=samptot)
   samptot <- nrow(indmap)
   indpc <- read.table("indPCOA/weighted_unifrac_pc.txt",sep="\t",skip=9,nrows=samptot)
   indMDS <- merge(indmap,indpc,by.x=c("X.SampleID"),by.y=c("V1"))
   #I need a varable that removes the missing RV.comb factor values so that it doesn't confuse the scatterplot text
   RV.plot.var <- as.factor(as.character(indMDS$hrv_type[!is.na(indMDS$hrv_type)]))
   #Begin the 3D scatterplot code. Many of the options just make the plot look better (angles, no labels,etc.)
   s3d <- scatterplot3d(indMDS[,c("viral_colldy","V2","V3")],cex.symbols=indMDS$domfam_bysamp_pct*5,angle=70,grid=F,y.ticklabs="",z.ticklabs="",xlim=c(0,200),color=as.numeric(indMDS$domtax.print),pch=16,scale.y=0.5)
   #In order to get the text into the plot (RV species), I need to convert from xyz to xy (found in a forum)
   s3d.coords <- s3d$xyz.convert(indMDS[,c("viral_colldy","V2","V3")])
   #Insert the RV text
   text(s3d.coords$x,s3d.coords$y,labels=RV.plot.var,pos=3)
   #Insert the indicators of exacerbation below the point (see "-0.02"), and change the point type (pch)
   s3d$points3d(indMDS$viral_colldy[indMDS$exacid>0],indMDS$V2[indMDS$exacid >0],indMDS$V3[indMDS$exacid>0]-0.02,pch=17,cex=1.5)
  # s3d$points3d(indMDS[,c("viral_colldy","V2","V3")],pch=indmap$RV.comb)
   #Include a legend for now (maybe take out later?)
   legend("topright",legend=unique(as.character(indMDS$domtax.print)),col=as.numeric(unique(indMDS$domtax.print)),pch=16,cex=0.9)
}
dev.off()


#Set up a fake plot to print out a legend with all dominant taxa in the individual plots
setwd("/data/Users/kmccauley/PROSE_NEW/Figures")
png(paste0("legend_",seedno,".png"),height=800,width=400)
par(mar=c(2.1,2.1,2.1,2.1))
plot(0,0,type="n",asp=5,xlab="",ylab="",col="white",col.axis="white",col.lab="white",col.sub="white",axes=FALSE)
legend("center",legend=unique(as.character(sort(map2$domtax.print))),col=as.numeric(unique(sort(map2$domtax.print))),pch=15,cex=1.2)
dev.off()

