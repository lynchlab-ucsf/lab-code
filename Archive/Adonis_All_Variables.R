#Code to Run Adonis from a QIIME distance matrix
   #Can be edited to run from an OTU table, but I wanted to keep the tutorial simple(ish)
   #OTU tables are harder, especially if we want flexibility for UniFrac as well as Canberra/Bray-Curtis
#Code incorporates ideas in scripts from Doug, Kei, and Meera (thank you!)
#Shows use of functions as well as for-loops
#Date: July 12, 2016
#Katie McCauley


rm(list=ls()) # I like to do this before every program I write. It removes everything from my environment
library(vegan) #initialize the 'vegan' package (you may need to install it first: it won't initialize if it's not downloaded)
library(parallel) # allows parallelization

#for troubleshooting/practice
#These files can be found in the Box folder
skipcol <- 5  # Skipcol is the number of columns to skip from meta-data file
datafile <- "BetaDiv_qiime/weighted_unifrac_dm.txt"
metadata <- "FINAL_MAP_FILE.txt"
workdir <- "/data/Users/kmccauley/PROSE_NEW/AnalysisData"
i <- 2

loopAdonis <- function(skipcol,datafile,metadata,workdir,perms=1000,sampvar="#SampleID", save=FALSE, outputfile,cores=1,skip_first_line=TRUE) {
   #First we set our working directory
   setwd(workdir)
   #Bring in the distance matrix
   dmat <- read.table(datafile,header = T, sep = "\t", check.names = F, comment.char="",row.names=1)
   #Bring in the metadata
   if(skip_first_line == FALSE) {
   MYmetaEF = read.table(metadata, header = T, sep = "\t", check.names = F,comment.char="") #Changed to check.names=F because I want to keep "#SampleID"
   } else {
   MYmetaEF = read.table(metadata, header = T, sep = "\t", check.names = F,comment.char="",skip=1) #Changed to check.names=F because I want to keep "#SampleID"
   }
   #Create the names of the columns for the output that we will generate
   col_names = list("parameter", "R2", "pvalue")
   #Create an empty matrix
   p = matrix(data=NA, nrow=length(colnames(MYmetaEF))-skipcol, ncol=3, dimnames=list(1:(length(colnames(MYmetaEF))-skipcol), col_names))
   #Start the for-loop (going from the first column you are interested in to the last column in the metadata file)
   for (i in (skipcol+1):length(colnames(MYmetaEF))){
      #Earlier, I tested this script on my PROSE dataset, and the script was taking so long that I wasn't sure it was running
      #print("test1")
      #Get a list of the ids where the metadata variable isn't missing
      ids <- as.character(MYmetaEF[!is.na(MYmetaEF[[i]]),1])
      sub_MYdat <- subset(dmat, select=ids)		# Subset the original untransformed data by selecting samples with no NAs
      MYarray2 <- t(sub_MYdat)				    # Transpose the subset data-matrix
      sub_MYdat2 <- subset(MYarray2, select=ids) # Again, select samples with no NAs from the other side of the distance matrix
      #print("test2") 
      set.seed(123) #Make sure your results are reproducible
      #Only do this for variables with more than one level
      MetaEFSubset <- MYmetaEF[as.character(MYmetaEF[,sampvar]) %in% as.character(ids),] #Create your subset dataset
      if(length(table(MetaEFSubset[[i]]))>1) {
         l <- adonis(as.dist(sub_MYdat2) ~  MetaEFSubset[[i]], data=MetaEFSubset, permutations = perms,parallel=cores) #Run ADONIS
         print (l) #Print out the ADONIS output for i'th variable
         m = l$aov.tab #Subset out the section with results that we're interested in
         #print("test3")
         p[i-skipcol,1] = colnames(MYmetaEF)[i] #print the name of the variable in the first column
         p[i-skipcol,2] = m[1,5] #print the R^2 value in the second column
         p[i-skipcol,3] = m[1,6] #print the p-value in the third column
         print(p) #print the p data-frame to be kept apprised of progress.
      } else {
         p[i-skipcol,1] = colnames(MYmetaEF)[i]
      }
   }
if (save==TRUE) { #If you would like to save your work, you can choose to do so from your function
   write.csv(p,file=outputfile) #this writes a csv, and will call it whatever you choose in your function (assuming you use save = TRUE)
}
}
#A simple run of the code
loopAdonis(5,"BetaDiv_qiime/weighted_unifrac_dm.txt","FINAL_MAP_FILE.txt","/data/Users/kmccauley/PROSE_NEW/AnalysisData", sampvar="#SampleID", perms=50,save=TRUE,outputfile="/data/Users/kmccauley/PROSE_NEW/AnalysisData/adonis_results.csv",cores=15)

