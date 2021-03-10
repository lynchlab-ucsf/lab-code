#Code to Run Adonis from Distance Matrix or OTU table
   #Code CANNOT compute PERMANOVAs from OTU tables for Weighted or Unweighted UniFrac. You *MUST* use a distance matrix as generated in QIIME
#Code incorporates ideas in scripts from Doug, Kei, and Meera (thank you!)
#Date Last Updated: November 26, 2016
#Katie McCauley


rm(list=ls()) #Clean the working environment
library(vegan) #initialize the 'vegan' package (you may need to install it first: it won't initialize if it's not downloaded)
library(parallel) # allows parallelization

loopAdonis <- function(skipcol,datafile,metadata,workdir,use_OTUtable=FALSE, method="bray", perms=1000,sampvar="#SampleID", save=FALSE, outputfile,cores=1,skip_first_line=TRUE,verbose=TRUE) {
   #First we set our working directory
   setwd(workdir)
   #Bring in the distance matrix/OTU table
   if(skip_first_line == TRUE) {
   dfile <- read.table(datafile,header = T, sep = "\t", check.names = F, comment.char="",row.names=1, skip=1)
   } else {
   dfile <- read.table(datafile, header=T, sep="\t", check.names=F, comment.char="", row.names=1)
   }
   #Bring in the metadata
   MYmetaEF = read.table(metadata, header = T, sep = "\t", check.names = F,comment.char="") #Changed to check.names=F because I want to keep "#SampleID"
   
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
      if(use_OTUtable == FALSE) {
         sub_MYdat <- subset(dfile, select=ids)		# Subset the original untransformed data by selecting samples with no NAs
         MYarray2 <- t(sub_MYdat)				    # Transpose the subset data-matrix
         sub_MYdat2 <- subset(MYarray2, select=ids) # Again, select samples with no NAs from the other side of the distance matrix
      } else {
         sub_MYdat2 <- subset(dfile, select=ids)
      }
      #print("test2") 
      set.seed(123) #Make sure your results are reproducible
      #Only do this for variables with more than one level
      MetaEFSubset <- MYmetaEF[as.character(MYmetaEF[,sampvar]) %in% as.character(ids),] #Create your subset dataset
      if(length(table(MetaEFSubset[[i]]))>1) {
         if(use_OTUtable == FALSE) {
         sub_MYdat2 <- as.dist(sub_MYdat2)
         l <- adonis(sub_MYdat2 ~  MetaEFSubset[[i]], data=MetaEFSubset, permutations = perms,parallel=cores) #Adonis from Distance Matrix
         } else {
         t.MYdat2 <- t(sub_MYdat2)
         l <- adonis(t.MYdat2 ~ MetaEFSubset[[i]], data=MetaEFSubset, method=method, permutations=perms, parallel=cores) #If using an OTU table
         }
         if(verbose==TRUE) {print (l)} #Print out the ADONIS output for i'th variable
         m = l$aov.tab #Subset out the section with results that we're interested in
         #print("test3")
         p[i-skipcol,1] = colnames(MYmetaEF)[i] #print the name of the variable in the first column
         p[i-skipcol,2] = m[1,5] #print the R^2 value in the second column
         p[i-skipcol,3] = m[1,6] #print the p-value in the third column
         if(verbose==TRUE) {print(p)} #print the p data-frame to be kept apprised of progress.
      } else {
         p[i-skipcol,1] = colnames(MYmetaEF)[i]
      }
   }
if (save==TRUE) { #If you would like to save your work, you can choose to do so from your function
   write.csv(p,file=outputfile) #this writes a csv, and will call it whatever you choose in your function (assuming you use save = TRUE)
}
}
#A simple run of the code
loopAdonis(skipcol=5,datafile="PROSE_OTUtable_subset.txt",metadata="FINAL_MAP_FILE_subset.txt",workdir="/data/Users/kmccauley/PROSE_NEW/AnalysisData", sampvar="#SampleID", perms=999,save=TRUE,outputfile="/data/Users/kmccauley/PROSE_NEW/AnalysisData/adonis_results.csv",cores=15, use_OTUtable=TRUE,method="canberra",skip_first_line=TRUE,verbose=FALSE)

