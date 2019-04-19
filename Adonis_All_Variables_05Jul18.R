#Code to Run Adonis from Distance Matrix or OTU table
   #Code CANNOT compute PERMANOVAs from OTU tables for Weighted or Unweighted UniFrac. You *MUST* use a distance matrix as generated in QIIME
#July 06: Code updated to use the new "vegan" package, and specifically 'adonis2'. Additionally, features to print out information about variable format and sample total have been included.
#Date Last Updated: July 06, 2018
#Katie McCauley
# To report issues with this script, go to: https://ucsf.box.com/s/300pq1tal5ottrhvund41qyjr0sgejqm

rm(list=ls()) #Clean the working environment
#library(vegan) #initialize the 'vegan' package (you may need to install it first: it won't initialize if it's not downloaded)
#library(parallel) # allows parallelization

#Loading packages differently, ensuring the most-recent version is installed (07/06):
# First confirm that pacman is available, and if not, install it:
if(!require("pacman")) install.packages("pacman")

if(pacman::p_version(vegan) == '2.5.2') {
   pacman::p_load(vegan, parallel)
} else {
   pacman::p_install(vegan)
   pacman::p_load(parallel)
}

loopAdonis <- function(datafile,metadata,workdir,use_OTUtable=FALSE, method="bray", perms=1000,sampvar="#SampleID", save=TRUE, outputfile="adonis_results.csv",cores=1, verbose=TRUE) {
   #First we set our working directory
   setwd(workdir)
   #Bring in the distance matrix/OTU table
   if(readLines(datafile, n=1) == "# Constructed from biom file") {
   dfile <- read.table(datafile,header = T, sep = "\t", check.names = F, comment="",row.names=1, skip=1)
   } else {
   dfile <- read.table(datafile, header=T, sep="\t", check.names=F, comment="", row.names=1)
   }

   #Bring in the metadata
   MYmetaEF = read.table(metadata, header = T, sep = "\t", check.names = F,comment="")
   var.tab <- if(class(try(read.table("Variable_Table.txt"),silent=FALSE)) %in% "try-error") {
      varlist <- names(MYmetaEF)[!names(MYmetaEF) %in% sampvar]
      print("Printing list of variables to be run. Edit the file to include only variables of interest.")
      write.table(varlist, "Variable_Table.txt", col.names=F, row.names=F, sep="\t", quote=F)
      stop()
   } else {
      varlist <- read.table("Variable_Table.txt", header=F, comment="", sep="\t", check.names=F)
      varlist <- as.vector(varlist[,1])
   }
      
   #Create the names of the columns for the output that we will generate
   col_names = list("samptot","R2", "pvalue")
   #Create an empty matrix
   p = matrix(data=NA, nrow=length(varlist), ncol=length(col_names), dimnames=list(varlist, col_names))
   #Start the for-loop (going from the first column you are interested in to the last column in the metadata file)
   for (i in varlist){

        #This code generates information about the variables that will be run through adonis
	if(verbose == TRUE) {
        	if(class(MYmetaEF[,i]) %in% "factor") {
        	   print(paste(i, "is a(n)", class(MYmetaEF[,i]), "variable with levels:", paste(levels(MYmetaEF[,i]), collapse=",")))
        	} else {
        	   print(paste(i, "is a(n)", class(MYmetaEF[,i]), "variable."))
        	}
	}

      #print("test1")
      #Get a list of the ids where the metadata variable isn't missing
      ids <- as.character(MYmetaEF[!is.na(MYmetaEF[[i]]), sampvar])
      if(use_OTUtable == FALSE) {
         sub_MYdat2 <- dfile[ids, ids]
      } else {
         sub_MYdat2 <- dfile[, ids]
      }
      #print("test2") 
      set.seed(123) #Make sure your results are reproducible
      #Only do this for variables with more than one level
      MetaEFSubset <- MYmetaEF[as.character(MYmetaEF[,sampvar]) %in% as.character(ids),] #Create your subset dataset
      if(length(table(MetaEFSubset[[i]]))>1) {
         if(use_OTUtable == FALSE) {
         sub_MYdat2 <- as.dist(sub_MYdat2)
         l <- adonis2(as.formula(paste0("sub_MYdat2 ~ ", i )), data=MetaEFSubset, permutations = perms,parallel=cores) #Adonis from Distance Matrix
         } else {
         t.MYdat2 <- t(sub_MYdat2)
         l <- adonis2(as.formula(paste0("t.MYdat2 ~ ", i)), data=MetaEFSubset, method=method, permutations=perms, parallel=cores) #If using an OTU table
         }
         #if(verbose==TRUE) {print (l)} #Print out the ADONIS output for i'th variable
         #m = l$aov.tab ## This is no longer needed since the output has been restructured in this updated version of "vegan"
         
        # p[i,1] = i #print the name of the variable in the first column #
         p[i,"samptot"] <- length(ids)
         p[i,"R2"] = l$R2[1] #print the R^2 value in the second column
         p[i,"pvalue"] = l$`Pr(>F)`[1] #print the p-value in the third column
      }
   }
   write.csv(p,file=outputfile) #this writes a csv, and will call it whatever you choose in your function (assuming you use save = TRUE)
}

#A simple run of the code when using an OTU table and not a distance matrix
loopAdonis(datafile="../../OTUtables/Gelatin_Added/rarefied_ureca_nasal_18131.txt",
	metadata="../FinalMapFile_URECA_Jun8_2018.txt",
	workdir="/data/Users/kmccauley/URECA_Nasal_Dust/Analysis/CommunityComp/",
	sampvar="#SampleID",
	perms=999, 
	outputfile="adonis_results.csv",
	cores=15,
	use_OTUtable=TRUE,
	verbose=TRUE, # Information about variables will now only print if verbose == TRUE. This way, you can have multiple loopAdonis functions and only print out the information for the first run.
	method="bray")

