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

if(pacman::p_version(vegan) == '2.5.3') {
   pacman::p_load(vegan, parallel, ape, ggplot2, lmerTest)
} else {
   pacman::p_install(vegan)
   pacman::p_load(parallel, ape, ggplot2, lmerTest)
}

loopAdonis <- function(datafile,metadata,workdir,use_OTUtable=FALSE, method="bray", perms=1000,sampvar="#SampleID", save=TRUE, outputfile="adonis_results.csv",cores=1, verbose=TRUE, subjid, print_plot=TRUE, subset_var=NULL, subset_group=NULL) {
   #First we set our working directory
   setwd(workdir)
   #Bring in the distance matrix/OTU table
   if(readLines(datafile, n=1) == "# Constructed from biom file") {
   dfile <- read.table(datafile,sep="\t",header = T, check.names = F, comment="",row.names=1, skip=1)
   } else {
   dfile <- read.table(datafile, header=T, sep="\t", check.names=F, comment="", row.names=1)
   }

   #Bring in the metadata
   MYmetaEF = read.table(metadata, header = T, sep = "\t", check.names = F,comment="")
	if(!is.null(subset_var)) MYmetaEF <- MYmetaEF[MYmetaEF[,subset_var] %in% subset_group,]
	if(use_OTUtable == FALSE) method <- NULL
      if(use_OTUtable == FALSE) {
        split <- strsplit(datafile,"/")[[1]]
        distname <- strsplit(split[length(split)],"_")[[1]][1]
        }
   var.tab <- if(class(try(read.table("Variable_Table.txt"),silent=FALSE)) %in% "try-error") {
      varlist <- names(MYmetaEF)[!names(MYmetaEF) %in% sampvar]
      classlist <- sapply(lapply(MYmetaEF[!names(MYmetaEF) %in% sampvar], class), paste)
      varlist <- cbind(varlist, as.character(classlist))
      print("Printing list of variables to be run. Edit the file to include only variables of interest.")
      write.table(varlist, "Variable_Table.txt", col.names=F, row.names=F, sep="\t", quote=F)
      stop()
   } else {
      varlist <- read.table("Variable_Table.txt", header=F, comment="", sep="\t", check.names=F)
   }
      
   #Create the names of the columns for the output that we will generate
   col_names = list("samptot","R2", "pvalue", "pval.pc1","pval.pc2")
   #Create an empty matrix
   p = matrix(data=NA, nrow=nrow(varlist), ncol=length(col_names), dimnames=list(varlist[,1], col_names))
	dir.create(paste0(workdir, "_",method, distname, "_", subset_var,subset_group,"/"))
	setwd(paste0(workdir, "_",method, distname, "_", subset_var,subset_group,"/"))
   for (i in varlist[,1]){
        call <- paste("as", varlist[varlist[,1] %in% i,2], sep = ".")
        MYmetaEF[,i] <- do.call(call, list(MYmetaEF[[i]]))
        #This code generates information about the variables that will be run through adonis
	if(verbose == TRUE) {
        	if(class(MYmetaEF[,i]) %in% "factor") {
        	   print(paste(i, "is a(n)", class(MYmetaEF[,i]), "variable with levels:", paste(levels(MYmetaEF[,i]), collapse=",")))
        	} else {
		 if(class(MYmetaEF[,i]) %in% "integer" | class(MYmetaEF[,i]) %in% "numeric") {
        	   print(paste(i, "is a(n)", class(MYmetaEF[,i]), "variable with values: ", paste(unique(MYmetaEF[,i]), collapse=",")))
        	} else {
		  print(paste(i, "is a(n)", class(MYmetaEF[,i]), "variable"))
	}
	}
	}

      #print("test1")
      #Get a list of the ids where the metadata variable isn't missing
      ids <- as.character(MYmetaEF[!is.na(MYmetaEF[[i]]), sampvar])
      if(use_OTUtable == FALSE) {
         sub_MYdat2 <- dfile[ids[ids %in% names(dfile)], ids[ids %in% names(dfile)]]
      } else {
         sub_MYdat2 <- dfile[, ids[ids %in% names(dfile)]]
      }
      #print("test2") 
      set.seed(123) #Make sure your results are reproducible
      #Only do this for variables with more than one level
      MetaEFSubset <- MYmetaEF[as.character(MYmetaEF[,sampvar]) %in% names(sub_MYdat2),] #Create your subset dataset
      if(length(table(factor(MetaEFSubset[[i]])))>1 & sum(dim(sub_MYdat2)>2)==2 ) {
         if(use_OTUtable == FALSE) {
	# For use of a distance matrix:
         sub_MYdat2 <- as.dist(sub_MYdat2)
         l <- adonis2(as.formula(paste0("sub_MYdat2 ~ ", i )), data=MetaEFSubset, permutations = perms,parallel=cores) #Adonis from Distance Matrix
         p[i,"samptot"] <- dim(MetaEFSubset)[1]
	 pcs <- pcoa(sub_MYdat2)$vectors[,c("Axis.1","Axis.2")]
	 MetaEFSubset.pcs <- merge(MetaEFSubset, pcs, by.x=sampvar, by.y=0)
	plt <- ggplot(MetaEFSubset.pcs, aes(x=Axis.1, y=Axis.2, color=MetaEFSubset.pcs[,i])) + 
		geom_point() +
		guides(color=guide_legend(title=paste(i)))
	if(print_plot == TRUE) ggsave(paste0(i,"_",distname, "_PCOAplot.pdf"), plt)
	if(sum(duplicated(MetaEFSubset.pcs[,subjid]))>0) {
	 lme.pc1 <- lmer(as.formula(paste0("Axis.1 ~", i , "+ (1|", subjid, ")")), data=MetaEFSubset.pcs)
	 lme.pc2 <- lmer(as.formula(paste0("Axis.2 ~", i , "+ (1|", subjid, ")")), data=MetaEFSubset.pcs)
	 p[i, "pval.pc1"] <- anova(lme.pc1)[,"Pr(>F)"]
	 p[i, "pval.pc2"] <- anova(lme.pc2)[,"Pr(>F)"]
	} else {
         p[i,"R2"] = l$R2[1] #print the R^2 value in the second column
         p[i,"pvalue"] = l$`Pr(>F)`[1] #print the p-value in the third column
	}
         } else { #For the use of an OTU table
         t.MYdat2 <- vegdist(t(sub_MYdat2), method=method)
         l <- adonis2(as.formula(paste0("t.MYdat2 ~ ", i)), data=MetaEFSubset, permutations=perms, parallel=cores)
         p[i,"samptot"] <- dim(MetaEFSubset)[1]
	pcs <- pcoa(t.MYdat2)$vectors[,c("Axis.1","Axis.2")]
	MetaEFSubset.pcs <- merge(MetaEFSubset.pcs, by.x=sampvar, by.y=0)
	if(sum(duplicated(MetaEFSubset.pcs[,subjid])>0)) {
         lme.pc1 <- lmer(as.formula(paste0("Axis.1 ~", i , "+ (1|", subjid, ")")), data=MetaEFSubset.pcs)
         lme.pc2 <- lmer(as.formula(paste0("Axis.2 ~", i , "+ (1|", subjid, ")")), data=MetaEFSubset.pcs)
         p[i, "pval.pc1"] <- anova(lme.pc1)[,"Pr(>F)"]
         p[i, "pval.pc2"] <- anova(lme.pc2)[,"Pr(>F)"]
	} else {
         p[i,"R2"] = l$R2[1] #print the R^2 value in the second column
         p[i,"pvalue"] = l$`Pr(>F)`[1] #print the p-value in the third column
	}
	plt <- ggplot(MetaEFSubset.pcs, aes(x=Axis.1, y=Axis.2, color=MetaEFSubset.pcs[,i])) +
		geom_point() +
		guides(color=guide_legend(title=paste(i)))
	if(print_plot ==TRUE) ggsave(paste0(i,"_", method, "_PCOAplot.pdf"), plt)
         }
      }
   }
   write.csv(p,file=paste0("results_", method, distname,"_", subset_var,subset_group,".csv")) #this writes a csv, and will call it whatever you choose in your function (assuming you use save = TRUE)
}

loopAdonis(datafile="/data/Users/kmccauley/CREW/DADA2_Comb/bdiv_dada2_22582/bray_curtis_dm.txt",
        metadata="/data/Users/kmccauley/CREW/CREW_Data_for_Brandon_Julia/CREW_TwoCohort_Combined.txt",
        subset_var=NULL,
        subset_group=NULL,
        workdir="/data/Users/kmccauley/CREW/BetaDiversity/STOOL/",
        sampvar="Stool_SampleID",
        subjid="PairID",
        perms=999,
        cores=15,
        use_OTUtable=FALSE,
        print_plot=TRUE,
        verbose=TRUE, # Information about variables will now only print if verbose == TRUE. This way, you can have multiple loopAdonis functions and only print out the information for the first run.
        method="canberra")

#A simple run of the code when using an OTU table and not a distance matrix
