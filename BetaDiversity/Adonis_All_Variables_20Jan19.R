#Code to Run Adonis from Distance Matrix or OTU table
#July 06: Code updated to use the new "vegan" package, and specifically 'adonis2'. Additionally, features to print out information about variable format and sample total have been included.
#Date Last Updated: July 06, 2018
#Katie McCauley
# To report issues with this script, go to: https://ucsf.box.com/s/300pq1tal5ottrhvund41qyjr0sgejqm

#library(vegan) #initialize the 'vegan' package (you may need to install it first: it won't initialize if it's not downloaded)
#library(parallel) # allows parallelization

#Loading packages differently, ensuring the most-recent version is installed (07/06):
# First confirm that pacman is available, and if not, install it:
if(!require("pacman")) install.packages("pacman")

if(pacman::p_version(vegan) == '2.5.6') {
   pacman::p_load(phyloseq, vegan, parallel, ape, ggplot2, lmerTest)
} else {
   pacman::p_install(vegan)
   pacman::p_load(parallel, ape, ggplot2, lmerTest)
}

loopAdonis <- function(phy, workdir=".", method="bray", perms=999, save=TRUE, outputfile="adonis_results.csv", cores=1, verbose=TRUE, subjid, print_plot=TRUE, subset_var=NULL, subset_group=NULL) {
   #First we set our working directory
  require(phyloseq)
   setwd(workdir)
   #Bring in the distance matrix/OTU table
   phyobj <- phy
   print(phyobj)
   MYmetaEF <- phyobj@sam_data
   var.tab <- if(class(try(read.table("Variable_Table.txt"),silent=FALSE)) %in% "try-error") {
      varlist <- names(MYmetaEF)
      classlist <- sapply(lapply(MYmetaEF, class), paste)
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
	dir.create(paste0(workdir,method, "_", subset_var,subset_group,"/"), showWarnings = FALSE)
	setwd(paste0(workdir,method, "_", subset_var,subset_group,"/"))
   for (i in varlist[,1]){
        call <- paste("as", varlist[varlist[,1] %in% i,2], sep = ".")
        sample_data(phyobj)[[i]] <- do.call(call, list(data.frame(sample_data(phyobj))[,i]))
        #This code generates information about the variables that will be run through adonis
	if(verbose == TRUE) {
        	if(class(data.frame(sample_data(phyobj))[,i]) %in% "factor") {
        	   print(paste(i, "is a(n)", class(data.frame(sample_data(phyobj))[,i]), "variable with levels:", paste(levels(data.frame(sample_data(phyobj))[,i]), collapse=",")))
        	} else {
		 if(class(data.frame(sample_data(phyobj))[,i]) %in% "integer" | class(data.frame(sample_data(phyobj))[,i]) %in% "numeric") {
        	   print(paste(i, "is a(n)", class(data.frame(sample_data(phyobj))[,i]), "variable with values: ", paste(unique(data.frame(sample_data(phyobj))[,i]), collapse=",")))
        	} else {
		  print(paste(i, "is a(n)", class(data.frame(sample_data(phyobj))[,i]), "variable"))
	}
	}
	}

        phyobj.ss <- subset_samples(phyobj, !is.na(sample_data(phyobj)[[i]]) & sample_sums(phyobj) > 0)
        phyobj.ss
      #print("test1")
      #Get a list of the ids where the metadata variable isn't missing
      #print("test2") 
      set.seed(123) #Make sure your results are reproducible
      #Only do this for variables with more than one level

      if(nrow(unique(sample_data(phyobj.ss)[,i]))>1 & sum(duplicated(sample_data(phyobj.ss)[,i])) > 0 & nsamples(phyobj.ss) > 2) {
	# For use of a distance matrix:
        distobj <- phyloseq::distance(phyobj.ss, method=method)
           l <- adonis2(distobj ~ data.frame(sample_data(phyobj.ss))[,i], permutations=perms, parallel=cores)

           p[i,"samptot"] <- nsamples(phyobj.ss)
	 pcs <- pcoa(distobj)$vectors[,c("Axis.1","Axis.2","Axis.3")]
	 MetaEFSubset.pcs <- data.frame(merge(sample_data(phyobj.ss), pcs, by=0), check.names=FALSE)
	 
	plt <- ggplot(MetaEFSubset.pcs, aes(x=Axis.1, y=Axis.2, color=MetaEFSubset.pcs[,i])) + 
		geom_point() +
		guides(color=guide_legend(title=paste(i)))
	
	if(print_plot == TRUE) ggsave(paste0(i,"_",method, "_PCOAplot.pdf"), plt)
	
	if(sum(duplicated(sample_data(phyobj.ss)[,subjid]))>0) {
	 lme.pc1 <- lmer(as.formula(paste0("Axis.1 ~", i , "+ (1|", subjid, ")")), data=MetaEFSubset.pcs)
	 lme.pc2 <- lmer(as.formula(paste0("Axis.2 ~", i , "+ (1|", subjid, ")")), data=MetaEFSubset.pcs)
         lme.pc3 <- lmer(as.formula(paste0("Axis.3 ~", i , "+ (1|", subjid, ")")), data=MetaEFSubset.pcs)
	 p[i, "pval.pc1"] <- anova(lme.pc1)[,"Pr(>F)"]
	 p[i, "pval.pc2"] <- anova(lme.pc2)[,"Pr(>F)"]
         p[i, "pval.pc3"] <- anova(lme.pc3)[,"Pr(>F)"]
	} else {
         p[i,"R2"] = l$R2[1] #print the R^2 value in the second column
         p[i,"pvalue"] = l$`Pr(>F)`[1] #print the p-value in the third column
	}
  }
  }
write.csv(p,file=paste0("results_", method,"_", subset_var,subset_group,".csv")) #this writes a csv, and will call it whatever you choose in your function
}


#load("/data/Users/kmccauley/URECA_STOOL/bad_data/dada2_phy_init.rds")
#loopAdonis(phyobj=phy_filt,
#        subset_var=NULL,
#        subset_group=NULL,
#        workdir="/data/Users/kmccauley/URECA_STOOL/",
#        subjid="RecruitID",
#        perms=999,
#        cores=2,
#        print_plot=TRUE,
#        verbose=TRUE, # Information about variables will now only print if verbose == TRUE. This way, you can have multiple loopAdonis functions and only print out the information for the first run.
#        method="canberra")

