#Code to run the three-model OTU Differential-Expression Analysis
#Written by : Katie McCauley, loosely based on GitHub scripts by Ali Faruqi

#In the example below, I am looking for OTUs that are differentially-abundant between people who ever had an exacerbation, compared to those who did not have an exacerbation, using the three-model approach. I also have the ability to look at subsets of my data (ie only within Xolair) if I want to.

#In the PROSE data, we had a lot of trouble with the mixed-effects models converging due to too many zeros. This issue was addressed using a few different data-driven decisions in consultation with Rho World. First, we dropped any results from models that did not converge. We saw that many models that produced warnings often also produced incorrect AICs that were so low they were chosen as the model of best fit, despite being completely divergent from other models that R determined fit correctly. Therefore, we dropped (and did not consider) any models that did not converge. In addition, before the OTU could go through the model-building process, at least a certain proportion of the cells (here, 15%) had to contain non-zero values. We determined this cutoff by plotting the log of the total reads for each OTU against the proportion of samples with zero-counts. The plot should result in something that looks like a half parabola, in which the cutoff should be the point at which the parabola's edge starts exponentially increasing (Code to make this plot is in the code below, but commented out). Rho used a different method, but we independently chose the same cutoff.

#To run this script in the background, type the following (without pound sign):
   # nohup Rscript {filename}.R &
      # In order to continue to use the command promt, type CTRL-C. The command above will also open up a filename called "nohup.out" where all of the comments are printed.
   #Also:
   # nohup Rscript {filename}.R &> myout.out &
      #This will return you to the command line and place all command output into a file called "myout.out"
#########June 5th Update################################################
#I added the ability to define what your SampleID variable looks like.

########################################################################
#install.packages(c("nlme","doMC","foreach","MASS","pscl"))
#To install "glmmADMB", see: http://glmmadmb.r-forge.r-project.org/
        ##Type:
#install.packages("R2admb")
#install.packages("glmmADMB", 
#    repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#            getOption("repos")),
#    type="source")

#Set up potentially-used libraries:
require(pscl)
require(MASS)
require(foreach)
require(doMC)
require(lme4)
require(glmmADMB)


rm(list=ls())
workdir <- "." #Directory where you want to save data
otutable <- "OTU_rho_hdf5_realtax.txt" #OTU table in text file from BIOM format with "taxonomy" as a column at some point (not necessarily at the end)
mapfile <- "ureca.map_040416_R.txt" #Map file, typically in QIIME format, with "#SampleID" as the first variable
treatment <- "asthma_agey7" #Variable (from map file) of comparison
trt1 <- "No" #Reference group in this analysis
trt2 <- "Yes"
outputname <- "ThreeModel_URECA_dup.csv" #CSV filename where results will be sent to
imagename <- "rawdata_ureca" #Because computation time can be extensive, and the model-selection step can break, I save the intermediate working enviornment under this name to be able to pull it up and troubleshoot after the parallel part
subset <- FALSE # Do you want to perform subset analysis? In the prose analysis, I wanted results among Placebo, ICS, and Xolair separately
ss_var <- "group" # Group that you will subset by (ie, only want the results from a subset of individuals)
ss_group <- "Xolair" # The specific group that you want results for

mixed_effects <- FALSE #Do you want to run a mixed-effects (ME) analysis (ie, theoretically-correlated samples)
ind_id <- "studyid" #Individual ID: Set to NULL if not using ME models

drop_warnings <- TRUE #If an OTU's model produces a convergence warning/error, don't include the results -- we're making an assumption that if the OTU doesn't converge, then the resulting coefficient and p-value are not biologically relevant

cores <- 15 #Number of server cores to use (max on lynchserver2 is 80, but it's adviseable not to go that high)
cutoff <- 0.15

diag_plot <- TRUE

####### Additional Commands
choose_by_BIC <- FALSE #choose the best of the three models using the AIC or BIC. If choose_by_BIC == FALSE, then the script will use the AIC to choose the best model.
skip_first <- TRUE #Skip the first line when reading in your OTU table

sampid <- "#SampleID"#Sometimes different datasets use different sampleID variables

###################
# A function for weighted means by each individual
wght_mean_byind <- function(num,id_var) {
dat <- data.frame(num,id_var)
ind.wt <- 1/table(id_var)[table(id_var) != 0]
num.ind <- length(ind.wt)
indwt.df <- data.frame(ind.wt)
wt.dat <- merge(dat,indwt.df,by.x="id_var",by.y=0)
wt.mean <- sum(wt.dat$num*wt.dat$ind.wt)/num.ind
wt.mean 
}

# A function to print warnings
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

###################
#Start Code

setwd(workdir)
if(skip_first == TRUE) {
MYdata <- read.table(otutable,header = T,sep = "\t", check.names = F, comment.char= "", quote="",row.names=1,skip=1)
} else {
MYdata <- read.table(otutable,header = T,sep = "\t", check.names = F, comment.char= "", quote="",row.names=1)
}


#The PROSE data had headers, was tab separated, used a "#" sign in front of "SampleID", used quotes, and the rownames (ie, OTUIDs) were in the first column and I didn't want to consider them in the data frame as a separate column
taxa.name <- data.frame(row.names(MYdata), as.character(MYdata$taxonomy))
#Separate out the taxonomy from the data, but bring the OTUIDS with the taxonomies (to correctly merge later)

#Drop taxonomy from the MYdata dataframe
MYdata$taxonomy <- NULL

#Back when the taxonomy pull didn't work, I used this to check the taxonomy of the OTU
#taxa.name[row.names(MYdata) == "_16039"]

#Read in the map file
MYmeta <- read.table(mapfile,header = T, sep = "\t", check.names = F, comment.char= "")
if(subset == TRUE) {
MYmeta <- MYmeta[MYmeta[,ss_var] == ss_group,]
}
allcols <- length(colnames(MYdata)) #Get the number of OTUs in your data
matrix1.subset <- !is.na(MYmeta[,treatment]) & MYmeta[,treatment]==trt1 #TRUE/FALSE statement for non-missing treatment variable, and being in group 1
matrix1 <- MYdata[c(as.character(MYmeta[,sampid][matrix1.subset]))] # pull out only the SampleIDs that have data on your variable of interest, based on T/F
matrix2.subset <- !is.na(MYmeta[,treatment]) & MYmeta[,treatment]==trt2 #TRUE/FALSE statement for non-missing treatment variable and being in group 2
matrix2 <- MYdata[c(as.character(MYmeta[,sampid][matrix2.subset]))] # pull out only the SampleIDs that have data on your variable of interest, based on T/F
matrix3 <- cbind(matrix1,matrix2) #bring two matrices together
mat.array <- t(matrix3) #transform so OTUs are columns
if(mixed_effects == TRUE) {
MYmeta2 <- MYmeta[,c(sampid,ind_id,treatment)] #The non-OTU metadata needed for the upcoming loop
both <- merge(mat.array,MYmeta2,by.x=0,by.y=sampid) #Add to OTU data
both[,ind_id] <- factor(as.character(both[,ind_id])) #LME models need studyid to be a factor
} else {
MYmeta2 <- MYmeta[,c(sampid,treatment)] #The non-OTU metadata needed for the upcoming loop
both <- merge(mat.array,MYmeta2,by.x=0,by.y=sampid) #Add to OTU data
}
both$Row.names <- NULL #by.x=0 (above) creates a sudo variable

#PROSE Data had a few issues, primarily that we had LOTS of zero counts, so we needed to determine what proportion of zeros was acceptable to be modeled
#Diagnostic Plot
if(diag_plot == TRUE) {
prev0 <- apply(both[1:length(rownames(MYdata))],2,function(x) {sum(x == 0)/length(x)})
tot.abund <- tot.abund <- apply(both[1:length(rownames(MYdata))],2,function(x) {sum(x)})
pdf("OTU_Cutoff_Diagnostic.pdf")
plot(log(tot.abund),prev0,ylab="Proportion of Samples with Zero Counts",xlab="Log of the Total Reads")
dev.off()
}

registerDoMC(cores)

if(drop_warnings == TRUE) {
options(warn=2) #If the model produces a warning (ie, doesn't converge), then it will treat it like an error and not produce results
}

if(mixed_effects == TRUE) {
   usePackages <- c('lme4','glmmADMB')
} else {
   usePackages <- NULL
}

all_data <- foreach(i=1:(length(rownames(MYdata))),.combine='rbind',.packages=usePackages) %dopar% {
  if(sum(both[,i] != 0)/length(both[,i]) > cutoff) {
OTUname <- colnames(both)[i]
if(mixed_effects == TRUE) {
  formula1 <- as.formula(paste0("both[,",i,"] ~ ",treatment," + (1|",ind_id,")"))
  formula2 <- as.formula(paste0("both[,",i,"] ~ ",treatment))
  if(drop_warnings == TRUE) {
    result.lm <- tryCatch(lmer(formula1, data=both),error=function(e) NA)
    result.pois <- tryCatch(glmer(formula1, family=poisson, data = both),error=function(e) NA)
    result.nb <- tryCatch(glmer.nb(formula1, data = both),error=function(e) NA)
    result.zinb <- tryCatch(glmmadmb(formula2, random=~1|studyid, data = both,link="log",family="nbinom",zeroInflation=TRUE),error=function(e) NA)
  } else {
    lm <- myTryCatch(lmer(formula1, data=both))
    pois <- myTryCatch(glmer(formula1, family=poisson, data = both))
    nb <- myTryCatch(glmer.nb(formula1, data = both))
    zinb <- myTryCatch(glmmadmb(formula2, random=~1|studyid, data = both,link="log",family="nbinom",zeroInflation=TRUE))
    result.lm <- lm$value
    result.pois <- pois$value
    result.zinb <- zinb$value
    result.nb <- nb$value
    warn.lm <- ifelse(is.null(lm$warning),0, as.character(lm$warning))
    warn.pois <- ifelse(is.null(pois$warning),0,as.character(pois$warning))
    warn.zinb <- ifelse(is.null(zinb$warning),0,as.character(zinb$warning))
    warn.nb <- ifelse(is.null(nb$warning),0,as.character(nb$warning))
    error.pois <- ifelse(is.null(pois$error),0,as.character(pois$error))
    error.zinb <- ifelse(is.null(zinb$error),0,as.character(zinb$error))
    error.nb <- ifelse(is.null(nb$error),0,as.character(nb$error))
}
  zinb.coeff <- tryCatch(exp(summary(result.zinb)$coefficients[2,1]),error=function(e) NA)
  zinb.pval <- tryCatch(summary(result.zinb)$coefficients[2,4],error=function(e) NA)
  aic.zinb <- tryCatch(AIC(result.zinb),error=function(e) NA)
  bic.zinb <- tryCatch(AIC(result.zinb,k=log(length(both[,i]))),error=function(e) NA)
  nb.disp <- tryCatch(getME(result.nb,"glmer.nb.theta"), error=function(e) NA) # Negative-Binomial dispersion parameter (from NB model)
  zinb.disp  <- tryCatch(summary(result.zinb)$alpha,error=function(e) NA) #Negative-Binomial dispersion parameter (from ZINB model)
  wgt_mean_trt1 <- tryCatch(wght_mean_byind(both[both[,treatment] == trt1,i],both[both[,treatment] == trt1,ind_id]),error=function(e) NA)
  wgt_mean_trt2 <- tryCatch(wght_mean_byind(both[both[,treatment] == trt2,i],both[both[,treatment] == trt2,ind_id]),error=function(e) NA)
  wgt_mean_diff <- tryCatch((wgt_mean_trt1 - wgt_mean_trt2),error=function(e) NA)
} else {
  formula1 <- as.formula(paste("both[,i] ~ ",treatment," | 1",sep=""))
  formula2 <- as.formula(paste("both[,i] ~ ",treatment,sep=""))
  if(drop_warnings == TRUE) {
    result.lm <- tryCatch(lm(formula2,data=both),error=function(e) NA)
    result.pois <- tryCatch(glm(formula2, family="poisson", data = both),error=function(e) NA)
    result.zinb <- tryCatch(zeroinfl(formula1, data = both, dist = "negbin"),error=function(e) NA)
    result.nb <- tryCatch(glm.nb(formula2, data = both),error=function(e) NA)
  } else {
    lm <- myTryCatch(lm(formula2, data=both))
    pois <- myTryCatch(glm(formula2, family="poisson", data = both))
    zinb <- myTryCatch(zeroinfl(formula1, data = both, dist = "negbin"))
    nb <- myTryCatch(glm.nb(formula2, data = both))
    result.lm <- lm$value
    result.pois <- pois$value
    result.zinb <- zinb$value
    result.nb <- nb$value
    warn.lm <- ifelse(is.null(lm$warning),0,as.character(lm$warning))
    warn.pois <- ifelse(is.null(pois$warning),0,as.character(pois$warning))
    warn.zinb <- ifelse(is.null(zinb$warning),0,as.character(zinb$warning))
    warn.nb <- ifelse(is.null(nb$warning),0,as.character(nb$warning))
    error.lm <- ifelse(is.null(lm$error),0,as.character(lm$error))
    error.pois <- ifelse(is.null(pois$error),0,as.character(pois$error))
    error.zinb <- ifelse(is.null(zinb$error),0,as.character(zinb$error))
    error.nb <- ifelse(is.null(nb$error),0,as.character(nb$error))
  }
  zinb.coeff <- tryCatch(exp(summary(result.zinb)$coef$count[2,1]),error=function(e) NA)
  zinb.pval <- tryCatch(summary(result.zinb)$coef$count[2,4],error=function(e) NA)
  aic.zinb <- tryCatch(AIC(result.zinb),error=function(e) NA)
  bic.zinb <- tryCatch(AIC(result.zinb,k=log(length(both[,i]))),error=function(e) NA)
}
  mean_trt1 <- tryCatch(mean(both[both[,treatment] == trt1,i]),error=function(e) NA)
  mean_trt2 <- tryCatch(mean(both[both[,treatment] == trt2,i]),error=function(e) NA)
  mean_diff <- tryCatch((mean_trt1 - mean_trt2),error=function(e) NA)
  lm.coeff <- tryCatch(exp(summary(result.lm)$coef[2,1]),error=function(e) NA)
  pois.coeff <- tryCatch(exp(summary(result.pois)$coefficients[2,1]),error=function(e) NA)
  pois.pval <- tryCatch(summary(result.pois)$coefficients[2,4],error=function(e) NA)
  nb.coeff <- tryCatch(exp(summary(result.nb)$coefficients[2,1]),error=function(e) NA)
  nb.pval <- tryCatch(summary(result.nb)$coefficients[2,4],error=function(e) NA)
  aic.lm <- tryCatch(AIC(result.lm),error=function(e) NA)
  aic.pois <- tryCatch(AIC(result.pois),error=function(e) NA)
  aic.nb <- tryCatch(AIC(result.nb),error=function(e) NA)
  bic.lm <- tryCatch(AIC(result.lm,k=log(length(both[,i]))),error=function(e) NA)
  bic.pois <- tryCatch(AIC(result.pois,k=log(length(both[,i]))),error=function(e) NA)
  bic.nb <- tryCatch(AIC(result.nb,k=log(length(both[,i]))),error=function(e) NA)
  variance <- var(both[,i])
  trt1vals <- tryCatch(as.numeric(both[,i][which(as.vector(both[,treatment])==trt1)]),error=function(e) NA)
  trt2vals <- tryCatch(as.numeric(both[,i][which(as.vector(both[,treatment])==trt2)]),error=function(e) NA)
  zerotrt1 <- tryCatch(sum(trt1vals == 0),error=function(e) NA)
  zerotrt2 <- tryCatch(sum(trt2vals == 0),error=function(e) NA)
  nonzerotrt1 <- tryCatch(length(trt1vals) - zerotrt1,error=function(e) NA)
  nonzerotrt2 <- tryCatch(length(trt2vals) - zerotrt2,error=function(e) NA)
  totaltrt1 <- tryCatch(sum(trt1vals),error=function(e) NA)
  totaltrt2 <- tryCatch(sum(trt2vals),error=function(e) NA)

if(mixed_effects == TRUE) {
  if(drop_warnings == TRUE) {
  finaldata <- c(OTUname, lm.coeff,lm.pval, pois.coeff,pois.pval,nb.coeff,nb.pval,zinb.coeff,zinb.pval,mean_trt1,mean_trt2,mean_diff,wgt_mean_trt1,wgt_mean_trt2,wgt_mean_diff,aic.lm,aic.pois,aic.nb,aic.zinb,bic.lm,bic.pois,bic.nb,bic.zinb,nb.disp,zinb.disp,zerotrt1,zerotrt2, nonzerotrt1, nonzerotrt2, totaltrt1, totaltrt2,variance)
  } else {
  finaldata <- c(OTUname, lm.coeff,lm.pval, pois.coeff,pois.pval,nb.coeff,nb.pval,zinb.coeff,zinb.pval,mean_trt1,mean_trt2,mean_diff,wgt_mean_trt1,wgt_mean_trt2,wgt_mean_diff,aic.lm,aic.pois,aic.nb,aic.zinb,bic.lm,bic.pois,bic.nb,bic.zinb,nb.disp,zinb.disp,zerotrt1,zerotrt2, nonzerotrt1, nonzerotrt2, totaltrt1, totaltrt2,variance,error.lm,error.pois,error.nb,error.zinb,warn.lm,warn.pois,warn.nb,warn.zinb)
  }
} else {
  if(drop_warnings == TRUE) {
  finaldata <- c(OTUname, lm.coeff,lm.pval,pois.coeff,pois.pval,nb.coeff,nb.pval,zinb.coeff,zinb.pval,mean_trt1,mean_trt2,mean_diff,aic.lm,aic.pois,aic.nb,aic.zinb,bic.lm,bic.pois,bic.nb,bic.zinb,zerotrt1,zerotrt2, nonzerotrt1, nonzerotrt2, totaltrt1, totaltrt2,variance)
  } else {
  finaldata <- c(OTUname, lm.coeff,lm.pval,pois.coeff,pois.pval,nb.coeff,nb.pval,zinb.coeff,zinb.pval,mean_trt1,mean_trt2,mean_diff,aic.lm,aic.pois,aic.nb,aic.zinb,bic.lm,bic.pois,bic.nb,bic.zinb,zerotrt1,zerotrt2, nonzerotrt1, nonzerotrt2, totaltrt1, totaltrt2,variance, error.lm,error.pois,error.nb,error.zinb,warn.lm,warn.pois,warn.nb,warn.zinb)
  }
}
}
}
all_data <- data.frame(all_data)
if(mixed_effects == TRUE) {
   if(drop_warnings == TRUE) {
names(all_data) <- c("OTUname","lm.coeff","lm.pval","pois.coeff","pois.pval","nb.coeff","nb.pval","zinb.coeff","zinb.pval","mean_trt1","mean_trt2","mean_diff","wgt_mean_trt1","wgt_mean_trt2","wgt_mean_diff","aic.lm","aic.pois","aic.nb","aic.zinb","bic.lm","bic.pois","bic.nb","bic.zinb","nb.disp","zinb.disp","zerotrt1","zerotrt2","nonzerotrt1","nonzerotrt2","totaltrt1","totaltrt2","variance")
   } else {
names(all_data) <- c("OTUname","lm.coeff","lm.pval","pois.coeff","pois.pval","nb.coeff","nb.pval","zinb.coeff","zinb.pval","mean_trt1","mean_trt2","mean_diff","wgt_mean_trt1","wgt_mean_trt2","wgt_mean_diff","aic.lm","aic.pois","aic.nb","aic.zinb","bic.lm","bic.pois","bic.nb","bic.zinb","nb.disp","zinb.disp","zerotrt1","zerotrt2","nonzerotrt1","nonzerotrt2","totaltrt1","totaltrt2","variance","error.lm","error.pois","error.nb","error.zinb","warn.lm","warn.pois","warn.nb","warn.zinb")
   }
} else {
   if(drop_warnings == TRUE) {
names(all_data) <- c("OTUname","lm.coeff","lm.pval","pois.coeff","pois.pval","nb.coeff","nb.pval","zinb.coeff","zinb.pval","mean_trt1","mean_trt2","mean_diff","aic.lm","aic.pois","aic.nb","aic.zinb","bic.lm","bic.pois","bic.nb","bic.zinb","zerotrt1","zerotrt2","nonzerotrt1","nonzerotrt2","totaltrt1","totaltrt2","variance")
   } else {
names(all_data) <- c("OTUname","lm.coeff","lm.pval","pois.coeff","pois.pval","nb.coeff","nb.pval","zinb.coeff","zinb.pval","mean_trt1","mean_trt2","mean_diff","aic.lm","aic.pois","aic.nb","aic.zinb","bic.lm","bic.pois","bic.nb","bic.zinb","zerotrt1","zerotrt2","nonzerotrt1","nonzerotrt2","totaltrt1","totaltrt2","variance","error.lm","error.pois","error.nb","error.zinb","warn.lm","warn.pois","warn.nb","warn.zinb")
   }
}

save.image(file=paste0(workdir,imagename,".Rdata"))

colnames(taxa.name) <- c("OTUname","taxonomy")
all_data <- merge(all_data, taxa.name, by=c("OTUname"))

attach(all_data)

if(choose_by_BIC == TRUE) {
result <- t(sapply(seq(nrow(all_data)), function(i)  {
   best.mod <- c(1,2,3,4)
   bic <- c(as.vector(bic.lm[i]),as.vector(bic.pois[i]),as.vector(bic.nb[i]),as.vector(bic.zinb[i]))
   bic.best <- ifelse(sum(is.na(bic))!=max(best.mod), best.mod[which(bic == min(as.numeric(bic),na.rm=TRUE))],NA)
   mod.name <- c("LinearMod","Poisson","NegBin","ZI-NegBin")[bic.best]
   coefs <- c(as.vector(lm.coeff[i]),as.vector(pois.coeff[i]),as.vector(nb.coeff[i]),as.vector(zinb.coeff[i]))
   best.coef <- tryCatch(coefs[bic.best],error=function(e) NA)
   qvals <- c(as.vector(lm.coeff[i]),as.vector(pois.pval[i]),as.vector(nb.pval[i]),as.vector(zinb.pval[i]))
   best.pval <- tryCatch(qvals[bic.best],error=function(e) NA)
   list(mod.name,as.numeric(as.character(best.coef)),as.numeric(as.character(best.pval)))
}))
} else {
result <- t(sapply(seq(nrow(all_data)), function(i)  {
   best.mod <- c(1,2,3,4)
   aic <- c(as.vector(aic.lm[i]),as.vector(aic.pois[i]),as.vector(aic.nb[i]),as.vector(aic.zinb[i]))
   aic.best <- ifelse(sum(is.na(aic))!=4, best.mod[which(aic == min(as.numeric(aic),na.rm=TRUE))],NA)
   mod.name <- c("LinearMod","Poisson","NegBin","ZI-NegBin")[aic.best]
   coefs <- c(as.vector(aic.coeff[i]),as.vector(pois.coeff[i]),as.vector(nb.coeff[i]),as.vector(zinb.coeff[i]))
   best.coef <- tryCatch(coefs[aic.best],error=function(e) NA)
   qvals <- c(as.vector(aic.pval[i]),as.vector(pois.pval[i]),as.vector(nb.pval[i]),as.vector(zinb.pval[i]))
   best.pval <- tryCatch(qvals[aic.best],error=function(e) NA)
   list(mod.name,as.numeric(as.character(best.coef)),as.numeric(as.character(best.pval)))
}))
}
detach(all_data)
result.mat <- matrix(do.call("rbind",result),ncol=3,nrow=dim(all_data)[1])
colnames(result.mat) <- c("best.mod","best.coef","best.pval")
final.OTUDiff <- data.frame(all_data,result.mat)
final.OTUDiff <- final.OTUDiff[!is.na(final.OTUDiff$best.pval),]
final.OTUDiff$qval.best <- p.adjust(as.vector(final.OTUDiff$best.pval), method="fdr")

setwd(workdir)
write.csv(final.OTUDiff,outputname,row.names=FALSE)


