## New Many-Model Script:

library(phyloseq)
library(foreach)
library(doParallel)
library(tidyverse)
library(MASS)
library(stringr)
library(pscl)

## These are a set of functions to be used for the "many-model" script (currently supporting Poisson, Negative Binomial and Zero-Inflated Negative Binomial, but planning to expand)
## Example usage is : 
##    source("ManyModelScript_July2020.R")
##    final_results <- many_model_script(myotu, mydata2, sampleid=0, model="~ Case.or.Control.Status.Full.Cohort")

## You can either provide a phyloseq object through phy=phyobj, or OTU table + sample data as in the above example.
## Filtering is done and described. You can also change the filtering parameters if your data needs different parameters.


## Data imported for testing (Keeping for troubleshooting purposes)
#myotu <- read.table("/data/Users/kmccauley/MUPPITS/OTUtables/MUPPITS_OTUtable_initial.txt", header=T, check.names=F, sep="\t",comment="", row.names=1, skip=1)
#mydata <- read.csv("/data/Users/kmccauley/MUPPITS/Analysis_Data/MUPPITS_Merged_Mapping_UCSF_withDomGen.csv", row.names=1)
#myphy <- readRDS("/data/Users/kmccauley/MUPPITS/OTUtables/usearch_nonrare_phyloseq_amend.RData")
#mydata2 <- mydata[!mydata$Case.or.Control.Status.Full.Cohort %in% "", ]

make_data <- function(otu=NULL, data=NULL, phy=NULL, sampleid=NULL) {
  ## If a phyloseq object isn't provided, but an OTU table and sample data *is* provided.
  if(is.null(phy) & !is.null(otu) & !is.null(data) & !is.null(sampleid)) {
    all_data <- merge(t(otu), data, by.x=0, by.y=sampleid)
    print(paste("Before merging, your sample data had", dim(data)[1], "samples and", dim(data)[2], "variables, while your OTU table had", dim(otu)[1], "OTUs and", dim(otu)[2],"samples."))
    print(paste("After merging, you have", dim(all_data)[1], "samples."))
    taxa_list <- rownames(otu)
    all_data[,taxa_list] <- lapply(all_data[,taxa_list], function(x) as.numeric(as.character(x)))
  } else {
    if(is.null(phy)) print("Either your OTU table, sample data, or sample ID descriptor are missing or incorrect")
  }
  if(!is.null(phy)) {
    print(phy)
    all_data <- merge(sample_data(phy), t(otu_table(phy)), by=0)
    taxa_list <- taxa_names(phy)
  }
  container <- list(all_data=all_data, taxa_list=taxa_list)
  return(container)
}

#tokeep <- make_data(myotu, mydata2, sampleid=0)
#tokeepP <- make_data(phy=myphy)

## Will help determine what filtering parameters to use on the data.
filter_params <- function(tokeep, pct_pres=NULL, log_reads=NULL) {
  reads <- colSums(tokeep$all_data[,tokeep$taxa_list])
  pct_pres_init <- colSums(tokeep$all_data[,tokeep$taxa_list]>0)/nrow(tokeep$all_data[,tokeep$taxa_list])
  
  if(is.null(pct_pres) & is.null(log_reads)) {
    plot(log(reads), pct_pres_init)
    abline(h=median(pct_pres_init), v=median(log(reads)),col="red")
    print("Red lines indicate median of proportion of samples the sequence is found in, and the median log of the total reads.")
    print("By default, the upper-right quadrant will be kept. Alternative filtering can be used with 'pct_pres' and 'log_reads' options.")
    droplist <- names(pct_pres_init)[!pct_pres_init>median(pct_pres_init) & !log(reads)>median(log(reads))]
    filt_data <- tokeep$all_data[, !names(tokeep$all_data) %in% droplist]
    filt_taxa <- tokeep$taxa_list[!tokeep$taxa_list %in% droplist]
    print(paste(length(droplist), "of", length(tokeep$taxa_list), "OTUs were filtered"))
    return(list(all_data=filt_data, taxa_list=filt_taxa))
  } else {
    print("At least one filtering level has been chosen manually.")
    newpres <- ifelse(!is.null(pct_pres), pct_pres, median(pct_pres_init))
    newread <- ifelse(!is.null(log_reads), log_reads, median(log(reads)))
    droplist <- names(pct_pres_init)[!pct_pres_init>newpres & !log(reads)>newread]
    filt_data <- tokeep$all_data[, !names(tokeep$all_data) %in% droplist]
    filt_taxa <- tokeep$taxa_list[!tokeep$taxa_list %in% droplist]
    print(paste(length(droplist), "of", length(tokeep$taxa_list), "OTUs were filtered"))
    return(list(all_data=filt_data, taxa_list=filt_taxa))
  }
}

#anly_data <- filter_params(tokeep)

#written_formula <- "~ Case.or.Control.Status.Full.Cohort"

all_models <- function(anly_data, model, feature, main=NULL) {
  pois <- tryCatch(glm(as.formula(paste0(feature, model)), data=anly_data$all_data, family="poisson"), error=function(e) NULL)
  negbin <- tryCatch(MASS::glm.nb(as.formula(paste0(feature, model)), data=anly_data$all_data), error=function(e) NULL)
  zinfl <- tryCatch(pscl::zeroinfl(as.formula(paste0(feature, model, "| 1")), data=anly_data$all_data, dist="negbin"), error=function(e) NULL)
  models <- list(pois=pois, negbin=negbin, zinfl=zinfl)
  gomod <- tryCatch(which.min(flatten(lapply(models, AIC))))
  if(length(gomod)>0) {
    if(is.null(main)) main <- stringr::str_split(model, " +")[[1]][2]
    if(names(gomod) %in% "zinfl") return(c(OTU=feature, win.model="zinfl", 
                                           results=summary(zinfl)$coef$count[grepl(main, rownames(summary(zinfl)$coef$count)),],
                                           AIC=unlist(lapply(models, function(x) ifelse(is.null(AIC(x)), NA, AIC(x))))))
    if(names(gomod) %in% "negbin") return(c(OTU=feature, win.model="negbin",
                                            results=summary(negbin)$coef[grepl(main, rownames(summary(negbin)$coef)),],
                                            AIC=unlist(lapply(models, function(x) ifelse(is.null(AIC(x)), NA, AIC(x))))))
    if(names(gomod) %in% "pois") return(c(OTU=feature, win.model="pois", 
                                          results=summary(pois)$coef[grepl(main, rownames(summary(pois)$coef)),],
                                          AIC=unlist(lapply(models, function(x) ifelse(is.null(AIC(x)), NA, AIC(x))))))
  }
}

many_model_script <- function(otu=NULL, data=NULL, phy=NULL, sampleid=NULL, pct_pres=NULL, log_reads=NULL, model, main=NULL) {
  anly_data <- make_data(otu, data, phy, sampleid)
  filt_data <- filter_params(anly_data, pct_pres, log_reads)
registerDoParallel(cores = 10)
test <- foreach(i=filt_data$taxa_list, .combine=rbind) %dopar%
  all_models(filt_data, model, feature=i)

taxon_results <- data.frame(test) %>% 
   mutate(p.fdr = p.adjust(as.numeric(as.character(results.Pr...z..)), method="fdr")) %>%
   arrange(p.fdr) %>% 
   mutate_all(unlist)
return(taxon_results)
}


