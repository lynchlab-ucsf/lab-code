## New Many-Model Script:

pacman::p_load(gdata, pscl, stringr, MASS, tidyverse, foreach, phyloseq, parallel, tweedie, glmmTMB)
## These are a set of functions to be used for the "many-model" script (currently supporting Poisson, Negative Binomial, Zero-Inflated Negative Binomial and Tweedie/CPLM, but planning to expand)
## Example usage is : 
##    source("ManyModelScript_July2020.R")
##    final_results <- many_model_script(myotu, mydata2, sampleid=0, model="~ Case.or.Control.Status.Full.Cohort")

## You can either provide a phyloseq object through phy=phyobj, or OTU table + sample data as in the above example.
## Filtering is done and described. You can also change the filtering parameters if your data needs different parameters.

## Parameters in the function:
### phy    - your phyloseq object
### otu    - OTU table if not providing a phyloseq object
### data   - dataset if not providing a phyloseq object
### model  - your statistical model which must start with a curly line. First variable will always be considered your main effect by default
### sampleid - a variable that matches your sample names to your data (if you didn't provide a phyloseq object)
### subjectid - a variable that identifies subjects (this is how repeated measures analysis is initiated)
### pct_pres - The percent prevalence threshold for your analysis (a value is picked by default, and you can change the value here)
### log_reads - The read-depth of OTUs, log transformed, which act as a read depth threshold (a value is picked by default, but it can be changed here)
### main   - The main effect of your analysis if it is not the first variable in your model
### cores  - The number of computational cores to use for your analysis (default=1)
### ref    - Reference group for your analysis (for instance, if you have Cases and Controls, R might choose Cases as your Reference group, but really you want controls because you want positive numbers to be increased in cases)
### run_models - Models to analyze (options include: "pois","negbin","zinfl","tweedie" for Poisson, Negative Binomial, Zero-Inflated Negative Binomial, and Tweedie, respectively). You can also pick and choose these.


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
    if(is.null(phy)) print("Either your OTU table, sample data, or sample ID descriptor are missing or incorrect.")
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
    print("Red lines indicate median of proportion of samples the taxon/feature is found in, and the median log of the total reads of the taxon/feature.")
    print("By default, the upper-right quadrant will be kept. Alternative filtering can be used with 'pct_pres' and 'log_reads' options.")
    droplist <- names(pct_pres_init)[!pct_pres_init>median(pct_pres_init) | !log(reads)>median(log(reads))]
    filt_data <- tokeep$all_data[, !names(tokeep$all_data) %in% droplist]
    filt_taxa <- tokeep$taxa_list[!tokeep$taxa_list %in% droplist]
    print(paste(length(droplist), "of", length(tokeep$taxa_list), "OTUs were filtered"))
    return(list(all_data=filt_data, taxa_list=filt_taxa))
  } else {
    print("At least one filtering level has been chosen manually. Figure represents modified parameters.")
    newpres <- ifelse(!is.null(pct_pres), pct_pres, median(pct_pres_init))
    newread <- ifelse(!is.null(log_reads), log_reads, median(log(reads)))
    plot(log(reads), pct_pres_init)
    abline(h=median(pct_pres), v=median(log_reads),col="red")
    droplist <- names(pct_pres_init)[!pct_pres_init>newpres | !log(reads)>newread]
    filt_data <- tokeep$all_data[, !names(tokeep$all_data) %in% droplist]
    filt_taxa <- tokeep$taxa_list[!tokeep$taxa_list %in% droplist]
    print(paste(length(droplist), "of", length(tokeep$taxa_list), "OTUs were filtered"))
    return(list(all_data=filt_data, taxa_list=filt_taxa))
  }
}

#anly_data <- filter_params(tokeepP)

#model <- "~ Case.or.Control.Status.Full.Cohort"

variable_info <- function(data, model, main, ref=ref) {
  print(paste0("Your model is: ", model, " and your variable of interest is: ", main))
  init_ref <- levels(anly_data$all_data[,main])[1]
  if(!is.na(ref) & !init_ref == ref) {
    anly_data$all_data[,main] <- relevel(anly_data$all_data[,main], ref=ref)
    print(paste0("The reference group changed and is now: ", ref))
  } else {
    print(paste0("Your current reference group is: ", init_ref, "; if you would like to change this, use the 'ref' parameter."))
  }
}


all_models <- function(anly_data, model, feature, run_models=run_models, main=NULL) {
  if("pois" %in% run_models) { pois <- tryCatch(glm(as.formula(paste0(feature, model)), data=anly_data$all_data, family="poisson"), error=function(e) NULL, warning=function(w) NULL)} else {pois <- NA}
  if("negbin" %in% run_models) { negbin <- tryCatch(MASS::glm.nb(as.formula(paste0(feature, model)), data=anly_data$all_data), error=function(e) NULL, warning=function(w) NULL)} else {negbin <- NA}
  if("zinfl" %in% run_models) { zinfl <- tryCatch(pscl::zeroinfl(as.formula(paste0(feature, model, "| 1")), data=anly_data$all_data, dist="negbin"), error=function(e) NULL, warning=function(w) NULL)} else {zinfl <- NA}
  if("tweedie" %in% run_models) { tweedie <- tryCatch(glm(as.formula(paste0(feature, model)), data=anly_data$all_data, family=tweedie(var.power=1.5)), error=function(e) NULL, warning=function(w) NULL) } else {tweedie <- NA}
  models <- sapply(run_models, function(x) ifelse(is.null(get(x)), NA, ifelse(is.na(AIC(get(x))), AICtweedie(get(x)), AIC(get(x)))))
  gomod <- tryCatch(which.min(models))
  if(length(gomod)>0) {
    if(is.null(main)) main <- stringr::str_split(model, " +")[[1]][2]
    if(names(gomod) %in% "zinfl") return(c(OTU=feature, win.model="zinfl",
                                           results=if(length(levels(anly_data$all_data[,main]))>2) {
                                             gdata::unmatrix(summary(zinfl)$coef$count[grepl(main, rownames(summary(zinfl)$coef$count)),], byrow=TRUE)
                                             } else {
                                               summary(zinfl)$coef$count[grepl(main, rownames(summary(zinfl)$coef$count)),]
                                              },
                                           AIC=unlist(models)))
    if(names(gomod) %in% "negbin") return(c(OTU=feature, win.model="negbin",
                                            results=if(length(levels(anly_data$all_data[,main]))>2) {
                                              gdata::unmatrix(summary(negbin)$coef[grepl(main, rownames(summary(negbin)$coef)),], byrow=TRUE)
                                            } else {
                                              summary(negbin)$coef[grepl(main, rownames(summary(negbin)$coef)),]
                                            },
                                            AIC=unlist(models)))
    if(names(gomod) %in% "pois") return(c(OTU=feature, win.model="pois", 
                                          results=if(length(levels(anly_data$all_data[,main]))>2) {
                                            gdata::unmatrix(summary(pois)$coef[grepl(main, rownames(summary(pois)$coef)),], byrow=TRUE)
                                          } else {
                                            summary(pois)$coef[grepl(main, rownames(summary(pois)$coef)),]
                                          },
                                          AIC=unlist(models)))
  if(names(gomod) %in% "tweedie") return(c(OTU=feature, win.model="tweedie", 
                                        results=if(length(levels(anly_data$all_data[,main]))>2) {
                                          gdata::unmatrix(summary(tweedie)$coef[grepl(main, rownames(summary(tweedie)$coef)),], byrow=TRUE)
                                        } else {
                                          summary(tweedie)$coef[grepl(main, rownames(summary(tweedie)$coef)),]
                                        },
                                        AIC=unlist(models)))
  }
}

#### This would contain mixed-effects models
mixed_models <- function(anly_data, model, feature, run_models=run_models, main=NULL, subjid) {
  if("pois" %in% run_models) { pois <- tryCatch(glmmTMB(as.formula(paste0(feature, model, "+ (1|", subjid, ")")), data=anly_data$all_data, family="poisson"), error=function(e) NULL, warning=function(w) NULL)} else {pois <- NA}
  if("negbin" %in% run_models) { negbin <- tryCatch(glmmTMB(as.formula(paste0(feature, model, "+ (1|", subjid, ")")), data=anly_data$all_data, family="nbinom1"), error=function(e) NULL, warning=function(w) NULL)} else {negbin <- NA}
  if("zinfl" %in% run_models) { zinfl <- tryCatch(glmmTMB(as.formula(paste0(feature, model, "+ (1|", subjid, ")")), data=anly_data$all_data, family="nbinom1", ziformula=~1), error=function(e) NULL, warning=function(w) NULL)} else {zinfl <- NA}
  if("tweedie" %in% run_models) { tweedie <- tryCatch(glmmTMB(as.formula(paste0(feature, model, "+ (1|", subjid, ")")), data=anly_data$all_data, family="tweedie"), error=function(e) NULL, warning=function(w) NULL) } else {tweedie <- NA}
  models <- sapply(run_models, function(x) ifelse(is.null(get(x)), NA, ifelse(is.na(AIC(get(x))), AICtweedie(get(x)), AIC(get(x)))))
  gomod <- tryCatch(which.min(models))
  if(length(gomod)>0) {
    if(is.null(main)) main <- stringr::str_split(model, " +")[[1]][2]
    if(names(gomod) %in% "zinfl") return(c(OTU=feature, win.model="zinfl",
                                           results=if(length(levels(anly_data$all_data[,main]))>2) {
                                             gdata::unmatrix(summary(zinfl)$coef$cond[grepl(main, rownames(summary(zinfl)$coef$cond)),], byrow=TRUE)
                                           } else {
                                             summary(zinfl)$coef$cond[grepl(main, rownames(summary(zinfl)$coef$cond)),]
                                           },
                                           AIC=unlist(models)))
    if(names(gomod) %in% "negbin") return(c(OTU=feature, win.model="negbin",
                                            results=if(length(levels(anly_data$all_data[,main]))>2) {
                                              gdata::unmatrix(summary(negbin)$coef$cond[grepl(main, rownames(summary(negbin)$coef$cond)),], byrow=TRUE)
                                            } else {
                                              summary(negbin)$coef$cond[grepl(main, rownames(summary(negbin)$coef$cond)),]
                                            },
                                            AIC=unlist(models)))
    if(names(gomod) %in% "pois") return(c(OTU=feature, win.model="pois", 
                                          results=if(length(levels(anly_data$all_data[,main]))>2) {
                                            gdata::unmatrix(summary(pois)$coef$cond[grepl(main, rownames(summary(pois)$coef$cond)),], byrow=TRUE)
                                          } else {
                                            summary(pois)$coef$cond[grepl(main, rownames(summary(pois)$coef$cond)),]
                                          },
                                          AIC=unlist(models)))
    if(names(gomod) %in% "tweedie") return(c(OTU=feature, win.model="tweedie", 
                                             results=if(length(levels(anly_data$all_data[,main]))>2) {
                                               gdata::unmatrix(summary(tweedie)$coef$cond[grepl(main, rownames(summary(tweedie)$coef$cond)),], byrow=TRUE)
                                             } else {
                                               summary(tweedie)$coef$cond[grepl(main, rownames(summary(tweedie)$coef$cond)),]
                                             },
                                             AIC=unlist(models)))
  }
}


many_model_script <- function(otu=NULL, data=NULL, phy=NULL, sampleid=NULL, subjectid=NULL, pct_pres=NULL, log_reads=NULL, model, main=NULL, cores=1, ref=NA, run_models=c("pois","negbin","zinfl","tweedie")) {
  if(is.null(main)) main <- stringr::str_split(model, " +")[[1]][2]
  anly_data <- make_data(otu, data, phy, sampleid)
  filt_data <- filter_params(anly_data, pct_pres, log_reads)
  variable_info(filt_data, model, main, ref)
  if(is.null(subjectid)) { print("Please provide a subjectID before continuing.") ; stop() }
  if(sum(duplicated(filt_data$all_data[,subjectid])) == 0) {
  test <- mcmapply(all_models, feature=filt_data$taxa_list, MoreArgs = list(model=model, anly_data=filt_data, main=main, run_models=run_models), mc.cores=cores)
  } else {
  test <- mcmapply(mixed_models, feature=filt_data$taxa_list, MoreArgs = list(model=model, anly_data=filt_data, main=main, run_models=run_models, subjid=subjectid), mc.cores=cores)
  }
  taxon_results <- data.frame(t(data.frame(test))) %>%
  rename_at(vars(contains(main)), function(x) sub(main,"", x)) %>%
   mutate_at(vars(contains("Pr...")), list(p.fdr=function(x) p.adjust(as.numeric(as.character(x)), method="fdr"))) %>%
   mutate_all(unlist) %>% 
    print()
if(!is.null(phy)) taxon_results <- merge(taxon_results, tax_table(phy)@.Data, by.x="OTU",by.y=0)
return(taxon_results)
}

#myphy2 <- subset_samples(myphy, !Case.or.Control.Status.Full.Cohort %in% "")
#my_results <- many_model_script(phy=myphy2, model="~ Case.or.Control.Status.Full.Cohort", cores=10, ref="Control", run_models=c("pois","negbin","zinfl","tweedie"), subjectid = "Subject.Identifier.for.the.Study")
