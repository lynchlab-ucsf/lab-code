## New Many-Model Script:

pacman::p_load(gdata, pscl, stringr, MASS, tidyverse, foreach, phyloseq, parallel, tweedie, glmmTMB, cplm, pbmcapply)
## These are a set of functions to be used for the "many-model" script (currently supporting Linear Model, Compound Poisson Linear Model, Poisson, Negative Binomial, Zero-Inflated Negative Binomial and Tweedie)
## Example usage is : 
##    source("ManyModelScript_July2020.R")
##    final_results <- many_model_script(myotu, mydata2, sampleid=0, model="~ Case.or.Control.Status.Full.Cohort")

## You can either provide a phyloseq object through phy=phyobj, or OTU table + sample data as in the above example.
## Filtering is done and described. You can also change the filtering parameters if your data needs different parameters.

## Parameters in the function:
### phy        - your phyloseq object
### otu        - OTU table if not providing a phyloseq object
### data       - dataset if not providing a phyloseq object
### model      - your statistical model which must start with a curly line followed by a space. First variable will always be considered your main effect by default
### sampleid   - a variable that matches your sample names to your data (if you didn't provide a phyloseq object)
### subjectid  - a variable that identifies subjects (this is how repeated measures analysis is initiated)
### pct_pres   - The percent prevalence threshold for your analysis (a value is picked by default, and you can change the value here)
### log_reads  - The read-depth of OTUs, log transformed, which act as a read depth threshold (a value is picked by default, but it can be changed here)
### main       - The main effect of your analysis if it is not the first variable in your model
### cores      - The number of computational cores to use for your analysis (default=1)
### ref        - Reference group for your analysis (for instance, if you have Cases and Controls, R might choose Cases as your Reference group, but really you want controls because you want positive numbers to be increased in cases)
### run_models - Models to analyze (options include: "lm","cplm","pois","negbin","zinfl","tweedie" for Linear Model, Compund Poisson Linear Model, Poisson, Negative Binomial, Zero-Inflated Negative Binomial, and Tweedie, respectively). You can also pick and choose these.


## Data imported for testing (Keeping for troubleshooting purposes)
#myotu <- read.table("/data/Users/kmccauley/MUPPITS/OTUtables/MUPPITS_OTUtable_initial.txt", header=T, check.names=F, sep="\t",comment="", row.names=1, skip=1)
#mydata <- read.csv("/data/Users/kmccauley/MUPPITS/Analysis_Data/MUPPITS_Merged_Mapping_UCSF_withDomGen.csv", row.names=1)
#myphy <- readRDS("/data/Users/kmccauley/MUPPITS/OTUtables/usearch_nonrare_phyloseq_amend.RData")
#mydata2 <- mydata[!mydata$Case.or.Control.Status.Full.Cohort %in% "", ]

print("UPDATE as of 5/6: The output of this script has been significantly simplified. Your previous plotting script will no longer work. Review the output before plotting.")
print("UPDATE as of 5/6: Poisson has been removed from the list of models by default (and your results may change). If you would like to force it, specify run_models=run_models=c('lm','pois','negbin','zinfl','tweedie'), but this is not recommended.")
print("UPDATE as of 5/6: This script should theoretically work for metabolomics data too -- I think I've fixed the 'bug'. If you encounter issues, please submit a bug report.")

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
    if(taxa_are_rows(phy)) all_data <- merge(sample_data(phy), t(otu_table(phy)), by=0) else all_data <- merge(sample_data(phy), otu_table(phy), by=0)
    taxa_list <- taxa_names(phy)
  }
  container <- list(all_data=all_data, taxa_list=taxa_list)
  return(container)
}

#tokeep <- make_data(myotu, mydata2, sampleid=0)
#tokeep <- make_data(phy=myphy)

## Will help determine what filtering parameters to use on the data.
filter_params <- function(tokeep, pct_pres=NULL, log_reads=NULL) {
  reads <- colSums(tokeep$all_data[,tokeep$taxa_list])
  pct_pres_init <- colSums(tokeep$all_data[,tokeep$taxa_list]>0)/nrow(tokeep$all_data[,tokeep$taxa_list])
  tokeep$all_data$total_reads <- rowSums(tokeep$all_data[,tokeep$taxa_list])
  if(sum(pct_pres_init == 0)>1) { print("You have at least one taxon not present in any samples. Please filter your OTU table to present taxa before continuing. \n In Phyloseq: filter_taxa(phy, function(x) sum(x) > 0, TRUE)"); stop() }
  
  if(is.null(pct_pres) & is.null(log_reads)) {
    plot(log(reads), pct_pres_init)
    abline(h=median(pct_pres_init), v=median(log(reads)),col="red")
    print(paste0("Red lines indicate median of proportion of samples the taxon/feature is found in (", round(median(pct_pres_init),3), "), and the median log of the total reads (", round(median(log(reads)),3),") of the taxon/feature."))
    if(median(pct_pres_init) < 0.05) print("THE MEDIAN PREVALENCE OF YOUR TAXA MAY BE TOO LOW FOR DEFAULTS. If you receive an error after models are complete, consider increasing 'pct.pres'.")
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
#anly_data <- filter_params(tokeep)

#model <- "~ Case.or.Control.Status.Full.Cohort"

variable_info <- function(data, model, main, ref=ref) {
  print(paste0("Your model is: ", model, " and your variable of interest is: ", main))
  if(is.factor(data$all_data[,main])) {
  init_ref <- levels(data$all_data[,main])[1]
  print(table(data$all_data[,main]))
  if(!is.na(ref) & !init_ref == ref) {
    data$all_data[,main] <- relevel(data$all_data[,main], ref=ref)
    print(paste0("The reference group changed and is now: ", ref))
  } else {
    print(paste0("Your current reference group is: ", init_ref, "; if you would like to change this, use the 'ref' parameter."))
  }
  }
  return(data)
}


all_models <- function(anly_data, model, feature, run_models=run_models, main=NULL) {
  adj <- ifelse(length(unique(anly_data$all_data$total_reads))>1, "+ total_reads", "")
  if("lm" %in% run_models) { lm <- tryCatch(glm(as.formula(paste0(feature, model, adj)), data=anly_data$all_data, family="gaussian"), error=function(e) NULL, warning=function(w) NULL)} else {lm <- NA}
  if("cplm" %in% run_models) { cplm <- tryCatch(cpglm(as.formula(paste0(feature, model, adj)), data=anly_data$all_data), error=function(e) NULL, warning=function(w) NULL)} else {cplm <- NA}
  if("pois" %in% run_models) { pois <- tryCatch(glm(as.formula(paste0(feature, model, adj)), data=anly_data$all_data, family="poisson"), error=function(e) NULL, warning=function(w) NULL)} else {pois <- NA}
  if("negbin" %in% run_models) { negbin <- tryCatch(MASS::glm.nb(as.formula(paste0(feature, model, adj)), data=anly_data$all_data), error=function(e) NULL, warning=function(w) NULL)} else {negbin <- NA}
  if("zinfl" %in% run_models) { zinfl <- tryCatch(pscl::zeroinfl(as.formula(paste0(feature, model, adj, "| 1")), data=anly_data$all_data, dist="negbin"), error=function(e) NULL, warning=function(w) NULL)} else {zinfl <- NA}
  if("tweedie" %in% run_models) { tweedie <- tryCatch(glm(as.formula(paste0(feature, model, adj)), data=anly_data$all_data, family=statmod::tweedie(var.power=1.5)), error=function(e) NULL, warning=function(w) NULL) } else {tweedie <- NA}
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
    if(names(gomod) %in% "lm") return(c(OTU=feature, win.model="lm", 
                                             results=if(length(levels(anly_data$all_data[,main]))>2) {
                                               gdata::unmatrix(summary(lm)$coef[grepl(main, rownames(summary(lm)$coef)),], byrow=TRUE)
                                             } else {
                                               summary(lm)$coef[grepl(main, rownames(summary(lm)$coef)),]
                                             },
                                             AIC=unlist(models)))
    if(names(gomod) %in% "cplm") return(c(OTU=feature, win.model="cplm", 
                                             results=if(length(levels(anly_data$all_data[,main]))>2) {
                                               out <- capture.output(cplm <- invisible(summary(cplm)$coef))
                                               gdata::unmatrix(cplm[grepl(main, rownames(cplm)),], byrow=TRUE)
                                             } else {
                                               out <- capture.output(cplm <- invisible(summary(cplm)$coef))
                                               cplm[grepl(main, rownames(cplm)),]
                                             },
                                             AIC=unlist(models)))
  }
}

#### This would contain mixed-effects models
mixed_models <- function(anly_data, model, feature, run_models=run_models, main=NULL, subjid) {
  adj <- ifelse(length(unique(anly_data$all_data$total_reads))>1, "+ total_reads", "")
  if("lm" %in% run_models) { lm <- tryCatch(suppressWarnings(glmmTMB(as.formula(paste0(feature, model, adj, "+ (1|", subjid, ")")), data=anly_data$all_data, family="gaussian")), error=function(e) NULL)} else {lm <- NA}
  if("cplm" %in% run_models) { cplm <- tryCatch(suppressWarnings(glmmTMB(as.formula(paste0(feature, model, adj, "+ (1|", subjid, ")")), data=anly_data$all_data, family="compois")), error=function(e) NULL)} else {cplm <- NA}
  if("pois" %in% run_models) { pois <- tryCatch(suppressWarnings(glmmTMB(as.formula(paste0(feature, model, adj, "+ (1|", subjid, ")")), data=anly_data$all_data, family="poisson")), error=function(e) NULL)} else {pois <- NA}
  if("negbin" %in% run_models) { negbin <- tryCatch(suppressWarnings(glmmTMB(as.formula(paste0(feature, model, adj, "+ (1|", subjid, ")")), data=anly_data$all_data, family="nbinom1")), error=function(e) NULL)} else {negbin <- NA}
  if("zinfl" %in% run_models) { zinfl <- tryCatch(suppressWarnings(glmmTMB(as.formula(paste0(feature, model, adj, "+ (1|", subjid, ")")), data=anly_data$all_data, family="nbinom1", ziformula=~1)), error=function(e) NULL)} else {zinfl <- NA}
  if("tweedie" %in% run_models) { tweedie <- tryCatch(suppressWarnings(glmmTMB(as.formula(paste0(feature, model, adj, "+ (1|", subjid, ")")), data=anly_data$all_data, family="tweedie")), error=function(e) NULL) } else {tweedie <- NA}
  models <- sapply(run_models, function(x) ifelse(is.null(get(x)), NA, ifelse(is.na(AIC(get(x))), AIC(get(x)), AIC(get(x)))))
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
    if(names(gomod) %in% "lm") return(c(OTU=feature, win.model="lm", 
                                             results=if(length(levels(anly_data$all_data[,main]))>2) {
                                               gdata::unmatrix(summary(lm)$coef$cond[grepl(main, rownames(summary(lm)$coef$cond)),], byrow=TRUE)
                                             } else {
                                               summary(lm)$coef$cond[grepl(main, rownames(summary(lm)$coef$cond)),]
                                             },
                                             AIC=unlist(models)))
    if(names(gomod) %in% "cplm") return(c(OTU=feature, win.model="cplm", 
                                             results=if(length(levels(anly_data$all_data[,main]))>2) {
                                               gdata::unmatrix(summary(cplm)$coef$cond[grepl(main, rownames(summary(cplm)$coef$cond)),], byrow=TRUE)
                                             } else {
                                               summary(cplm)$coef$cond[grepl(main, rownames(summary(cplm)$coef$cond)),]
                                             },
                                             AIC=unlist(models)))
  }
}


many_model_script <- function(otu=NULL, data=NULL, phy=NULL, sampleid=NULL, subjectid=NULL, pct_pres=NULL, log_reads=NULL, model, main=NULL, cores=1, ref=NA, run_models=c("lm","negbin","zinfl","tweedie")) {
  if(is.null(main)) main <- stringr::str_split(model, " +")[[1]][2]
  anly_data <- make_data(otu, data, phy, sampleid)
  filt_data <- filter_params(anly_data, pct_pres, log_reads)
  filt_data <- variable_info(filt_data, model, main=main, ref)
  ref <- levels(filt_data$all_data[,main])[1]
  if(is.null(subjectid)) { print("Please provide a subjectID before continuing.") ; stop() }
  if(sum(duplicated(filt_data$all_data[,subjectid])) == 0) {
    print("Running Fixed-Effects Models")
    stat_anly <- pbmcmapply(all_models, feature=filt_data$taxa_list, MoreArgs = list(model=model, anly_data=filt_data, main=main, run_models=run_models), mc.cores=cores)
    stat_anly <- if(any(class(stat_anly) %in% "list")) stat_anly$value else stat_anly
    print("Models finished calculating...")
  } else {
    print("Running Mixed-Effects Models")
    stat_anly <- pbmcmapply(mixed_models, feature=filt_data$taxa_list, MoreArgs = list(model=model, anly_data=filt_data, main=main, run_models=run_models, subjid=subjectid), mc.cores=cores)
    stat_anly <- if(any(class(stat_anly) %in% "list")) stat_anly$value else stat_anly
  }
  if(is.factor(filt_data$all_data[,main])) {
    differences <- filt_data$all_data %>% 
      group_by(.data[[main]]) %>% 
      summarize_at(vars(starts_with(filt_data$taxa_list)), mean) %>% 
      column_to_rownames(var=main) %>%
      t() %>% data.frame(check.names = F) %>% 
      rownames_to_column("OTU") %>% 
      mutate_at(vars(-starts_with(c(!!ref, "OTU"))), list(mean_diff = ~ . - .data[[ref]]))
    taxon_results <- data.frame(t(data.frame(stat_anly))) %>%
      rename_at(vars(contains(main)), function(x) sub(main,"", x)) %>%
      dplyr::select(!contains(c("Std..Error",".value"))) %>%
      rename_with(~ gsub("Pr......","Pvalue", .x)) %>% 
      rename_with(~ gsub("results.","", .x)) %>% 
      mutate_at(vars(contains("results")), list(function(x) as.numeric(as.character(x)))) %>%
      mutate_at(vars(contains("Pvalue")), list(p.fdr=function(x) p.adjust(x, method="fdr"))) %>%
      mutate_at("OTU", as.character) %>% 
      mutate_all(unlist) %>% 
      left_join(differences, by="OTU")
  } else {
    taxon_results <- data.frame(t(data.frame(stat_anly))) %>%
      rename_at(vars(contains(main)), function(x) sub(main,"", x)) %>%
      dplyr::select(!contains(c("Std..Error",".value"))) %>%
      rename_with(~ gsub("Pr......","Pvalue", .x)) %>% 
      rename_with(~ gsub("results.","", .x)) %>% 
      mutate_at(vars(contains("results")), list(function(x) as.numeric(as.character(x)))) %>%
      mutate_at(vars(contains("Pvalue")), list(p.fdr=function(x) p.adjust(x, method="fdr"))) %>%
      mutate_all(unlist)
  }
  if(!is.null(phy) & !is.null(tax_table(phy, errorIfNULL=FALSE))) taxon_results <- merge(taxon_results, tax_table(phy)@.Data, by.x="OTU",by.y=0)
  if(is.null(phy) & is.null(otu$taxonomy)) {
    otu.tax <- cbind(otu=rownames(otu), taxonomy=as.character(otu$taxonomy))
    taxon_results <- merge(taxon_results, otu.tax, by.x="OTU", by.y="otu")
  }
  return(taxon_results)
}

#myphy2 <- subset_samples(myphy, !Case.or.Control.Status.Full.Cohort %in% "" & Analysis.Visit %in% "Visit 1a")
#myphy3 <- filter_taxa(myphy2, function(x) sum(x) > 0, TRUE)
#my_results <- many_model_script(phy=myphy3, model="~ Case.or.Control.Status.Full.Cohort + Age.in.years", cores=10, ref="Control", run_models=c("lm","cplm","pois","negbin","zinfl","tweedie"), subjectid = "Subject.Identifier.for.the.Study")
