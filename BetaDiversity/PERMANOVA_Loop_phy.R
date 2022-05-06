## Create a beta diversity loop script for the lab...
## Written January 6, 2022
## Katie McCauley

## You can run this by *just* providing your phyloseq object (which will run on all variables) or you can provide the specific variables/distance matrices you want the script to use.
## Therefore your code can look like:
#### source("PERMANOVA_Loop_phy.R")
#### beta_diversity_loop(myphy)
## OR
#### source("PERMANOVA_Loop_phy.R")
#### my_variables <- c("SampleGroup1","TestingGroup2")
#### beta_diversity_loop(myphy, dists=c("bray","canberra"), var_list=my_variables)

## It will return a list of PERMANOVA for every variable and every distance matrix. 


var_anly <- function(var, indata, dist_name=NULL, dist_obj=NULL) {
  require(phyloseq)
  print(var)
  if(!is.null(dist_name)) {
    make_dist <- phyloseq::distance(indata, method=dist_name)
  }
  if(!is.null(dist_obj)) {
    make_dist <- dist_obj
  }
  ## This is theoretically the best way to subset a phyloseq object within a function
  samp_to_keep <- rownames(indata@sam_data[!is.na(indata@sam_data[,var]) & !indata@sam_data[,var] %in% ""])
  if(nrow(unique(indata@sam_data[,var]))>1) {
    ss.phy <- prune_samples(samp_to_keep, indata)
    ss.dist <- as.matrix(make_dist)[sample_names(ss.phy), sample_names(ss.phy)]
    permanova <- vegan::adonis2(as.formula(paste0("ss.dist ~ ", var)), data=data.frame(ss.phy@sam_data))
    ## Now I need to make the permanova results look pretty
    return(list(n=length(samp_to_keep), R2=permanova$R2[[1]], P=permanova$`Pr(>F)`[[1]]))
  } else {
    return(list(n=NA, R2=NA, P=NA))
  }
}

beta_diversity_loop <- function(phy, dists=c("bray","wunifrac","canberra","unifrac"), var_list=NULL) {
  ## First, create the distance matrices requested...
  print("Generating All Distance Matrices....")
  all_dists <- sapply(dists, function(x) phyloseq::distance(phy, x), simplify = F)
  
  ## get variable names of interest if they weren't provided
  if(is.null(var_list)) {
    var_list <- names(sample_data(phy))
  }
  all_results <- list()
  for(i in 1:length(all_dists)) {
    print(paste("Running PERMANOVAs for Distance Matrix", i))
    all_results[[dists[i]]] <- t(sapply(var_list, var_anly, indata=phy, dist_obj=all_dists[[i]]))
  }
  return(all_results)
}
## My test inputs
#myphy <- readRDS("/mnt/kmccauley/MUPPITS/MUPPITS_phy_20722.rds")
#beta_diversity_loop(myphy, dists=c("bray","canberra"), var_list=c("PrimerPlate","SampleType"))
