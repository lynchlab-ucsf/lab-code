## Alpha Diversity Script

## This script will perform alpha diversity analyses on your sample data.
### First argument is your phyloseq object,
### Second argument is the subject ID in the dataset
### Third argument asks whether you want to use specific variables in your analysis 
      ### If you don't, no need to specify variables -- all samples in your sample data will be used
      ### If you do, then you will need to specify the variables you want run. An example is below.

#### Example Usage:
## source("/data/Users/kmccauley/LabCode/AlphaDiversity_script.R")
## (read in phyloseq object)
## phy_tree(phy_filt_tree) <- ape::root(phy_tree(phy_filt_tree), 1, resolve.root=TRUE)
## alpha_diversity_analysis(phy_filt_tree, subjectid = "subject", use_all_variables=c("ethnicity","asthmastatus","ige_values"))

## Steps: input would be phyloseq object
alpha_diversity_analysis <- function(phy, subjectid, specified_variables=TRUE) {
  diversities <- c("Faith's Diversity"="PD_whole_tree","Richness"="chao1","Evenness"="equitability")
  require(vegan)
  require(dplyr)
  require(ggplot2)
  require(lmerTest)
  
  if(specified_variables == TRUE) vars <- names(phy@sam_data)
  if(specified_variables != TRUE) vars <- specified_variables

    ## Generate Alpha Diversity Values
    sample_data(phy)$equitability <- diversity(t(otu_table(phy)))/log(specnumber(t(otu_table(phy))))
    sample_data(phy)$chao1 <- t(estimateR(t(otu_table(phy))))[,2]
    sample_data(phy)$PD_whole_tree <- picante::pd(t(otu_table(phy)), phy_tree(phy))[,2]
    
    mydata <- data.frame(phy@sam_data)
    cont.alpha.div <- function(myvar,alpha,mydata) {
      z <- round(summary(lm(as.numeric(mydata[,alpha]) ~ as.numeric(mydata[,myvar])))$coef[-1,c(1,4)], 3)
      return(z)
    }
    rm.alpha.div <- function(myvar,alpha,mydata, subjid=subjectid) { ## Edited for repeated measures
      z <- round(summary(lmer(as.numeric(mydata[,alpha]) ~ as.numeric(mydata[,myvar]) + (1|mydata[,subjid]), data=mydata))$coef[-1,c(1,4)], 3)
      return(z)
    }
    vars.to.plot <- list()
    
    alpha_plots <- function(y, x, mydata) {
      ggplot(mydata, aes_string(x=x, y=y)) + 
        geom_point() +
        geom_smooth(method="glm") +
        ylab(x) +
        ggtitle(y) +
        xlab(y)
    }
    for (i in diversities) {
      alpha.name <- names(diversities[diversities %in% i])
      
      
      if(sum(duplicated(mydata[,subjectid]))==0) init.obj <- do.call(rbind, sapply(vars, cont.alpha.div, alpha = i, mydata=mydata)) %>% 
        data.frame(check.names = FALSE) %>% 
        tibble::rownames_to_column("Variable") %>% 
        arrange(`Pr(>|t|)`) %>% 
        filter(!is.na(`Pr(>|t|)`))
      if(sum(duplicated(mydata[,subjectid]))>0) init.obj <- do.call(rbind, sapply(vars, rm.alpha.div, alpha = i, mydata=mydata)) %>% 
          data.frame(check.names = FALSE) %>% 
          tibble::rownames_to_column("Variable") %>% 
          arrange(`Pr(>|t|)`) %>% 
          filter(!is.na(`Pr(>|t|)`))
      
      vars.to.plot <- init.obj$Variable[init.obj[, "Pr(>|t|)"] < 0.05]
      
      myHeader <- c(alpha.name = 3)
      names(myHeader) <- alpha.name
      
      init.obj %>% 
        mutate(`Pr(>|t|)` = cell_spec(`Pr(>|t|)`, "html", color = ifelse(`Pr(>|t|)` <  0.05, "red", "grey80"))) %>%
        kable(format = "html", escape = F) %>% 
        kable_styling("striped", full_width = F) %>% 
        add_header_above(header = myHeader) %>% 
        print()
      
      for (j in vars.to.plot) {
        alpha_var_plots <- alpha_plots(i, j, mydata=mydata)
        print(alpha_var_plots)
        ggsave(paste0("alpha_vars_", i, "_", j, ".pdf"), alpha_var_plots, device = cairo_pdf, 
               dpi = 300)
      }
    }
}
