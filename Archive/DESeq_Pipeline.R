#DESeq Script
#Katie McCauley
# Script to run DESeq on an OTU table and map file.
# For any issues, please note them here: https://ucsf.box.com/s/vek039bed2krlwj8714lxmhnzdmu6dcf

# Most recent update: July 17, 2018
## 	Added reference group/control group options and output options
##	Also made a temporary change that makes this code available for use with the installable version of DESeq2 (hope to change this soon!)

# Code based on : https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html

rm(list=ls())

## Install packages:
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
pacman::p_load(DESeq2,phyloseq)


my_map <- read.table("Location/of/Map_File.txt", header=T, check.names=F, comment="", sep="\t")
my_otu_table <- read.table("/Location/of/OTU_table.txt", header=TRUE, check.names=F, comment="", sep="\t",row.names = 1)
var_of_int <- "my_comparison_variable"
#Specify preferred significance level:
alpha <- 0.1
ref_group <- "reference_group" #These can be null, and it will just compare the first against the last level of the variable (which is now evident in the printed results, which include the comparison being made). Eventually, my goal will be to do every possible combination if these are null (but this is not yet implemented).
cont_group <- "contrast_group"
output <- "/Location/to/put/output.csv" # Change to a CSV file name if printed results are desired. Can also be NULL to print to the console.



##### Functions Needed for Later:

split.tax <- function(tax.dat) {
  taxanames <- strsplit(as.character(tax.dat[,2]),"; ")
  mat <- t(sapply(taxanames,
                  function(x,m) c(x,rep(NA,m-length(x))),
                  max(rapply(taxanames,length))))
  
  newnames <- substr(mat,4,35)
  newnames <- as.matrix(newnames)
  colnames(newnames) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  row.names(newnames) <- tax.dat[,1]
  return(newnames)
}

print_res <- function(x, phy, alpha=0.1, sig.only=TRUE, ref= NULL,  cont=NULL, var=var_of_int, file=NULL) {
  if(is.null(ref) & is.null(cont)) {
    # This is where the code would be modified to print all possible comparisons...
    res <- results(x, pAdjustMethod = "fdr")
    varlevs <- levels(sample_data(phy)[[var_of_int]])
    res$comparison <- paste0(varlevs[1], "_vs_", varlevs[length(varlevs)])
  } else {
    res <- results(x, contrast=c(var, cont, ref), pAdjustMethod = "fdr")
    res$comparison <- paste0(ref, "_vs_", cont)
  }
  res = res[order(res$padj, na.last=NA), ]
  if(sig.only==TRUE) {
    sigtab = res[(res$padj < alpha), ]
    if(nrow(sigtab) > 0) {
      sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy)[rownames(sigtab), ], "matrix"))
      sigtab <- tibble::rownames_to_column(as.data.frame(sigtab), var="OTUname")
      if(is.null(file)) {
        return(sigtab)
      } else {
        write.csv(x=sigtab, file=file)
      }
    } else {
      return(paste("No significant results at p(adj) <", alpha))
    }
  } else {
    if(is.null(file)) {
      return(print(res))
    } else {
      write.csv(x=res, file=file)
    }
  }
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

make_phyloseq <- function(otu_table, map_file, sampleid="#SampleID") {
  #This is a function read in above to create the taxonomy data needed for phylsoeq
  split <- split.tax(cbind(OTUID=rownames(otu_table), taxonomy=as.character(otu_table[,"taxonomy"])))
  # Now remove the taxonomy variable from your OTU table because it confuses phyloseq:
  otu_table[,"taxonomy"] <- NULL
  print(paste("Assuming", sampleid, "is your sample identifier. To change this, search for the 'make_phyloseq' function and use the 'sampleid' option."))
  if(sum(is.na(map_file[,sampleid]))>0) {
  new_map_file <- map_file[!is.na(map_file[,sampleid]),]
  print("Missing IDs found. Mapping file will be subset to non-missing Sample IDs.")
  } else {
    new_map_file <- map_file
  }
      rownames(new_map_file) <- new_map_file[,sampleid]
      phy.obj <- phyloseq(otu_table(otu_table, taxa_are_rows = TRUE), tax_table(split), sample_data(new_map_file))
      print(phy.obj)
      return(phy.obj)
}

######

#Now that the data and functions have been read in, analysis begins below:

#####
phy <- make_phyloseq(my_otu_table, my_map)

phy.prune <- prune_taxa(taxa_sums(phy) > 0, phy)
#This removes any missing values from your variable/data frame (causes an error when running phyloseq_to_deseq2):
phy.prune <- prune_samples(as.vector(!is.na(sample_data(phy.prune)[,var_of_int])), phy.prune)
phy.prune
my_formula <- as.formula(paste("~", var_of_int))
deseqdat <- phyloseq_to_deseq2(phy.prune, my_formula)
deseqdat2 = estimateSizeFactors(deseqdat, type="poscounts")
deseqdat2 <- DESeq(deseqdat2, fitType = "local")
print_res(deseqdat2, phy.prune, ref=ref_group, cont=cont_group, file=output, alpha=alpha)

