#Functions Needed for DESeq Analysis:

## Last Update: Nov 22, 2019
#### I determined that I needed to add a require(ggplot2) for the plot_res function.

split.tax <- function(tax.dat) {
  taxanames <- strsplit(as.character(tax.dat[,2]),"; ")
  mat <- t(sapply(taxanames,
                  function(x,m) c(x,rep(NA,m-length(x))),
                  max(rapply(taxanames,length))))
  
  newnames <- gsub("_","",mat)
  newnames <- as.matrix(newnames)
  colnames(newnames) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  row.names(newnames) <- tax.dat[,1]
  newnames[,6][newnames[,6] %in% c("") | is.na(newnames[,6])] <- newnames[,5][newnames[,6] %in% c("") | is.na(newnames[,6])]
  newnames[,6][newnames[,6] %in% c("") | is.na(newnames[,6])] <- newnames[,4][newnames[,6] %in% c("") | is.na(newnames[,6])]
  return(newnames)
}

split.tax_SILVA <- function(tax.dat) {
  taxanames <- strsplit(as.character(tax.dat[,2]),"; ")
  mat <- t(sapply(taxanames,
                  function(x,m) c(x,rep(NA,m-length(x))),
                  max(rapply(taxanames,length))))
  
  newnames <- gsub("_","",mat)
  newnames <- as.matrix(newnames)
  colnames(newnames) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  row.names(newnames) <- tax.dat[,1]
  newnames[,6][grepl("uncultured", newnames[,6][newnames[,6]] ) | is.na(newnames[,6])] <- newnames[,5][grepl("uncultured", newnames[,6][newnames[,6]]) | is.na(newnames[,6])]
  newnames[,6][grepl("uncultured", newnames[,6][newnames[,6]] ) | is.na(newnames[,6])] <- newnames[,4][grepl("uncultured", newnames[,6][newnames[,6]]) | is.na(newnames[,6])]
  return(newnames)
}

print_res <- function(x, phy, alpha=0.1, sig.only=TRUE, ref= NULL,  cont=NULL, var=var_of_int, file=NULL) {
  if(is.null(ref) & is.null(cont)) {
    # This is where the code would be modified to print all possible comparisons...
    res <- results(x, pAdjustMethod = "fdr")
    if(class(sample_data(phy)[[var]]) %in% "factor") {
    varlevs <- levels(sample_data(phy)[[var]])
    res$comparison <- paste0(varlevs[1], "_vs_", varlevs[length(varlevs)])
    }
  } else {
    res <- results(x, contrast=c(var, cont, ref), pAdjustMethod = "fdr")
    #res$baseMeanA <- 2^(mcols(x)[, startsWith(names(mcols(x)), var)])
    #res$baseMeanB <- 2^(mcols(x)[, "Intercept"] + mcols(x)[, startsWith(names(mcols(x)), var)])
    res$comparison <- paste0(ref, "_vs_", cont)
    res$ref <- ref
    res$cont <- cont
  }
  res = res[!is.na(res$padj), ]
  res = res[order(res$log2FoldChange), ]
  res = cbind(as(res, "data.frame"), as(tax_table(phy)[rownames(res), ], "matrix"))
  res <- tibble::rownames_to_column(as.data.frame(res), var="OTUname")
  if(sig.only==TRUE) {
    sigtab = res[(res$padj < alpha), ]
    if(nrow(sigtab) > 0) {
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
      return(res)
    } else {
      write.csv(x=res, file=file)
    }
  }
}

plot_res <- function(x, phy=NULL, alpha=0.05, ref= NULL,  cont=NULL, var=var_of_int, diff_cut=0, file=NULL, color=TRUE) {
  require(ggplot2)
  if(!is.null(phy)) { #If you have a phyloseq object, do this:
    if(is.null(ref) & is.null(cont)) {
      res <- results(x, pAdjustMethod = "fdr")
      if(class(sample_data(phy)[[var]]) %in% "factor") {
      varlevs <- levels(sample_data(phy)[[var]])
      res$comparison <- paste0(varlevs[1], "_vs_", varlevs[length(varlevs)])
      diffs <- sapply( levels(x[,var]), function(lvl) rowMeans( counts(x,normalized=TRUE)[,x[,var] == lvl, drop=F] ) )
      res <- cbind(res, diffs)
      }
    } else {
      res <- results(x, contrast=c(var, cont, ref), pAdjustMethod = "fdr", minmu=1)
      res$comparison <- paste0(ref, "_vs_", cont)
      diffs <- sapply( levels(x[[var]]), function(lvl) rowMeans( counts(x,normalized=TRUE)[,x[[var]] == lvl, drop=F] ) )
      res <- cbind(res, diffs)
    }
  } else { # X object is dataframe with significant taxa (also assuming the bacteria names are already included)
    res <- x
  }
  res = res[order(res$padj, na.last=NA), ]
  #x$myvar <- c(sample_data(phy)[,var])[[1]]
  #baseMeanPerLevel <- sapply(levels(x$myvar), function(lvl) rowMeans( counts(x, normalized=TRUE)[,x$myvar %in% lvl]))
  #delta = baseMeanPerLevel[,cont] - baseMeanPerLevel[,ref]
  sigtab = res[!is.na(res$padj), ]
  #sigtab <- merge(sigtab, delta, by=0)
  #names(sigtab)[names(sigtab) %in% "y"] <- "delta"
  #rownames(sigtab) <- sigtab$Row.names
  #sigtab$Row.names <- NULL
  sigtab$delta <- ifelse(sigtab$log2FoldChange > 0, 2^sigtab$log2FoldChange * sigtab$baseMean, -(2^sigtab$log2FoldChange * sigtab$baseMean))
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy)[rownames(sigtab), ], "matrix"))
  sigtab <- tibble::rownames_to_column(as.data.frame(sigtab), var="OTUname")
  sigtab$id.group_tag <- sigtab$padj < alpha & abs(sigtab$delta) > diff_cut & sigtab$OTUname %in% sigtab[order(sigtab$padj),][1:20,]$OTUname
  sigtab$id.group_sig <- sigtab$padj < alpha
  sigtab$Genus.color <- NA
  sigtab$Genus.color[sigtab$id.group_sig] <- as.character(sigtab$Genus[sigtab$id.group_sig])
  ## Make the dashed line look less like points
  ## Color code by the dominant genus (same color code as used in the paper)
  if(color==TRUE) {
  fig <- ggplot(sigtab, aes(x=log2FoldChange, y=-log(padj), color=Genus.color, size=baseMean)) + 
    geom_point()  +
    geom_hline(yintercept=-log(alpha), color="grey", linetype="longdash", size=1) +
    annotate("label", x=max(sigtab$log2FoldChange)-(abs(max(sigtab$log2FoldChange))*.2), y=0, label=cont) +
    annotate("label", x=min(sigtab$log2FoldChange)+(abs(min(sigtab$log2FoldChange))*.2), y=0, label=ref) +
    ggrepel::geom_text_repel(data=sigtab[sigtab$id.group_tag, ], aes(label=OTUname), color="black",
                              min.segment.length=0.01,size=3, show.legend=FALSE, point.padding = 1, max.iter = 1e5, parse=TRUE, force=10) +
    guides(size=FALSE) +
    ylab("Q-value (-log)") +
    xlab("Log2 Fold Change")
  fig
  } else {
    fig <- ggplot(sigtab, aes(x=log2FoldChange, y=-log(padj))) + 
      geom_point(aes(color=id.group_sig))  +
      scale_color_manual("", values=c("grey50","red"), labels=c("Not Significant","FDR < 0.05")) +
      geom_hline(yintercept=-log(alpha), color="grey", linetype="longdash", size=1) +
      annotate("label", x=max(sigtab$log2FoldChange)-(abs(max(sigtab$log2FoldChange))*.2), y=0, label=cont) +
      annotate("label", x=min(sigtab$log2FoldChange)+(abs(min(sigtab$log2FoldChange))*.2), y=0, label=ref) +
      ggrepel::geom_text_repel(data=sigtab[sigtab$id.group_tag, ], aes(label=paste(Genus, OTUname)), color="black",
                               min.segment.length=0.01,size=3, point.padding = 0.2, max.iter = 1e5, force=10) +
      ylab("Q-value (-log)") +
      xlab("Log2 Fold Change")
    fig
  }
  if(!is.null(file)) {
    ggsave(file, plot=fig, device="pdf", height=6, width=6)
  }
  if(is.null(file)) {
    return(fig)
  }
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

make_phyloseq <- function(otu_table, map_file, sampleid="#SampleID", is.SILVA=FALSE) {
  require(phyloseq)
  #This is a function read in above to create the taxonomy data needed for phylsoeq
  if(is.SILVA == FALSE) {
  split <- split.tax(cbind(OTUID=rownames(otu_table), taxonomy=as.character(otu_table[,"taxonomy"])))
  } else {
    split <- split.tax_SILVA(cbind(OTUID=rownames(otu_table), taxonomy=as.character(otu_table[,"taxonomy"])))
  }
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
