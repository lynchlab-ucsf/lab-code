plot_res <- function(x, phy, alpha=0.05, ref= NULL,  cont=NULL, var=var_of_int, diff_cut=0, file=NULL) {
  if(!is.null(phy)) { #If you have a phyloseq object, do this:
  if(is.null(ref) & is.null(cont)) {
    res <- results(x, pAdjustMethod = "fdr")
    varlevs <- levels(sample_data(phy)[[var_of_int]])
    res$comparison <- paste0(varlevs[1], "_vs_", varlevs[length(varlevs)])
  } else {
    res <- results(x, contrast=c(var, cont, ref), pAdjustMethod = "fdr")
    res$comparison <- paste0(ref, "_vs_", cont)
  }
  } else { # X object is dataframe with significant taxa (also assuming the bacteria names are already included)
    res <- x
  }
  res = res[order(res$padj, na.last=NA), ]
  sigtab = res[!is.na(res$padj), ]
  sigtab$delta <- sigtab$log2FoldChange*sigtab$baseMean
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy)[rownames(sigtab), ], "matrix"))
  sigtab <- tibble::rownames_to_column(as.data.frame(sigtab), var="OTUname")
  sigtab$id.group <- sigtab$padj < alpha & abs(sigtab$delta) > diff_cut
  fig <- ggplot(sigtab, aes(x=delta, y=-log(padj), color=id.group)) + 
    geom_point() +
    geom_label(aes(x=max(delta)-(abs(max(delta))*.2), y=0, label=cont), color="black", show.legend=FALSE) +
    geom_label(aes(x=min(delta)+(abs(min(delta))*.2), y=0, label=ref), color="black", show.legend=FALSE) +
    ggrepel::geom_text_repel(data=sigtab[sigtab$id.group, ], aes(label=paste(Genus, OTUname)), size=2, show.legend=FALSE) +
    scale_color_discrete("", labels=c("Adjusted P < 0.01", "Adjusted P > 0.01")) +
    ylab("Adjusted P-value (-log)") +
    xlab("Estimated Delta")
  fig
  if(!is.null(file)) {
    ggsave(file, plot=fig, device="pdf", height=6, width=6)
  }
  if(is.null(file)) {
    return(fig)
  }
}

