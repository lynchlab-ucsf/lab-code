# Cluster Functions
## Katie McCauley
## Last Updated: Feb 05, 2019

## These are just functions. DO NOT COPY/MODIFY THIS FILE. SEE USAGE EXAMPLE BELOW!
## Developed to identify clusters in data. Only two functions are needed. The rest are "helper functions".

## Function details:
# hclust_func <- function(distmat, phy, clustno=which.max(avg.sil), result_prefix = NULL, run.RR=FALSE)
# pam_func <- function(distmat, phy, clustno, result_prefix=NULL, run.RR=FALSE)
# distmat = distance matrix. I suggest generating from phyloseq. Code has only been tested with the use of a phyloseq object. I am concerned that the use of an outside distance matrix will not match with the phyloseq object
# phy = your phyloseq object.
# clustno= the number of clusters you want to enforce on your data. For hclust, if you don't specify a cluster number, it will choose the number with the largest silhouette statistic. For PAM, you will need to set a cluster number initially, review the gap statistic, and then choose an "official" cluster number. A future iteration may provide only the gap statistics when a clustno is not provided. However, this is not implemented here
# result_prefix = if you would like to print the resulting figures, choose a prefix that identifies your output. If NULL, I generate an automatic name given the "distmat" name. If NA, plots will be printed to the console and not saved.
# run.RR = Whether or not you want to run Relative Risk analysis. This was study-specific. KEEP AS FALSE UNLESS YOU REVIEW THE "RR.func" FUNCTION! I will not support updates or modifications to the function at this moment.

## Usage in your own code:
# source("/data/Users/kmccauley/LabCode/Clustering_Functions.R")
# {Read in your phyloseq object}
# {Generate Distance Matrices from Phyloseq Object}
# hclust_func(mydist, myphy, 3, result_prefix=NA)

pacman::p_load(ape, ggplot2, cluster, phyloseq, data.table, vegan)

## From: https://github.com/joey711/phyloseq/issues/418
plot_clusgap = function(clusgap, title="Gap Statistic calculation results"){
require("ggplot2")
gstab = data.frame(clusgap$Tab, k=1:nrow(clusgap$Tab))
p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size=5)
p = p + geom_errorbar(aes(ymax=gap+SE.sim, ymin=gap-SE.sim))
p = p + ggtitle(title)
print(p)
}
pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

summarize_taxa = function(physeq, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq),
                     var1 = get_variable(physeq, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt = mdt[, list(meanRA = mean(RelativeAbundance),
                         sdRA = sd(RelativeAbundance),
                         minRA = min(RelativeAbundance),
                         maxRA = max(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}

RR.func <- function(func_dat=temp.V12, var="cluster", out="Case.or.Control.Status.Original") {
  data <- data.frame(sample_data(func_dat))
  for(i in levels(factor(data[[var]]))) {
    print(i)
    print(summary(glm(data$cluster %in% i ~ factor(data[[out]]), data=data, family=binomial))$coef[,c(1,4)])
  }
}

hist.plot <- function(physeq, level="Genus", myPalette=cbPalette) {
  require(stringr)
  if(level != "Genus.Species") {
    summs <- summarize_taxa(physeq, level)
  } else {
    tax_table(physeq)[tax_table(physeq) %in% "NA"] <- ""
    Genus.Species <- paste(tax_table(physeq)[,"Genus"], tax_table(physeq)[,"Species"], sep=" ")
    tax_table(physeq) <- cbind(tax_table(physeq), Genus.Species)
    summs <- summarize_taxa(physeq, "Genus.Species")
  }
  phy.merge.samp <- merge_samples(physeq, "cluster")
  genus.list <- c(summs[order(summs$maxRA, decreasing = TRUE),][1:length(myPalette)][[level]])
  tax_table(phy.merge.samp)[!tax_table(phy.merge.samp)[,level] %in% genus.list, level] <- NA
  phy.merge.samp <- transform_sample_counts(phy.merge.samp, function(x) x / sum(x))
  plt <- plot_bar(phy.merge.samp, fill=level) + 
    scale_fill_manual(values=myPalette) + 
    xlab("Cluster") +
    theme(legend.title=element_blank(), legend.text=element_text(size=7))
  return(plt)
}


hclust_func <- function(distmat, phy, clustno=which.max(avg.sil), result_prefix = NULL, run.RR=FALSE, 
                        level="Genus.Species", colorPalette=cbPalette, return_dataframe=FALSE) {
  if(is.null(result_prefix)) result_prefix <- paste0(deparse(substitute(distmat)),"_hclust")
  clust.dat <- distmat
  phy <- phy
  myhclust <- hclust(as.dist(clust.dat), method="ward.D2")
  cluster.opts <- 2:18
  avg.sil <- NULL
    for(i in cluster.opts) {
      cuts <- cutree(myhclust, k=i)
      avg.sil[i] <- summary(silhouette(cuts, as.dist(clust.dat)))$si.summary["Mean"]
    }
  if(!is.na(result_prefix)) {
    ggsave(paste0(result_prefix, "_silhouetteplot.pdf"),  plot(avg.sil), device="pdf", height=5, width=7)
  } else {
    plot(avg.sil)
  }

  choose.clust <- cutree(myhclust, k=clustno)
  pcs <- pcoa(clust.dat)$vectors[,1:2]
  pcs.dat <- merge(pcs, choose.clust, by=0)
  rownames(pcs.dat) <- pcs.dat$Row.names
  pcs.dat$Row.names <- NULL

  temp <- merge(sample_data(phy), pcs.dat, by=0)
  rownames(temp) <- temp$Row.names
  temp$Row.names <- NULL
  temp$cluster <- factor(temp$y, labels=paste0("Cluster",1:length(unique(temp$y))))
  plt <- ggplot(temp, aes(x=Axis.1,y=Axis.2, color=cluster)) + 
    geom_point() + 
    stat_ellipse(level=0.65)
  if(!is.na(result_prefix)) {
    ggsave(paste0(result_prefix, "_PCoAplot.pdf"),  plt, device="pdf", height=5, width=7)
  } else {
    print(plt)
  }
  print(table(temp$cluster))

  temp <- temp[labels(clust.dat),]
  print(adonis2(clust.dat ~ cluster, data=temp))
if(run.RR == TRUE) {
  ## Risk Ratios for each cluster:
  temp.V12 <- temp[!temp$Analysis.Visit %in% "Visit 0",]
  temp.V0 <- temp[temp$Analysis.Visit %in% "Visit 0" & !temp$Case.or.Control.Status.Full.Cohort %in% "",]

  print("Exacerbation Risk by Cluster")
  RR.func(temp.V12, var="cluster","Case.or.Control.Status.Original")
  print("Exacerbation Risk based on Baseline by Cluster")
  RR.func(temp.V0, var="cluster", "Case.or.Control.Status.Full.Cohort")
  
  print("Viral Detection Risk by Cluster")
  RR.func(temp.V12, var="cluster","Virus_pos")
  print("Viral Detection at Baseline by Cluster")
  RR.func(temp.V0, var="cluster", "Virus_pos")
}
  phy2 <- merge_phyloseq(phy, sample_data(temp))
  if(!is.na(result_prefix)) {
    ggsave(paste0(result_prefix, "_clusterplot.pdf"),  hist.plot(phy2, level=level, myPalette=colorPalette), device="pdf", width=10, height=5)
  } else {
    hist.plot(phy2, level=level, myPalette=colorPalette)
  }
  if(return_dataframe != FALSE) {
    return(sample_data(phy2))
  }
}

pam_func <- function(distmat, phy, clustno, result_prefix=NULL, run.RR=FALSE, level="Genus.Species",
                     colorPalette=cbPalette, return_dataframe=FALSE) {
  if(is.null(result_prefix)) result_prefix <- paste0(deparse(substitute(distmat)),"_pam")
  phy <- phy
  clust.dat <- distmat
  pcs <- pcoa(clust.dat)$vectors[,1:15]
  find.clust <- clusGap(pcs, FUN=pam1, K.max=18, B=50)
  
  if(!is.na(result_prefix)) {
    ggsave(paste0(result_prefix, "_gapplot.pdf"),  plot_clusgap(find.clust), device="pdf", height=5, width=7)
  } else {
    plot_clusgap(find.clust)
  }

  clusters <- pam1(pcs, clustno)
  pcs.dat <- merge(pcs, clusters, by.x=0, by.y=0)
  pcs.dat$cluster <- factor(paste0("Cluster", pcs.dat$cluster))
  plt <- ggplot(pcs.dat, aes(x=Axis.1, y=Axis.2, color=cluster)) + 
    geom_point() + 
    stat_ellipse(level=0.68)
  if(!is.na(result_prefix)) {
    ggsave(paste0(result_prefix, "_PCoAplot.pdf"),  plt, device="pdf", height=5, width=7)
  } else {
    print(plt)
  }

  row.names(pcs.dat) <- pcs.dat$Row.names
  pcs.dat$Row.names <- NULL
  print(table(pcs.dat$cluster))

  pcs.dat <- pcs.dat[labels(clust.dat),]
  print(adonis2(clust.dat ~ cluster, data=pcs.dat))

  phy.dat <- merge_phyloseq(phy, sample_data(data.frame(pcs.dat)))
if(run.RR==TRUE) {
  temp.V12 <- subset_samples(phy.dat, !Analysis.Visit %in% "Visit 0")
  temp.V0 <- subset_samples(phy.dat, Analysis.Visit %in% "Visit 0" & !Case.or.Control.Status.Full.Cohort %in% "")

  print("Exacerbation Risk by Cluster")
  RR.func(temp.V12, var="cluster","Case.or.Control.Status.Original")
  print("Exacerbation Risk based on Baseline by Cluster")
  RR.func(temp.V0, var="cluster", "Case.or.Control.Status.Full.Cohort")
  
  print("Viral Detection Risk by Cluster")
  RR.func(temp.V12, var="cluster","Virus_pos")
  print("Viral Detection at Baseline by Cluster")
  RR.func(temp.V0, var="cluster", "Virus_pos")
}
  if(!is.na(result_prefix)) {
    ggsave(paste0(result_prefix, "_clusterplot.pdf"),  hist.plot(phy.dat, level=level, myPalette=colorPalette), device="pdf", width=10, height=5)
  } else {
    hist.plot(phy.dat, level=level, myPalette=colorPalette)
  }
  if(return_dataframe != FALSE) {
    return(sample_data(phy.dat))
  }
}

