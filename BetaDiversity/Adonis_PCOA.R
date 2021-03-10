adonis.pcoa <- function(phyloseq, distmethod, variable) {
  distmat <- phyloseq::distance(phyloseq, method=distmethod)
  sampledf <- data.frame(sample_data(phyloseq))
  sampledf <- sampledf[labels(distmat),]
  r2stat <- adonis2(as.formula(paste0("distmat ~", variable)), data=sampledf)
  tolab <- paste("R^2=", round(r2stat$R2[1],3), "\n","P=", round(r2stat$`Pr(>F)`[1],3))
  mb.ord <- ordinate(phyloseq, method="PCoA", dist=distmat)
  p1.dat <- plot_ordination(phyloseq, mb.ord, type="sample", justDF=TRUE)
  p1 <- ggplot(p1.dat, aes(x=Axis.1, y=Axis.2, fill=p1.dat[,variable])) + 
    geom_point(shape=21, size=3, color="black") +
    stat_ellipse(level=0.68) +
    geom_text(aes(x=max(Axis.1), y=min(Axis.2), label=tolab), hjust=1, vjust=-0.5, inherit.aes = FALSE, check_overlap = TRUE) +
    xlab(paste("Axis 1 [", round(mb.ord$values$Eigenvalues[1],2),"%]")) +
    ylab(paste("Axis 2 [", round(mb.ord$values$Eigenvalues[2],2),"%]")) +
    scale_fill_manual(values=ucsfColSch) +
    ggtitle(distmethod)
  return(p1)
}
