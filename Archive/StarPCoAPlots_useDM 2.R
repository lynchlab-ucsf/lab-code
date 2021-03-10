# Function to Generate Star Plots (PCOA plots with line segments to the centroid of each group)
## Written by Katie McCauley on Sept 14, 2017

#First argument is your distance matrix, second argument is the variable you want to generate the distance to, third argument is the data frame that contains the samples and variable for mapping.
## After running the function, type the following to use it.
### star.figs(my.dist.matrix, "Variable.For.Grouping", my.dataset)

#If you have your own color scheme that you would prefer to use, you can use the term mypalette to pass a string of colors for use.

star.figs <- function(dm, variable, dat, mypalette=NULL, sampid=0) {
  set.seed(123)
  if(!is.null(mypalette)) {
     tempPal <- mypalette
     names(tempPal) <- levels(factor(dat[,variable]))
  } else {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    tempPal <- gg_color_hue(length(levels(dat[,variable])))
    names(tempPal) <- levels(dat[,variable])
  }
  if(sampid == 0) {
    sub.dm <- dm[rownames(dat), rownames(dat)]
  } else {
    sub.dm <- dm[as.character(dat[,sampid]), as.character(dat[,sampid])]
  }
  
  over.ad <- adonis(as.dist(sub.dm) ~ dat[,variable])
  over.R2 <- round(over.ad$aov.tab[1,5],3)
  over.P <- round(over.ad$aov.tab[1,6],3)
  mypcs <- pcoa(sub.dm)
  pcs <- data.frame(mypcs$vectors)[,1:2]
  pcs.meta <- merge(pcs, dat, by.x=0, by.y=sampid)
  centroids <- aggregate(cbind(Axis.1,Axis.2)~pcs.meta[,variable], data=pcs.meta, mean)
  gg <- merge(pcs.meta, centroids, by.x=variable, by.y=1, suffixes=c("",".centroid"))
  myplot <- list()
  pcs.meta1 <- data.frame("PC1"=gg$Axis.1, "PC2"=gg$Axis.2, "var" = pcs.meta[,variable])
  myplot <- ggplot(gg, aes(x=Axis.1, y=Axis.2, color=gg[,variable])) +
    geom_point() +
    ggtitle("Overall") +
    scale_color_manual(" ", values=tempPal) +
    annotate("text", x=-Inf,y=Inf, vjust=1.2,hjust=0, label=paste("italic(R)^2 == ", over.R2),parse=TRUE) +
    annotate("text", x=-Inf,y=Inf, vjust=3.2,hjust=0, label=paste("italic(P) == ", over.P), parse=TRUE) +
    geom_segment(aes(x=Axis.1.centroid,y=Axis.2.centroid, xend=Axis.1, yend=Axis.2)) +
    geom_point(aes(x=Axis.1.centroid, y=Axis.2.centroid), size=5, shape=21,fill="white", stroke=5, show.legend=FALSE) +
    stat_ellipse(type="norm", level=0.68)
  print(myplot)
}

#star.figs(bray.dm, "pheno_y7", dat=mapping.pheno.1yr, sampid="#SampleID")
