# DM is a distance matrix, Factor is the variable you want to split by, dat is your mapping file, and name is the name of the figure file to be generated
# Uses user-defined cbPalette, which can be found here: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# You may also need to change all instances of "#SampleID" to your sample variable.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pairwise.figs <- function(dm, factor, dat, name, sampid="#SampleID") {
  pacman::p_load(vegan, ape, ggplot2, gridExtra)
  pair <- NULL
  pair <- t(combn(levels(dat[,factor]),2))
  tempPal <- cbPalette[c(1:(length(levels(dat[,factor]))+1))]
  names(tempPal) <- levels(dat[,factor])
  set.seed(123)
  over.ad <- adonis2(as.formula(paste0("dm ~", factor)), data=dat)
  over.R2 <- round(over.ad$R2[1],3)
  over.P <- over.ad$`Pr(>F)`[1]
  mypcs <- pcoa(dm)
  pcs <- data.frame(mypcs$vectors)[,1:3]
  pcs.meta <- merge(pcs, dat, by.x=0, by.y=sampid)
  myplot <- list()
  pcs.meta1 <- data.frame("PC1"=pcs.meta$Axis.1, "PC2"=pcs.meta$Axis.2, "var" = pcs.meta[,factor])
  myplot[[1]] <- ggplot(pcs.meta1, aes(x=PC1, y=PC2, color=var)) +
    geom_point() +
    ggtitle("Overall") +
    annotate("text", x=-Inf,y=Inf, vjust=1.2,hjust=0, label=paste("italic(R)^2 == ", over.R2),parse=TRUE) +
    annotate("text", x=-Inf,y=Inf, vjust=3.2,hjust=0, label=paste("italic(P) == ", over.P), parse=TRUE) +
    scale_color_manual(" ", values=tempPal) +
    stat_ellipse(level=0.68)
  for(i in 2:(nrow(pair)+1)) {
    sub.dat <- dat[dat[,factor] %in% c(pair[(i-1),1],pair[(i-1),2]),]
    if(sampid != 0) {
    sub.dm <- dm[as.character(sub.dat[,sampid]),as.character(sub.dat[,sampid])]
    } else {
      sub.dm <- dm[rownames(sub.dat), rownames(sub.dat)]
    }
    set.seed(123)
    ad <- adonis2(as.formula(paste0("sub.dm ~", factor)), data=sub.dat)
    R2 <- round(ad$R2[1],3)
    P <- ad$`Pr(>F)`[1]
    mypcs <- pcoa(sub.dm)
    pcs <- data.frame(mypcs$vectors)[,1:3]
    pcs.meta <- merge(pcs, sub.dat, by.x=0,by.y=sampid)
    pcs.meta1 <- data.frame("PC1" = pcs.meta$Axis.1, "PC2" = pcs.meta$Axis.2, "var" = pcs.meta[,factor])
    p <- ggplot(data=pcs.meta1, aes(x=PC1,y=PC2,color=var)) +
      geom_point() +
      scale_color_manual(values=tempPal) + 
      guides(color=FALSE) + 
      ggtitle(paste(pair[i-1,1], "vs", pair[i-1,2])) +
      theme(title = element_text(size=7)) +
      stat_ellipse(level=0.68) +
      annotate("text", x=-Inf,y=Inf, vjust=1.2,hjust=0, label=paste("italic(R)^2 == ", R2),parse=TRUE) +
      annotate("text", x=-Inf,y=Inf, vjust=3.2,hjust=0, label=paste("italic(P) == ", P), parse=TRUE)
    myplot[[i]] <- p
  }
myplot[[i + 1]]<- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + theme_void()
if(i < 6) {
  lay.mat <- rbind(rep(1,(i-1)), rep(1,(i-1)), c(2:i))
} else {
  # This will need to be modified for other's plots.
  lay.mat <- rbind(rep(1,(4)),rep(1,(4)), rep(1,(4)), c(2:5), c(6:9), c(12,10:11,12))
}
ggsave(paste0(name,".pdf"), grid.arrange(grobs=myplot, layout_matrix=lay.mat), device="pdf", width=8, height=11)
}
