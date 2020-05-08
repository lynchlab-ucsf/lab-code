# Code to remove negative controls.
# Katie McCauley
# Is it possible to develop this code so that it plays well with the pipeline so far.
#Maybe use the input CSV to determine which samples belong to each run and do run-specific negative control filtering?

setwd("/data/Users/kmccauley/CREW/DADA2_Comb/")
otu <- read.table("qiime_otu_table_filt.txt", header=TRUE, check.names=F, sep="\t", comment="", skip=1, row.names=1)

neg.dat <- otu[, grepl("NTC", names(otu)) & !names(otu) %in% c("taxonomy")]
split_tax <- strsplit(as.character(otu$taxonomy), "; ")
genus_level <- lapply(split_tax, function(x) {
  x <- x[-7]
  x[length(x)]
})
genus_name <- sapply(genus_level, paste)
otu2 <- cbind(otu, genus_name)
otu2[1:10,1:10]

negsA <- rowSums(neg.dat>0)/ncol(neg.dat)
#otu3 <- otu2[otu2$genus_name %in% common_contaminant,]
samp.dat <- otu2[,!grepl("NTC", names(otu2)) & !names(otu2) %in% c("#OTU ID","taxonomy","genus_name")]
#The proportion of samples with counts for each OTU
sampsA <- rowSums(samp.dat>0)/ncol(samp.dat)
dat <- merge(negsA, sampsA, by=0)
names(dat) <- c("#OTU ID","negs","samps")
dat2 <- merge(dat, otu2[,c("taxonomy","genus_name")], by.x="#OTU ID",by.y=0)
levels(dat2$genus_name)[levels(dat2$genus_name) %in% "NA"] <- "Other" #Most-frequent OTU classifications

library(ggplot2)
pre.cleaning <- ggplot(dat2, aes(negs,samps)) + 
  geom_point() + 
  ylab("Proportion of Samples") + 
  xlab("Proportion of Negative Controls") + 
  ggtitle("") + 
  ggrepel::geom_label_repel(data=subset(dat2, negs > 0.1), aes(label=genus_name)) +
  labs(tag="A")
#ggsave("Pre_Cleaning_Figure.pdf", pre.cleaning, device="pdf", height=6, width=7)

outright.drop <- dat2$`#OTU ID`[dat2$negs > 0.15 & dat2$samps < 0.15]

#ord = outright.drop
otu.ord <- otu2[!rownames(otu2) %in% outright.drop,]

#otu.sub
neg.dat2 <- neg.dat[rowSums(neg.dat) >0,]
mean.in.NTC <- ceiling(apply(neg.dat2[!rownames(neg.dat2) %in% outright.drop,],1, mean)) # Want a round number, so I am taking the ceiling.
max.in.NTC <- apply(neg.dat2[!rownames(neg.dat2) %in% outright.drop,],1, max)

otu.notax <- otu.ord[, !names(otu.ord) %in% c("taxonomy","genus_name")]


otu.notax.mean <- t(otu.notax)
for(i in names(mean.in.NTC)) {
  otu.notax.mean[, i] <- otu.notax.mean[, i] - mean.in.NTC[i]
  otu.notax.mean[, i][otu.notax.mean[, i] < 0] <- 0
}
otu.notax.mean <- data.frame(t(otu.notax.mean), check.names=F)

otu.notax.max <- t(otu.notax)
for(i in names(max.in.NTC)) {
  otu.notax.max[, i] <- otu.notax.max[, i] - max.in.NTC[i]
  otu.notax.max[, i][otu.notax.max[,i] < 0] <- 0
}
otu.notax.max <- data.frame(t(otu.notax.max), check.names=F)

neg.dat <- otu.notax.mean[, grepl("NTC", names(otu.notax.mean))]
samp.dat <- otu.notax.mean[, !grepl("NTC", names(otu.notax.mean))]
negsA <- rowSums(neg.dat>0)/ncol(neg.dat)
sampsA <- rowSums(samp.dat>0)/ncol(samp.dat)
dat <- merge(negsA, sampsA, by=0)
names(dat) <- c("#OTU ID","negs","samps")
dat2 <- merge(dat, otu2[,c("taxonomy","genus_name")], by.x="#OTU ID",by.y=0)
library(ggplot2)
post.clean <- ggplot(dat2, aes(negs,samps)) + 
  geom_point() + 
  ylab("Proportion of Samples") + 
  xlab("Proportion of Negative Controls") + 
  ggtitle("") + 
  ggrepel::geom_label_repel(data=subset(dat2, negs > 0.05), aes(label=paste(genus_name, `#OTU ID`))) +
  scale_color_brewer(palette="Paired") +
  labs(tag="B")
library(gridExtra)
ggsave("Combined_Cleaning_Figure.pdf", grid.arrange(pre.cleaning, post.clean, nrow=1), device="pdf", height=5, width=14)

if(identical(rownames(samp.dat), rownames(otu.ord))) samp.dat.tax <- cbind(samp.dat, taxonomy=otu.ord$taxonomy)
write.table(samp.dat.tax, "otu_table_noneg.txt", col.names=TRUE, row.names=TRUE, sep="\t")
