# DADA2 Post-Processing.
#Takes the files created in DADA2 and makes them useful/qiime compatible. The names of these files are consistent with the names given in the output files from my "DADA2_Script.R" code.
# I also rename sequences as "OTUs" here, and save the sequences elsewhere. 
# I suggest running this with "nohup", since the "optimal tree" generation can take up to 12 hours. However, the OTU table and a temporary tree should print out relatively quickly.

# Author: Katie McCauley

rm(list=ls())
setwd(".")
dada2_otu <- read.table("dada2_otutable.txt", header=TRUE, check.names=F, comment="", sep="\t", row.names=1)
dada2_tax <- read.table("dada2_tax_table.txt", header=TRUE, check.names=F, comment="", sep="\t")
dada2_seqs <- read.table("dada2_seqs.txt", header=F, comment="", sep="\t")

#Make the sample names actual sample names instead of file names:
listnames <- strsplit(names(dada2_otu), split="_", fixed=TRUE)
names(dada2_otu) <- sapply(lapply(listnames, function(x) x[-c(length(x),length(x)-1) ]), paste, collapse="_")
names(dada2_otu) <- sub("run1_|run2_|run3_", "", names(dada2_otu))

#Put taxonomy at the end of the data (in QIIME format)
dada2_tax$taxonomy <- apply(dada2_tax, 1, paste0, collapse="; ")
dada2_otu$taxonomy <- dada2_tax$taxonomy

pacman::p_load(dplyr, tibble)
dada2_otuA <- dada2_otu %>% tibble::rownames_to_column("#OTUID")

# Rename the sequences as OTUs:
seqs <- as.character(dada2_seqs[,1])
if(sum(rownames(dada2_otu) != seqs) == 0) {
rownames(dada2_otu) <- paste0("SV_", 1:nrow(dada2_otu))
} else {
print("Sequence file does not match SV row names.")
stop()
}

pacman::p_load(data.table)
tab <- data.table(V1=paste0(">",rownames(dada2_otu)), V2=dada2_seqs, order=1:nrow(dada2_seqs))
header <- data.frame(V1=tab$V1, V2=tab$order)
reads <- data.frame(V1=tab$V2, V2=tab$order)
comb <- rbind(header, reads)
comb <- comb[order(comb$V2),]
write.table(comb$V1, "dada2_fasta.fa", quote=F, row.names=F, sep="\t", col.names=F)
tab2 <- data.frame(V1=rownames(dada2_otu), V2=dada2_seqs)
write.table(tab2, "dada2_seq_OTU_tab.txt", quote=F, row.names=F, sep="\t", col.names=F)

pacman::p_load(dplyr, tibble)
dada2_otu <- dada2_otu %>% tibble::rownames_to_column("#OTUID")

write.table(dada2_otu, "qiime_otu_table.txt", col.names=TRUE, quote=F, sep="\t", row.names=F)
write.table(dada2_otuA, "seqs_otu_table.txt", col.names=TRUE, quote=F, sep="\t", row.names=F)
print("OTU Table Complete.")

seqs <- as.character(dada2_seqs[,1])
names(seqs) <- paste0("'",dada2_otu[,1],"'")
library(phangorn)
library(msa)
library(DECIPHER)
print("Running Alignment")
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
print("Running phangorn")
phang.align <- as.phyDat(as(alignment, "matrix"), type="DNA")
print("Making distance matrix")
dm <- dist.ml(phang.align)
print("Making tree from distance matrix")
treeNJ <- NJ(dm)
print("Fitting")
fit <- pml(treeNJ, data=phang.align)
write.tree(fit$tree, "dada2_first_tree.tre")
save.image("save_tree.RData")
fitGTR <- update(fit, k=4, inv=0)
print("Optimizing Fit")
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
          rearrangement = "stochastic", control = pml.control(trace = 0))
write.tree(fitGTR$tree, "dada2_optim_tree.tre")

