## Make faux fasta for Picrust

argv <- commandArgs(TRUE)
library(phyloseq)
library(data.table)
phy <- readRDS(argv[1])

tab <- data.table(V1=paste0(">", taxa_names(phy)), V2=data.frame(refseq(phy))[,1], order=1:length(taxa_names(phy)))
header <- data.frame(V1=tab$V1, V2=tab$order)
reads <- data.frame(V1=tab$V2, V2=tab$order)
comb <- rbind(header, reads)
comb <- comb[order(comb$V2),]
write.table(comb$V1, "dada2_fasta.fa", quote=F, row.names=F, sep="\t", col.names=F)
