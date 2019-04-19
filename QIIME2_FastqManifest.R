# Make a "fastq manifest" file:
# Assumes that all of your forward reads are in one folder and reverse reads are in the other folder.

fwd_path <- "/data/Users/kmccauley/PROSE_NEW/for_dada2/forward_reads"
rev_path <- "/data/Users/kmccauley/PROSE_NEW/for_dada2/reverse_reads"
fwd_reads <- sort(list.files(fwd_path, pattern=".1.fastq.gz", full.names=TRUE))
rev_reads <- sort(list.files(rev_path, pattern=".2.fastq.gz", full.names=TRUE))

# The problem is that this code works for *my* samples right now.. Not sure how to fix this...
sample.names <- sapply(strsplit(basename(fwd_reads), "[.]"), function(x) x[-length(x)])[1,]

fwd.test <- data.frame('sample-id'=sample.names,'absolute-filepath'= fwd_reads, direction="forward", check.names=F)
rev.test <- data.frame('sample-id'=sample.names,'absolute-filepath'= rev_reads, direction="reverse", check.names=F)
full.df <- rbind(fwd.test, rev.test)

write.csv(full.df, "fastq_manifest.csv", row.names=F, col.names=T)

