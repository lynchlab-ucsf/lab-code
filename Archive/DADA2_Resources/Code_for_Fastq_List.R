# This is a code to generate a navigator for the DADA2 pipeline. It basically will tell DADA2 where to find the sample data it is looking for.
# If you provide an argument, it will be a NextSeq Run ID (ie, 190826_NS500170_0076_AH772CBGXC). This will pull only one run, and will pull the information from the "Processed" directory. From there, you can remove any sample data that you don't want to consider, though I am considering adding a "valid_states" option, much like QIIME's filter_fasta.py.
## Katie McCauley
## Last Modified September 4, 2019

rm(list=ls())
pacman::p_load(optparse)

option_list <- list(make_option(c("-r", "--run"), type="character", default=NULL, help="Provide the run name", metavar="character"),
                    make_option(c("-m", "--map"), type="character", default=NULL, help="Map file for subsetting the files", metavar="character"),
                    make_option(c("-s", "--valid_states"), type="character", default=NULL, help="Valid states for subsetting", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#Arguments:
## run_name (give the name of the run present and unzipped in the processed folder
## map (give the name of the map file that you might want to use to filter samples from this file
## valid_states (give the variable and the value that you want to keep) -- see filter_fasta.py for examples

if(is.null(opt$run)) {
path <- getwd()

a <- list.dirs(".",  recursive=FALSE, full.names=TRUE)
directories <- a[grepl("NS500170", a) & !grepl("script_pipeline_output", a)]
print(paste("R sees", length(directories), "runs in this folder"))
} else {

path <- paste0("/data/NextSeq_data/NextSeq_Processed/", opt$run)
directories <- path
}
locations <- sort(list.files(paste0(directories,"/submission"), pattern="^R.",full.names=TRUE))
print(paste("R will obtain fastq files from these locations:"))
print(locations)

if(!is.null(opt$map)) {
if(grepl(".csv", opt$map)) sampmap <- read.csv(paste0(directories,"/", opt$map), header=TRUE, check.names=F)
if(grepl(".txt", opt$map)) sampmap <- read.table(paste0(directories,"/", opt$map), header=TRUE, check.names=F, sep="\t", comment="")
var <- strsplit(opt$valid_states, ":")
subsetmap <- sampmap[sampmap[,var[[1]][1]] %in% var[[1]][2],]
print(head(subsetmap))
samplename.loc <- which(grepl("sampleid", names(sampmap), ignore.case=TRUE))
print(samplename.loc)
initfastqs <- list.files(locations, pattern="[.]fastq", full.names=TRUE)
fastqs <- initfastqs[grepl(paste(subsetmap[,samplename.loc], collapse="|"), initfastqs)]
print(head(fastqs))
} else {
fastqs <- list.files(locations, pattern="[.]fastq", full.names=TRUE)
}

sample.names <- sapply(strsplit(basename(fastqs), "[_]"), function(x) x[1]) 
# My sample names were kinda complicated, so I used the above code to extract the sample names and create cleaner sample name files for those that get filtered. This can also be an incorrect file name that is corrected after the fastq_file_list is printed.

run.source <- sapply(strsplit(fastqs, "[/]"), function(x) x[2])

comb <- data.frame(fastqs, sample.names, run.source)

comb$direction <- ifelse(grepl("R1", comb$fastqs), "R1","R2")
print("Sample Totals:")
table(comb$direction,useNA="a")

if(sum(table(comb$sample.names) < 2)>0) {
print(paste("Paired R1/R2 files are not available for", sum(table(comb$sample.names) < 2)," samples."))
print(names(table(comb$sample.names))[table(comb$sample.names) < 2])
print("Carefully consider if this is expected before proceeding.")
}
if(sum(table(comb$sample.names) > 2)>0) {
print("Some samples have more than two files associated with them.")
samps <- names(table(comb$sample.names))[table(comb$sample.names) > 2]
print(paste("This is the case for ", length(samps), "samples."))
print(comb[comb$sample.names %in% samps,])
print("Carefully consider if this is expected before proceeding.")
}
write.csv(comb, "fastq_file_list.csv",row.names=F)
