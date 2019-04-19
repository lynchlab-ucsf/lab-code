#This code takes a dataset and prepares it for use in QIIME
#Make QIIME file

make.qiime <- function(inputfile,outputfile,idvar) {

intable <- read.table(inputfile, header=T, sep="\t",comment="", check.names=F)
#intable$BarcodeSequence <- "CCGTGACAACTC"
#intable$LinkerPrimerSequence <- "GTGTGCCAGCMGCCGCGGTAA"
intable$Description <- intable[,idvar]

init.cols <- c(idvar)
last.cols <- c("Description")
mid.cols <- names(intable)[-which(names(intable) %in% c(init.cols,last.cols))]



towrite <- intable[,c(init.cols,names(intable)[-which(names(intable) %in% c(init.cols,last.cols))],last.cols)]
names(towrite)[names(towrite) %in% idvar] <- "#SampleID"
write.table(towrite, outputfile, sep="\t",row.names=F,quote=TRUE)
}

make.qiime("/data/Users/kmccauley/PROSE_NEW/KMcCauley01/dat/der/FINAL_MAP_FILE_subset_colds_RVprop.txt","/data/Users/kmccauley/PROSE_NEW/KMcCauley01/dat/der/FINAL_MAP_FILE_subset_colds_RVprop2.txt", idvar="#SampleID")
