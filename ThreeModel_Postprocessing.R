# Post-processing of CSV file from "ThreeModel_LabCode" script.
## It will extract all OTUs with p-value less than 0.05 (because you may need to decide on your preferred q-value cut-off), and keep only columns that are important
### The make_data function at the bottom takes in the file path to the "full" CSV file, and the second option is the name and location of the "pruned" CSV file

make_data <- function(csvfile,outfile) {
   read.data <- read.csv(csvfile)
   sig.data <- read.data[read.data$best.pval < 0.05,]
   keep.vars <- c("OTUname","pois.pval","nb.pval","zinb.pval","best.mod","qval.best","mean_diff","zerotrt1","zerotrt2","nonzerotrt1","nonzerotrt2","totaltrt1","totaltrt2","taxonomy")
   fin.data <- sig.data[keep.vars]
   fin.data <- fin.data[order(fin.data$qval.best),]
   write.csv(fin.data,outfile)
}

make_data("/data/Users/kmccauley/PROSE_NEW/AnalysisData/TESTING_ThreeModel.csv","/data/Users/kmccauley/PROSE_NEW/AnalysisData/TESTING_PROCESSED_ThreeModel.csv")
