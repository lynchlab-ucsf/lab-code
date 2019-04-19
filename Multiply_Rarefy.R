#############################################################################################################
######  This is a script to calculate a representative rarefied OTU table from an unrarefied OTU table, ######
######  an alterative to single rarefied tables that stabilizes the random sampling and results         ######
######  in a rarefied table that may be more consistent with the original data than a single            ######  	     
######  rarefied table. Briefly, many single-rarefied OTU tables are calculated, and the distance       ######
######  between the subject-specific rarefied vectors is calculated.  The rarefied vector that is       ######
######  the minimum average (or median) distance from itself to all other rarefied vectors              ######
######  is considered the most representative for that subject and built into the new                   ######
######  rarefied table.  This work is currently in progress by HFHS investigators                       ######
######  (Sitarik A, Levin A, Havstad S) and UCSF investigators (Fujimura K, Lynch S, Faruqi A).         ######
##############################################################################################################


#################################

# library(GUniFrac) # (don't need the package if you call Rarefy below)
library(vegan)

# (Rarefy function from GUniFrac package)
Rarefy <- function (otu.tab, depth = min(rowSums(otu.tab))) 
{
    otu.tab <- as.matrix(otu.tab)
    ind <- (rowSums(otu.tab) < depth)
    sam.discard <- rownames(otu.tab)[ind]
    otu.tab <- otu.tab[!ind, ]
    rarefy <- function(x, depth) {
        y <- sample(rep(1:length(x), x), depth)
        y.tab <- table(y)
        z <- numeric(length(x))
        z[as.numeric(names(y.tab))] <- y.tab
        z
    }
    otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
    rownames(otu.tab.rff) <- rownames(otu.tab)
    colnames(otu.tab.rff) <- colnames(otu.tab)
    return(list(otu.tab.rff = otu.tab.rff, discard = sam.discard))
}


################################################################################
###############			Parameters 		################################
################################################################################

# specify the raw OTU count table, with samples = rows, taxa = columns
# rawtab = otu_tab_t

# specify the depth you would like to rarefy your tables to
# the default is to just use the minimum sequencing depth
# raredepth = min(rowSums(rawtab))

# specify the number of rarefied tables you would like to generate 
# to calculate your representatiave rarefied table from 
# ntables = 100

# specify the distance measure to use to calculate distance between rarefied data sets, for each subject
# can be any of the methods available in the vegdist function of vegan
# distmethod = "euclidean"

# specify the method to summarize across distances
# if mean distance, then summarymeasure = mean
# if median distance, then summarymeasure = median
# summarymeasure = mean

# specify the seed start for the rarefied tables
# for each subsequent table, 1 will be added that the previous seed
# for reproducibility, always save your seedstart value (or just use the default for simplicity).
# seedstart = 500

# specify if you want progress updates to be printed
# verbose = TRUE

### returns a representative rarefied OTU table of class matrix. 

################################################################################
################################################################################
################################################################################

reprare <- function(rawtab=otu_tab_t, raredepth = min(rowSums(otu_tab_t)), ntables=100, distmethod="euclidean", 
summarymeasure=mean, seedstart=500, verbose=TRUE) {

raretabs = list()
for (z in 1:ntables) {
if (verbose==TRUE) {
print(paste("calculating rarefied table number", z, sep=" "))
}
set.seed(seedstart + z)
raretabs[[z]] = Rarefy(rawtab, depth = raredepth)[[1]]
}

raretabsa = array(unlist(raretabs), dim = c(nrow(raretabs[[z]]), ncol(rawtab), ntables))

final_tab = c()
for (y in 1:nrow(raretabs[[z]])) {
if (verbose==TRUE) {
print(paste("determining rep rarefied vector for subject number", y, sep=" "))
}
distmat = as.matrix(vegdist(t(raretabsa[y,,]), method=distmethod)) # distance across reps for subject y
distsummary = apply(distmat, 2, summarymeasure) 
whichbestrep = which(distsummary == min(distsummary))[1]  # the best rep is the one with the minimum average/median distance to all other reps. (in case of ties, just select the first)
bestrep = raretabsa[y,,whichbestrep] # select that rep only for subject y
final_tab = rbind(final_tab, bestrep) # build that rep for subject y into final table
}
rownames(final_tab) = rownames(raretabs[[z]])
colnames(final_tab) = colnames(rawtab)

return(final_tab)
}



###### example runs of the function:  ######
runexample=FALSE
if (runexample==TRUE) {

### dummy data set for example ###
ntaxa = 200
nsubj = 50
set.seed(444)
dummyOTU <- matrix(sample(0:500, ntaxa*nsubj, prob=c(0.7,0.1,0.1,rep(0.1/498, 498)), replace=TRUE), ncol=ntaxa)
colnames(dummyOTU) = paste("OTU", 1:ntaxa, sep="")
rownames(dummyOTU) = paste("subj", 1:nsubj, sep="")
sort(rowSums(dummyOTU)) # sequencing depth is uneven

# specify the minimum depth
repraretable = reprare(rawtab=dummyOTU, raredepth=min(rowSums(dummyOTU)), ntables=100, distmethod="euclidean", 
summarymeasure=mean, seedstart=500, verbose=TRUE)
dim(repraretable)
sort(rowSums(repraretable)) # sequencing depth is now even

# specify a depth other than the minimum
repraretable = reprare(rawtab=dummyOTU, raredepth=3380, ntables=100, distmethod="euclidean", 
summarymeasure=mean, seedstart=500, verbose=TRUE)
dim(repraretable) # subjects with less than the minimum are no longer in the table
sort(rowSums(repraretable)) # sequencing depth is now even

}

#repraretable = reprare(rawtab=closed.array, raredepth=min(rowSums(closed.array)), ntables=100, distmethod="euclidean", 
                       #summarymeasure=mean, seedstart=500, verbose=TRUE)

#rep.array<-as.data.frame(t(repraretable))

setwd("/data/Users/kmccauley/PROSE_OTUPICK")
MYdata=read.table("filtered_otutable_finalfilt_0_001_tax.txt",sep="\t",header=T,comment="",check.names=F)
taxonomy <- MYdata[,c("#OTU ID","taxonomy")]
rownames(MYdata) <- MYdata[,"#OTU ID"]
MYdata[,c("#OTU ID")] <- NULL
MYdata$taxonomy <- NULL

MYarray=as.data.frame(t(MYdata))


repraretable <- reprare(rawtab=MYarray, raredepth=2000, ntables=100, distmethod="bray", 
                       summarymeasure=mean, seedstart=123, verbose=TRUE)

rep.array<-as.data.frame(t(repraretable))


dim(rep.array)
rep.array[,"#OTU ID"] <- row.names(rep.array)

rarefied2 <- merge(rep.array, taxonomy, by=c("#OTU ID"))
write.table(rarefied2,file="/data/Users/kmccauley/PROSE_OTUPICK/multiplyrarefy/rarefied_filtered_otutable.txt",sep="\t",quote=F,row.names=FALSE)

