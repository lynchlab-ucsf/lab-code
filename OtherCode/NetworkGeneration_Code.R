rm(list=ls())
set.seed(123)
#devtools::install_github("zdk123/SpiecEasi")
pacman::p_load(WGCNA, SpiecEasi, phyloseq)
library(SpiecEasi)
library(phyloseq)
library(igraph)

setwd("/data/Users/kmccauley/MUPPITS/Manuscript/NetworkGeneration/")
nonrare <- read.table("/data/Users/kmccauley/MUPPITS/OTUtables/MUPPITS_nonrare_OTUtable_initial.txt", header=TRUE, check.names=F, comment="", sep="\t", row.names=1)
#Keep the taxonomy data for later, but separate it so that I have a table of counts
nonrare.tax <- cbind(rownames(nonrare), as.character(nonrare$taxonomy))
#Remove the character-based taxonomy variable
nonrare$taxonomy <- NULL

#Clean up nonrare.tax to implement into phyloseq:
taxanames <- strsplit(as.character(nonrare.tax[,2]),"; ")
mat <- t(sapply(taxanames,
                function(x,m) c(x,rep(NA,m-length(x))),
                max(rapply(taxanames,length))))

newnames <- substr(mat,4,35)
names(newnames) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
newnames <- as.matrix(newnames)
row.names(newnames) <- nonrare.tax[,1]

#Make my phyloseq object with an OTU table and taxonomy (no need for a tree or metadata)
phy <- phyloseq(otu_table(nonrare, taxa_are_rows=TRUE), tax_table(newnames))
phy

phy.filt <- filter_taxa(phy, function(x) sum(x >0) > (0.1*length(x)), TRUE)
#sparcc.nasal <- sparcc(t(otu_table(phy.filt)))
#save(sparcc.nasal, file="MUPPITSnasalSparCC.Rdata")
load("/data/Users/kmccauley/MUPPITS/NetworkAnalysis/MUPPITSnasalSparCC.Rdata")
sparcc.graph <- sparcc.nasal$Cor >= 0.5
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix::Matrix(sparcc.graph)

rownames(sparcc.graph) <- colnames(sparcc.graph) <- rownames(tax_table(phy.filt))
ig.mod <- graph.adjacency(sparcc.graph, mode='undirected')
min.ig.mod <- igraph::delete.vertices(simplify(ig.mod), degree(ig.mod)==0)
plot(min.ig.mod, vertex.color="green") # we can see all the attributes and weights
modules =cluster_fast_greedy(min.ig.mod)
print(modules)

lapply(communities(modules), function(x) tax_table(phy.filt)[match(x, rownames(tax_table(phy.filt)))])

subcolors <- colors()[!grepl("grey|gray",colors())]
mycolors <- sample(subcolors, 20)

## Move into a table of module-level information
OTUtab <- as.data.frame(otu_table(phy.filt))
modules_assigned <- flatten(data.frame(communities(modules)))
names(modules_assigned) <- paste0("Mod", 1:length(modules_assigned))

tax_desc <- lapply(1:length(modules_assigned), function(x) {
  tax.dat <- data.frame(tax_table(phy.filt)[match(modules_assigned[[x]], rownames(tax_table(phy.filt)))]@.Data)
  tax.dat$module <- names(modules_assigned)[x]
  tax.dat$otu_no <- rownames(tax.dat)
  return(tax.dat)
})

taxa_in_modules <- bind_rows(tax_desc)
write.csv(x=taxa_in_modules, file="taxa_in_modules.csv")

## Basically assigning the modules to the OTUs
OTUtab$modules <- unlist(lapply(rownames(OTUtab), function(x) {
  val <- names(modules_assigned)[sapply(modules_assigned, function(z) is.element(x, z))]
  if(length(val) == 0) val <- "Unassigned"
  return(val)
}))

newOTUmods <- OTUtab %>% 
  group_by(modules) %>% 
  summarise_all(sum) %>%
  column_to_rownames(.,"modules") %>% 
  t()

ModuleSummary <- OTUtab %>% 
  group_by(modules) %>% 
  summarize_if(is.numeric, sum) %>%
  column_to_rownames("modules") %>% 
  t() %>% 
  colSums() %>% 
  sort()
ModuleSummary

write.csv(x = newOTUmods, file="ModuleAggregated_OTU_SparCCOnly_seed123.csv", row.names=TRUE)

