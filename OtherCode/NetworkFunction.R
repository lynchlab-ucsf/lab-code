make_network <- function(phyloseq, cutoff=0.25, fig.name="") {
    set.seed(123)
require(SpiecEasi)
require(tidyverse)
require(igraph)
require(RColorBrewer)
    
phyloseq <- filter_taxa(phyloseq, function(x) sum(x) > 0, TRUE)
results <- list()
print("Making network ...")
sparcc_res <- sparcc(t(otu_table(phyloseq)))
print("Network Complete! Preparing results ...")
sparcc_graph <- sparcc_res$Cor >= cutoff
diag(sparcc_graph) <- 0
sparcc_graph <- Matrix::Matrix(sparcc_graph)
rownames(sparcc_graph) <- colnames(sparcc_graph) <- taxa_names(phyloseq)
ig.mod <- graph.adjacency(sparcc_graph, mode='undirected')
V(ig.mod)$value <- tax_table(phyloseq)@.Data[,"Genus"]
min.ig.mod <- igraph::delete.vertices(simplify(ig.mod), degree(ig.mod)==0)

results$ig <- min.ig.mod
comm = cluster_fast_greedy(min.ig.mod)
results$comm <- comm
print(modularity(comm))

results$communities <- lapply(communities(comm), function(x) tax_table(phyloseq)[match(x, rownames(tax_table(phyloseq)))])

OTUtab <- as.data.frame(otu_table(phyloseq))
comm_assigned <- flatten(data.frame(communities(comm)))
names(comm_assigned) <- paste0("Mod", 1:length(comm_assigned))

tax_desc <- lapply(1:length(comm_assigned), function(x) {
  tax.dat <- data.frame(tax_table(phyloseq)[match(comm_assigned[[x]], rownames(tax_table(phyloseq)))]@.Data)
  tax.dat$module <- names(comm_assigned)[x]
  tax.dat$otu_no <- rownames(tax.dat)
  return(tax.dat)
})

results$taxa_in_comm <- dplyr::bind_rows(tax_desc)
#write.csv(x=taxa_in_comm, file="taxa_in_comm.csv")

## Basically assigning the comm to the OTUs
OTUtab$comm <- unlist(lapply(rownames(OTUtab), function(x) {
  val <- names(comm_assigned)[sapply(comm_assigned, function(z) is.element(x, z))]
  if(length(val) == 0) val <- "Unassigned"
  return(val)
}))

results$table <- OTUtab %>% 
  group_by(comm) %>% 
  summarise_all(sum) %>%
  tibble::column_to_rownames(.,"comm") %>% 
  t()
                                 
V(min.ig.mod)$size=4
V(min.ig.mod)$label=""
coul <- c(brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))
label.list <- names(sort(table(as.factor(V(min.ig.mod)$value)), decreasing=TRUE))[2:17]
my_color <- coul[as.numeric(factor(V(min.ig.mod)$value, levels=label.list))]

l <- layout_with_fr(min.ig.mod, niter=50)
results$plot <- plot(min.ig.mod, layout=l, vertex.color=my_color)
results$plot <- legend("bottomleft", legend=label.list, fill = coul, inset = c(-0.1, -0.1),pt.cex=3, cex=0.6)
results$plot <- title(fig.name)
return(results)
print("...done!")
}
                                 
#init_network_test <- make_network(wisc2mo_filt)