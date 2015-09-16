################################################################################
##                                Dependencies                                ##
################################################################################
#create uniq list of arabidopsis genes in metabolic and microarray sim-matrix
#cat ../Blast_Gold_standard/Results/genes_in_adjacency_matrix.list  gene_names21122.csv |sort|uniq -d
#wc of list: 8211
source("http://bioconductor.org/biocLite.R")
biocLite("GOSemSim")
biocLite("org.At.tair.db")
#installed "GO.db" dependency
library(GO.db)
library(GOSemSim)
library(stringr)

ls("package:GOSemSim")
ls("package:GO.db")


setwd("~/master_thesis/Arabidopsis/GO_similarity_mf/")
oldpar <-par()
load(".RData")
save.image(".RData")
################################################################################
##                           organismic net read node list                    ##
################################################################################

#read organismic network and extrac gene node names
organismicGOfiltered <- read.table("../Blast_Gold_standard/Results/Gen-Gen.edgelist_GS-GO_filtered.txt", col.names=c("Node1", "Node2"))

x <- as.data.frame(str_match(organismicGOfiltered$Node1, "^(.{9});(.*$)")[,-1])
y <- as.data.frame(str_match(organismicGOfiltered$Node2, "^(.{9});(.*$)")[,-1])

geneID <- unique(c(as.character(x[,1]), as.character(y[,1])))

node1 <- x
node2 <- y
colnames(node1) <- c("GeneID", "EC")
colnames(node2) <- c("GeneID", "EC")

rm(x ,y)
################################################################################
##              calculate similarity between genes for GO terms               ##
################################################################################
#geneSim: Given two genes, 
#this function will calculate the semantic similarity between them, 
#and return their semantic similarity and the corresponding GO terms

molecular_function_GO <- mgeneSim(genes=geneID, ont="MF", 
                                  organism="arabidopsis", measure="Wang")

cellular_component_GO <- mgeneSim(genes=geneID, ont="CC", 
                                   organism="arabidopsis", measure="Wang")

################################################################################
##              assign weight to edges in organismicGOfiltered                ##
################################################################################

x <- as.character(node1$GeneID)
y <- as.character(node2$GeneID)
z <- as.character(rownames(molecular_function_GO))

organismicGOfiltered$go_mf <- rep(0, times=length(organismicGOfiltered$Node1))
pb <- txtProgressBar(min = 0, max = length(organismicGOfiltered$go_mf), style = 3)
for(i in 1:length(organismicGOfiltered$go_mf)) {
  ifelse(x[i]%in%z & y[i]%in%z, 
         organismicGOfiltered$go_mf[i] <- molecular_function_GO[x[i], y[i]],
         organismicGOfiltered$go_mf[i] <- NA)
  setTxtProgressBar(pb, i)
}
close(pb)

write.table(organismicGOfiltered, "results/edgelist_weighted.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")



