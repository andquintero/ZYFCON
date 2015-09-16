################################################################################
##                                Dependencies                                ##
################################################################################

setwd("~/master_thesis/Arabidopsis/expression/")
load(".RData")
save.image(".RData")

library(stringr)
################################################################################
##                           read big similarity matrix                       ##
################################################################################

expression_sim <- read.csv("data/S.csv")
expression_sim <- expression_sim + t(expression_sim)
diag(expression_sim) <- 1

sim_geneID <- readLines("data/gene_names21122.csv")

rownames(expression_sim)  <- sim_geneID
colnames(expression_sim)  <- sim_geneID

################################################################################
##                           organismic net read node list                    ##
################################################################################

#read organismic network and extrac gene node names
organismicGOfiltered <- read.table("../Blast_Gold_standard/Results/Gen-Gen.edgelist_GS-GO_filtered.txt", col.names=c("Node1", "Node2"))

x <- as.data.frame(str_match(organismicGOfiltered$Node1, "^(.{9});(.*$)")[,-1])
y <- as.data.frame(str_match(organismicGOfiltered$Node2, "^(.{9});(.*$)")[,-1])

node1 <- x
node2 <- y
colnames(node1) <- c("GeneID", "EC")
colnames(node2) <- c("GeneID", "EC")

rm(x ,y)

################################################################################
##              assign weight to edges in organismicGOfiltered                ##
################################################################################

x <- as.character(node1$GeneID)
y <- as.character(node2$GeneID)
z <- as.character(rownames(expression_sim))

organismicGOfiltered$expr <- rep(0, times=length(organismicGOfiltered$Node1))
pb <- txtProgressBar(min = 0, max = length(organismicGOfiltered$expr), style = 3)
for(i in 1:length(organismicGOfiltered$expr)) {
  ifelse(x[i]%in%z & y[i]%in%z, 
         organismicGOfiltered$expr[i] <- expression_sim[x[i], y[i]],
         organismicGOfiltered$expr[i] <- NA)
  setTxtProgressBar(pb, i)
}
close(pb)

write.table(organismicGOfiltered, "results/edgelist_weighted.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
