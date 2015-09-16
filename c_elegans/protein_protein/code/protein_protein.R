################################################################################
##                                Dependencies                                ##
################################################################################

setwd("~/master_thesis/s_cerevisiae/protein_protein/")
load(".RData")
save.image(".RData")

library(stringr)
library(ade4)

################################################################################
##                           organismic net read node list                    ##
################################################################################

#read organismic network and extrac gene node names
organismicGOfiltered <- read.table("../Blast_Gold_standard/Results/Gen-Gen.edgelist_GS-GO_filtered.txt", col.names=c("Node1", "Node2"))

x <- as.data.frame(str_match(organismicGOfiltered$Node1, "^(.*);(.*$)")[,-1])
y <- as.data.frame(str_match(organismicGOfiltered$Node2, "^(.*);(.*$)")[,-1])

geneID <- unique(c(as.character(x[,1]), as.character(y[,1])))

node1 <- x
node2 <- y
colnames(node1) <- c("GeneID", "EC")
colnames(node2) <- c("GeneID", "EC")

rm(x ,y)

################################################################################
##                      read interactions and mapping files                   ##
################################################################################
#interactions
ppi <- read.table("data/mapped_ppi.txt", col.names=c("A", "B"))


#keep in column A only genes in the network
ppi2 <- data.frame("A" = c(as.character(ppi$A), as.character(ppi$B)),
                   "B" = c(as.character(ppi$B), as.character(ppi$A)))
x <- ppi2$A%in%geneID
ppi2 <- ppi2[x,]

################################################################################
##                           transform to binary matrix                       ##
################################################################################

ctable_ppi <- table(ppi2) #as table
ctable_ppi <- as.data.frame.matrix(ctable_ppi, stringsAsFactors=TRUE) #as data.frame
ctable_ppi <- data.matrix(ctable_ppi) #as matrix
ctable_ppi[ctable_ppi > 1] <- 1
dim(ctable_ppi)
ctable_ppi <- ctable_ppi[-1,-1]
diag(ctable_ppi) <- 1

#filter to keep only genes in organismic network
x <- rownames(ctable_ppi)%in%geneID
ctable_ppi <- ctable_ppi[x,]
dim(ctable_ppi)

################################################################################
##                        calculate binary similarities                       ##
################################################################################

dist.binary(ctable_ppi[1:10,], method = 1, diag = FALSE, upper = FALSE)
ppi_Jaccard <- dist.binary(ctable_ppi, method = 1, diag = FALSE, upper = FALSE)
ppi_Jaccard <- as.matrix(ppi_Jaccard)
ppi_Jaccard <- 1-ppi_Jaccard

ppi_Ochiai <- dist.binary(ctable_ppi, method = 7, diag = FALSE, upper = FALSE)
ppi_Ochiai <- as.matrix(ppi_Ochiai)
ppi_Ochiai <- 1-ppi_Ochiai

x <- rownames(ppi_Ochiai)
pb <- txtProgressBar(min = 0, max = length(x), style = 3)
for(i in 1:length(x)) {
  for(j in i:length(x)) {
    if(ctable_ppi[x[i],x[j]] > 0) {
      ppi_Ochiai[x[i],x[j]]  <- 1
      ppi_Ochiai[x[j],x[i]]  <- 1
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

################################################################################
##              assign weight to edges in organismicGOfiltered                ##
################################################################################

x <- as.character(node1$GeneID)
y <- as.character(node2$GeneID)
z <- as.character(rownames(ppi_Ochiai))

organismicGOfiltered$ppi_jac <- rep(0, times=length(organismicGOfiltered$Node1))
pb <- txtProgressBar(min = 0, max = length(organismicGOfiltered$ppi_jac), style = 3)
for(i in 1:length(organismicGOfiltered$ppi_jac)) {
  ifelse(x[i]%in%z & y[i]%in%z, 
         organismicGOfiltered$ppi_jac[i] <- ppi_Jaccard[x[i], y[i]],
         organismicGOfiltered$ppi_jac[i] <- NA)
  setTxtProgressBar(pb, i)
}
close(pb)

organismicGOfiltered$ppi_och <- rep(0, times=length(organismicGOfiltered$Node1))
pb <- txtProgressBar(min = 0, max = length(organismicGOfiltered$ppi_och), style = 3)
for(i in 1:length(organismicGOfiltered$ppi_och)) {
  ifelse(x[i]%in%z & y[i]%in%z, 
         organismicGOfiltered$ppi_och[i] <- ppi_Ochiai[x[i], y[i]],
         organismicGOfiltered$ppi_och[i] <- NA)
  setTxtProgressBar(pb, i)
}
close(pb)

write.table(organismicGOfiltered, "results/edgelist_weighted.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


################################################################################
##                               image similarity matrix                      ##
################################################################################
source("http://bioconductor.org/biocLite.R")
biocLite("Heatplus")
library(Heatplus)

#read gold standard genesIDs and its EC 
gogold <- read.table("../Gold_Standard/data/yeast_GS_EC-gene_GO.txt")

sim <- as.matrix(ppi_Ochiai[gogold$V2,gogold$V2])
dim(sim)
sim[is.na(sim)] <- 0

pdf("results/gold_sim_image_ppi.pdf", width=4, height=4)
heatmap_2(sim, Rowv=F, Colv=F, do.dendro=c(F, F), col=grey(seq(0, 1, length = 256)), scale="none")
dev.off()

tiff("results/gold_sim_image_ppi.tiff", width=4, height=4, units="in", res = 300)
heatmap_2(sim, Rowv=F, Colv=F, do.dendro=c(F, F), col=grey(seq(0, 1, length = 256)), scale="none")
dev.off()

heatmap_2(sim, col=grey(seq(0, 1, length = 256)), scale="none")
