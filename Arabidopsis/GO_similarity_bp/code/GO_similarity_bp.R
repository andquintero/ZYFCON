################################################################################
##                                Dependencies                                ##
################################################################################
#create uniq list of arabidopsis genes in metabolic and microarray sim-matrix
#cat ../Blast_Gold_standard/Results/genes_in_adjacency_matrix.list  gene_names21122.csv |sort|uniq -d
#wc of list: 8211
source("http://bioconductor.org/biocLite.R")
biocLite("GOSemSim")
#installed "GO.db" dependency
library(GO.db)
library(GOSemSim)
library(stringr)

ls("package:GOSemSim")
ls("package:GO.db")


setwd("~/master_thesis/Arabidopsis/GO_similarity_bp//")
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

biological_process_GO <- mgeneSim(genes=geneID, ont="BP", 
                                  organism="arabidopsis", measure="Wang")

cellular_component_GO <- mgeneSim(genes=geneID, ont="CC", 
                                   organism="arabidopsis", measure="Wang")

################################################################################
##              assign weight to edges in organismicGOfiltered                ##
################################################################################

x <- as.character(node1$GeneID)
y <- as.character(node2$GeneID)
z <- as.character(rownames(biological_process_GO))

organismicGOfiltered$go_bp <- rep(0, times=length(organismicGOfiltered$Node1))
pb <- txtProgressBar(min = 0, max = length(organismicGOfiltered$go_bp), style = 3)
for(i in 1:length(organismicGOfiltered$go_bp)) {
  ifelse(x[i]%in%z & y[i]%in%z, 
         organismicGOfiltered$go_bp[i] <- biological_process_GO[x[i], y[i]],
         organismicGOfiltered$go_bp[i] <- NA)
  setTxtProgressBar(pb, i)
}
close(pb)

write.table(organismicGOfiltered, "results/edgelist_weighted.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")









write.table(MFm,file="arabidopsis_MF.txt")
image(MFm, axes = FALSE, col = grey(seq(0, 1, length = 256)))

MF<-mgeneSim(genes=gen_names, ont="MF", organism="arabidopsis", measure="Wang")
BP<-mgeneSim(genes=gen_names, ont="BP", organism="arabidopsis", measure="Wang")
CC<-mgeneSim(genes=gen_names, ont="CC", organism="arabidopsis", measure="Wang")
par(mfrow=c(1,3))
image(MF, axes = FALSE, col = grey(seq(1, 0, length = 256)), xlab="MF")
image(BP, axes = FALSE, col = grey(seq(1, 0, length = 256)), xlab="BP")
image(CC, axes = FALSE, col = grey(seq(1, 0, length = 256)), xlab="CC")

write.table(BP,file="arabidopsis_BP.txt", sep="\t")

image((t(MF)[,nrow(MF):1]), axes = FALSE, col = grey(seq(1, 0, length = 256)), xlab="MF")
image((t(BP)[,nrow(BP):1]), axes = FALSE, col = grey(seq(1, 0, length = 256)), xlab="BP")
image((t(CC)[,nrow(CC):1]), axes = FALSE, col = grey(seq(1, 0, length = 256)), xlab="CC")

BPbig<-read.table("arabidopsis_MF.txt")
write.table(MF,file="arabidopsis_MF.txt", sep="\t")

