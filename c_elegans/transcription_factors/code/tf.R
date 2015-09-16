################################################################################
##                                Dependencies                                ##
################################################################################

setwd("~/master_thesis/s_cerevisiae/transcription_factors/")
load(".RData")
save.image(".RData")

library(stringr)
library(ade4)

################################################################################
##                         Load yeast regulatory network                      ##
################################################################################

reg_interactions <- readLines("data/RegulationTwoColumnTable_Documented_2013927.tsv")
reg_interactions <- as.data.frame(str_match(reg_interactions, "^(.*);(.*)$")[,-1])

reg_net <- data.frame("GeneID" = reg_interactions[,2],
                      "TF" = reg_interactions[,1])
x <- reg_net$GeneID == "N/A"
reg_net$GeneID <- ifelse(x, NA, as.character(reg_net$GeneID))
reg_net <- reg_net[complete.cases(reg_net),]

################################################################################
##                           transform to binary matrix                       ##
################################################################################

ctable_reg_net <- table(reg_net) #as table
ctable_reg_net <- as.data.frame.matrix(ctable_reg_net, stringsAsFactors=TRUE) #as data.frame
ctable_reg_net <- data.matrix(ctable_reg_net) #as matrix
ctable_reg_net[ctable_reg_net > 1] <- 1
dim(ctable_reg_net)
################################################################################
##                        calculate binary similarities                       ##
################################################################################

dist.binary(ctable_reg_net[1:10,], method = 1, diag = FALSE, upper = FALSE)
tf_Jaccard <- dist.binary(ctable_reg_net, method = 1, diag = FALSE, upper = FALSE)
tf_Jaccard <- as.matrix(tf_Jaccard)
tf_Jaccard <- 1-tf_Jaccard

tf_Ochiai <- dist.binary(ctable_reg_net, method = 7, diag = FALSE, upper = FALSE)
tf_Ochiai <- as.matrix(tf_Ochiai)
tf_Ochiai <- 1-tf_Ochiai

################################################################################
##                           organismic net read node list                    ##
################################################################################

#read organismic network and extrac gene node names
organismicGOfiltered <- read.table("../Blast_Gold_standard/Results/Gen-Gen.edgelist_GS-GO_filtered.txt", col.names=c("Node1", "Node2"))

x <- as.data.frame(str_match(organismicGOfiltered$Node1, "^(.*);(.*$)")[,-1])
y <- as.data.frame(str_match(organismicGOfiltered$Node2, "^(.*);(.*$)")[,-1])

#geneID <- unique(c(as.character(x[,1]), as.character(y[,1])))

node1 <- x
node2 <- y
colnames(node1) <- c("GeneID", "EC")
colnames(node2) <- c("GeneID", "EC")

rm(x ,y)

################################################################################
##              assign weight to edges in organismicGOfiltered                ##
################################################################################

x <- node1$GeneID
y <- node2$GeneID

organismicGOfiltered$tf_jac <- rep(0, times=length(organismicGOfiltered$Node1))
pb <- txtProgressBar(min = 0, max = length(organismicGOfiltered$tf_jac), style = 3)
for(i in 1:length(organismicGOfiltered$tf_jac)) {
  organismicGOfiltered$tf_jac[i] <- tf_Jaccard[x[i], y[i]]
  setTxtProgressBar(pb, i)
}
close(pb)

organismicGOfiltered$tf_och <- rep(0, times=length(organismicGOfiltered$Node1))
pb <- txtProgressBar(min = 0, max = length(organismicGOfiltered$tf_och), style = 3)
for(i in 1:length(organismicGOfiltered$tf_och)) {
  organismicGOfiltered$tf_och[i] <- tf_Ochiai[x[i], y[i]]
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

sim <- as.matrix(tf_Ochiai[gogold$V2,gogold$V2])
dim(sim)
sim[is.na(sim)] <- 0

pdf("results/gold_sim_image_tf.pdf", width=4, height=4)
heatmap_2(sim, Rowv=F, Colv=F, do.dendro=c(F, F), col=grey(seq(0, 1, length = 256)), scale="none")
dev.off()

tiff("results/gold_sim_image_tf.tiff", width=4, height=4, units="in", res = 300)
heatmap_2(sim, Rowv=F, Colv=F, do.dendro=c(F, F), col=grey(seq(0, 1, length = 256)), scale="none")
dev.off()


