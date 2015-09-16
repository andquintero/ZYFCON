################################################################################
##                                Dependencies                                ##
################################################################################

setwd("~/master_thesis/s_cerevisiae/integration/")
load(".RData")
save.image(".RData")

library(stringr)

################################################################################
##                           organismic net read node list                    ##
################################################################################

#read organismic network and extract gene node names
organismicGOfiltered <- read.table("../Blast_Gold_standard/Results/Gen-Gen.edgelist_GS-GO_filtered.txt", col.names=c("Node1", "Node2"))

################################################################################
##                        organismic net read edge weights                    ##
################################################################################

#expression
x <- read.table("../expression/results/edgelist_weighted.txt", header=TRUE)
organismicGOfiltered$expr <- x$expr
organismicGOfiltered$expr[organismicGOfiltered$expr > 1] <- 1

#GO molecular function
x <- read.table("../GO_similarity_mf/results/edgelist_weighted.txt", header=TRUE)
organismicGOfiltered$go_mf <- x$go_mf

#GO biological process
x <- read.table("../GO_similarity_bp/results/edgelist_weighted.txt", header=TRUE)
organismicGOfiltered$go_bp <- x$go_bp

#transcription factor
x <- read.table("../transcription_factors/results/edgelist_weighted.txt", header=TRUE)
organismicGOfiltered$tf_och <- x$tf_och

#pfam
x <- read.table("../pfam/results/edgelist_weighted.txt", header=TRUE)
organismicGOfiltered$pfam_och <- x$pfam_och

#protein-protein interaction
x <- read.table("../protein_protein/results/edgelist_weighted.txt", header=TRUE)
organismicGOfiltered$ppi_och <- x$ppi_och

rm(x)

################################################################################
##                           read gold standard geneID and ec                 ##
################################################################################

#read gold standard genesIDs and its EC 
gogold <- read.table("../Gold_Standard/data/yeast_GS_EC-gene_GO.txt")
gold_geneID_ec <- paste(gogold$V2,";",gogold$V1, sep="")


################################################################################
##                           subset gold net and putative net                 ##
################################################################################

#identify gold standards in organismic net
x <- organismicGOfiltered$Node1%in%gold_geneID_ec
y <- organismicGOfiltered$Node2%in%gold_geneID_ec

gold <- data.frame("Node1" = ifelse(x, 1, 0),
                   "Node2" = ifelse(y, 1, 0))
rm(x,y)

x <- gold$Node1 == 1 & gold$Node2 == 1

edges_gold <- organismicGOfiltered[x,]
edges_putative <- organismicGOfiltered[!x,]

################################################################################
##                           empirical cumulative distributions               ##
################################################################################

gecdf_expr <- ecdf(edges_gold$expr)
gecdf_go_mf <- ecdf(edges_gold$go_mf)
gecdf_go_bp <- ecdf(edges_gold$go_bp)
gecdf_tf_och <- ecdf(edges_gold$tf_och)
gecdf_pfam_och <- ecdf(edges_gold$pfam_och)
gecdf_ppi_och <- ecdf(edges_gold$ppi_och)

plot(gecdf_expr)
plot(gecdf_go_bp)
plot(gecdf_go_mf)
plot(gecdf_tf_och)
plot(gecdf_pfam_och)
plot(gecdf_ppi_och)
 
################################################################################
##     fuction to calculate ks.test p.val in all rows agaist a ecdf           ##
################################################################################

edge_pval <- function(data, gecdf) {
  #data = vector to apply
  #gecdf= empirical cumulative distribution fuction
  pb <- txtProgressBar(min = 0, max = length(data), style = 3)
  #x <- rep(NA, times=length(data[,col]))
  x <- c()
  for(i in 1:length(data)) {
    ifelse(is.na(data[i]),
           x[i] <- NA,
           x[i] <- ks.test(data[i], gecdf, alternative="less")$p.value)
    setTxtProgressBar(pb, i)
  }
  return(x) 
  close(pb)
}


################################################################################
##               calcule ks.test p.val to all edges and data types            ##
################################################################################

#create a dataframe to save results of test
ks_pval_gold_ecdf <- organismicGOfiltered
ks_pval_gold_ecdf[,3:length(colnames(ks_pval_gold_ecdf))] <- NA

# #BAD calcule ks.test p-values BAD
# ks_pval_gold_ecdf$expr <- edge_pval(data=organismicGOfiltered$expr, gecdf=gecdf_expr)
# ks_pval_gold_ecdf$go_mf <- edge_pval(data=organismicGOfiltered$go_mf, gecdf=gecdf_go_mf)
# ks_pval_gold_ecdf$go_bp <- edge_pval(data=organismicGOfiltered$go_bp, gecdf=gecdf_go_bp)
# ks_pval_gold_ecdf$tf_och <- edge_pval(data=organismicGOfiltered$tf_och, gecdf=gecdf_tf_och)
# ks_pval_gold_ecdf$pfam_och <- edge_pval(data=organismicGOfiltered$pfam_och, gecdf=gecdf_pfam_och)
# ks_pval_gold_ecdf$ppi_och <- edge_pval(data=organismicGOfiltered$ppi_och, gecdf=gecdf_ppi_och)

ks_pval_gold_ecdf$expr <- gecdf_expr(organismicGOfiltered$expr)
ks_pval_gold_ecdf$go_mf <- gecdf_go_mf(organismicGOfiltered$go_mf)
ks_pval_gold_ecdf$go_bp <- gecdf_go_bp(organismicGOfiltered$go_bp)
ks_pval_gold_ecdf$tf_och <- gecdf_tf_och(organismicGOfiltered$tf_och)
ks_pval_gold_ecdf$pfam_och <- gecdf_pfam_och(organismicGOfiltered$pfam_och)
ks_pval_gold_ecdf$ppi_och <- gecdf_ppi_och(organismicGOfiltered$ppi_och)

write.table(ks_pval_gold_ecdf, "results/edgelist_weighted_all.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

################################################################################
##                                   Write WOENs                              ##
################################################################################

colnames(ks_pval_gold_ecdf)
#exprs
write.table(ks_pval_gold_ecdf[,c(1,2,3)], "results/edgelist_weighted_exprs.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#gomf
write.table(ks_pval_gold_ecdf[,c(1,2,4)], "results/edgelist_weighted_gomf.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#gobp
write.table(ks_pval_gold_ecdf[,c(1,2,5)], "results/edgelist_weighted_gobp.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#tf
write.table(ks_pval_gold_ecdf[,c(1,2,6)], "results/edgelist_weighted_tf.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#pfam
write.table(ks_pval_gold_ecdf[,c(1,2,7)], "results/edgelist_weighted_pfam.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
#ppi
write.table(ks_pval_gold_ecdf[,c(1,2,8)], "results/edgelist_weighted_ppi.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


################################################################################
##                                Fisher p.value sum                          ##
################################################################################

fishersmethod <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
  return(p.val)
}

edges_pvaluesum <- data.frame("Node1" = organismicGOfiltered$Node1,
                              "Node2" = organismicGOfiltered$Node2,
                              "p.value" = rep(NA, times=length(organismicGOfiltered$Node1)),
                              "num.p" = rep(0, times=length(organismicGOfiltered$Node1)))


pb <- txtProgressBar(min = 0, max = length(edges_pvaluesum$p.value), style = 3)
for(i in 1:length(edges_pvaluesum$p.value)) {
  x <- na.omit(as.numeric(ks_pval_gold_ecdf[i,3:length(colnames(ks_pval_gold_ecdf))]))
  edges_pvaluesum$num.p[i] <- length(x)
  edges_pvaluesum$p.value[i] <- fishersmethod(x)
  setTxtProgressBar(pb, i)
}
close(pb)


#new data frame of ordered edges, first gold standard, second putative
edges_order_pvaluesum <- edges_pvaluesum
edges_order_pvaluesum$Node1  <- as.character(edges_order_pvaluesum$Node1)
edges_order_pvaluesum$Node2  <- as.character(edges_order_pvaluesum$Node2)

pb <- txtProgressBar(min = 0, max = length(edges_order_pvaluesum$p.value), style = 3)
for(i in 1:length(edges_pvaluesum$p.value)) {
  x <- as.character(organismicGOfiltered$Node1[i])
  y <- as.character(organismicGOfiltered$Node2[i])
  if(x%in%gold_geneID_ec) {
    edges_order_pvaluesum$Node1[i] <- x
    edges_order_pvaluesum$Node2[i] <- y
  } else {
    edges_order_pvaluesum$Node2[i] <- x
    edges_order_pvaluesum$Node1[i] <- y
  }
  setTxtProgressBar(pb, i)
}
close(pb)

write.table(edges_order_pvaluesum, "results/edgelist_weighted.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


################################################################################
##                     subset p.values of gold standards                      ##
################################################################################

#how many edges have a small p.value #153255 of 478003
sum(edges_pvaluesum$p.value < 0.05)
sum(edges_pvaluesum$p.value < 0.2197)
sum(edges_pvaluesum$p.value < 0.1006) #1st Qu. of gold standards

x <- gold$Node1 == 1 & gold$Node2 == 1

edges_pvaluesum_gold <- edges_pvaluesum[x,]
edges_putative_putative <- edges_pvaluesum[!x,]

#statistics of gold standard edges
#mean = 0.2197
min(edges_pvaluesum_gold$p.value)
sum(edges_pvaluesum_gold$p.value == 0)
summary(edges_pvaluesum_gold$p.value)


################################################################################
##                               image similarity matrix                      ##
################################################################################

tiff("results/ecdf_exprs.tiff", width=4, height=4, units="in", res = 300)
plot(gecdf_expr, ylab = "Cumulative probability", xlab = "Coexpression similarity", main = "Coexpression ECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
tiff("results/ecdf_gobp.tiff", width=4, height=4, units="in", res = 300)
plot(gecdf_go_bp, ylab = "Cumulative probability", xlab = "GO biological process similarity", main = "GO biological process\nECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
tiff("results/ecdf_gomf.tiff", width=4, height=4, units="in", res = 300)
plot(gecdf_go_mf, ylab = "Cumulative probability", xlab = "GO molecular function similarity", main = "GO molecular function\nECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
tiff("results/ecdf_tf.tiff", width=4, height=4, units="in", res = 300)
plot(gecdf_tf_och, ylab = "Cumulative probability", xlab = "TF regulation similarity", main = "Transcription factors\nregulation ECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
tiff("results/ecdf_pfam.tiff", width=4, height=4, units="in", res = 300)
plot(gecdf_pfam_och, ylab = "Cumulative probability", xlab = "Pfam similarity", main = "Pfam ECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
tiff("results/ecdf_ppi.tiff", width=4, height=4, units="in", res = 300)
plot(gecdf_ppi_och, ylab = "Cumulative probability", xlab = "PPI similarity", main = "Protein-protein\ninteraction ECDF", cex.main = 2, cex.lab = 1.3)
dev.off()

pdf("results/ecdf_exprs.pdf", width=4, height=4)
plot(gecdf_expr, ylab = "Cumulative probability", xlab = "Coexpression similarity", main = "Coexpression ECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
pdf("results/ecdf_gobp.pdf", width=4, height=4)
plot(gecdf_go_bp, ylab = "Cumulative probability", xlab = "GO biological process similarity", main = "GO biological process\nECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
pdf("results/ecdf_gomf.pdf", width=4, height=4)
plot(gecdf_go_mf, ylab = "Cumulative probability", xlab = "GO molecular function similarity", main = "GO molecular function\nECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
pdf("results/ecdf_tf.pdf", width=4, height=4)
plot(gecdf_tf_och, ylab = "Cumulative probability", xlab = "TF regulation similarity", main = "Transcription factors\nregulation ECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
pdf("results/ecdf_pfam.pdf", width=4, height=4)
plot(gecdf_pfam_och, ylab = "Cumulative probability", xlab = "Pfam similarity", main = "Pfam ECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
pdf("results/ecdf_ppi.pdf", width=4, height=4)
plot(gecdf_ppi_och, ylab = "Cumulative probability", xlab = "PPI similarity", main = "Protein-protein\ninteraction ECDF", cex.main = 2, cex.lab = 1.3)
dev.off()
