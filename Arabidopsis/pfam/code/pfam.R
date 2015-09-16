################################################################################
##                                Dependencies                                ##
################################################################################

source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
install.packages("ade4")

library(biomaRt)
library(stringr)
library(ade4)

setwd("~/master_thesis/Arabidopsis/pfam/")
load(".RData")
save.image(".RData")

################################################################################
##                           organismic net read node list                    ##
################################################################################

#read organismic network and extract gene node names
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
##                   Listing and choosing available databases                 ##
################################################################################

x <- listMarts();View(x)
ensemblPlants = useMart("plants_mart_26")
x <- listDatasets(ensembl);View(x)
ensemblPlants = useDataset("athaliana_eg_gene",mart=ensemblPlants)
#if one knows the database and the organism:
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

################################################################################
##                  Listing and choosing filters transcript ID                ##
################################################################################

filters <- listFilters(ensemblPlants)
attributes <- listAttributes(ensemblPlants)
filterOptions("biotype", ensemblPlants)
?getBM
ara_pfam <- getBM(filters = "biotype",
                   values = "protein_coding",
                   attributes = c("ensembl_gene_id",
                                  "ensembl_transcript_id",
                                  "pfam_pf"),
                   mart=ensemblPlants)

write.csv(ara_pfam, file="data/ara_pfam.csv")

################################################################################
##                  Listing and choosing filters gene ID                      ##
################################################################################

ara_pfam_geneID <- getBM(filters = "biotype",
                  values = "protein_coding",
                  attributes = c("ensembl_gene_id",
                                 "pfam_pf"),
                  mart=ensemblPlants)

#filter genes without fam annotation
x <- nchar(ara_pfam_geneID$pfam_pf) == 7
ara_pfam_geneID <- ara_pfam_geneID[x,]

################################################################################
##                           transform to binary matrix                       ##
################################################################################

ctable_pfam <- table(ara_pfam_geneID) #as table
ctable_pfam <- as.data.frame.matrix(ctable_pfam, stringsAsFactors=TRUE) #as data.frame
ctable_pfam <- data.matrix(ctable_pfam) #as matrix
ctable_pfam[ctable_pfam > 1] <- 1
dim(ctable_pfam)

#filter to keep only genes in organismic network
x <- rownames(ctable_pfam)%in%geneID
ctable_pfam <- ctable_pfam[x,]
dim(ctable_pfam)

################################################################################
##                        calculate binary similarities                       ##
################################################################################

dist.binary(ctable_pfam[1:10,], method = 1, diag = FALSE, upper = FALSE)
pfam_Jaccard <- dist.binary(ctable_pfam, method = 1, diag = FALSE, upper = FALSE)
pfam_Jaccard <- as.matrix(pfam_Jaccard)
pfam_Jaccard <- 1-pfam_Jaccard

pfam_Ochiai <- dist.binary(ctable_pfam, method = 7, diag = FALSE, upper = FALSE)
pfam_Ochiai <- as.matrix(pfam_Ochiai)
pfam_Ochiai <- 1-pfam_Ochiai



################################################################################
##              assign weight to edges in organismicGOfiltered                ##
################################################################################

x <- as.character(node1$GeneID)
y <- as.character(node2$GeneID)
z <- as.character(rownames(pfam_Jaccard))

organismicGOfiltered$pfam_jac <- rep(0, times=length(organismicGOfiltered$Node1))
pb <- txtProgressBar(min = 0, max = length(organismicGOfiltered$pfam_jac), style = 3)
for(i in 1:length(organismicGOfiltered$pfam_jac)) {
  ifelse(x[i]%in%z & y[i]%in%z, 
         organismicGOfiltered$pfam_jac[i] <- pfam_Jaccard[x[i], y[i]],
         organismicGOfiltered$pfam_jac[i] <- NA)
  setTxtProgressBar(pb, i)
}
close(pb)

organismicGOfiltered$pfam_och <- rep(0, times=length(organismicGOfiltered$Node1))
pb <- txtProgressBar(min = 0, max = length(organismicGOfiltered$pfam_och), style = 3)
for(i in 1:length(organismicGOfiltered$pfam_och)) {
  ifelse(x[i]%in%z & y[i]%in%z, 
         organismicGOfiltered$pfam_och[i] <- pfam_Ochiai[x[i], y[i]],
         organismicGOfiltered$pfam_och[i] <- NA)
  setTxtProgressBar(pb, i)
}
close(pb)

write.table(organismicGOfiltered, "results/edgelist_weighted.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

