################################################################################
##                                Dependencies                                ##
################################################################################

source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
library(stringr)

setwd("~/master_thesis/Arabidopsis/Gold_Standard/")
load(".RData")

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
##                          Listing and choosing filters                      ##
################################################################################

filters <- listFilters(ensemblPlants)
attributes <- listAttributes(ensemblPlants)
filterOptions("biotype", ensemblPlants)
?getBM
ara_go_ec <- getBM(filters = "biotype",
                   values = "protein_coding",
                   attributes = c("ensembl_gene_id",
                                  "ensembl_transcript_id",
                                  "go_accession","go_linkage_type", "go_namespace_1003",
                                  "kegg_enzyme"),
                   mart=ensemblPlants)

write.csv(ara_go_ec, file="Data/ara_go_ec.csv")

################################################################################
##                                     ec2go                                  ##
################################################################################

ec2go <- readLines("Data/ec2go.txt")
ec2go <- as.data.frame(str_match(ec2go[3:length(ec2go)], "^EC:(.*) > (.*); (.*)$")[,-1])
colnames(ec2go) <- c("EC", "description", "GO")

#delete general reaction descriptions
ec2go <- ec2go[sapply(strsplit(as.character(ec2go$EC), "\\."), length) == 4,]

#delete general reaction
x <- strsplit(as.character(ec2go$EC), split="\\.", perl=TRUE)
x <- as.data.frame(do.call(rbind, x))
ec2go <- ec2go[x$V4 != "-",]

################################################################################
##                  Extract only genes with validated EC numbers              ##
################################################################################

#only genes with GO annotations in ec2go mapping
ara_onlyEC <- ara_go_ec[ara_go_ec$go_accession%in%ec2go$GO,]

#only genes with validated EC numbers
x <- ara_onlyEC$go_linkage_type%in%c("EXP", "IDA", "IPI", "IMP", "IGI")
x <- ara_onlyEC$go_linkage_type == "EXP" | ara_onlyEC$go_linkage_type == "IDA"
ara_onlyEC <- ara_onlyEC[x,]

#create ec2go diccionary 
ecgo_dicc <- vector(mode="list", length=length(ec2go$GO))
names(ecgo_dicc) <- ec2go$GO
x <- as.character(ec2go$EC)
for(i in 1:length(ec2go$GO)) {
  ecgo_dicc[[i]] <- x[i]
}

#search GO in diccionary for exp validated genes
ec <- c()
for(i in 1:length(ara_onlyEC$go_accession)) {
  ec[i] <- ecgo_dicc[[as.character(ara_onlyEC$go_accession[i])]]
}

#if gene id will be used:
exp_valEC <- data.frame("EC" = ec,
                        "gene_id" = ara_onlyEC$ensembl_gene_id)
exp_valEC <- exp_valEC[!duplicated(exp_valEC),]

#if transcript id will be used:
exp_valEC <- data.frame("transcript_id" = ara_onlyEC$ensembl_transcript_id,
                        "EC" = ec)

write.table(exp_valEC, file="Data/arab_GS_EC-gene_GO.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

################################################################################
##                           Compare GS and organismic net                    ##
################################################################################

#after using script 2 of metabolic network, for creation of the net
#read the GS net and the arabidopsis organismic net filtered by GO
#to see if all nodes and edges of the GSnet are in the organismic net

gen_gen_list_onlyGO <- read.table("results/Gen-Gen.list_onlyexpGO.txt", col.names=c("Node1", "Node2"))
organismicGOfiltered <- read.table("../Blast_Gold_standard/Results/Gen-Gen.edgelist_GS-GO_filtered.txt", col.names=c("Node1", "Node2"))
goldNodes <- unique(c(as.character(gen_gen_list_onlyGO$Node1), as.character(gen_gen_list_onlyGO$Node2)))
imputedNodes <- unique(c(as.character(organismicGOfiltered$Node1), as.character(organismicGOfiltered$Node2)))

sum(goldNodes%in%imputedNodes)


#cheack agains blasted genes
blasted_genes <- read.table("../Blast_Gold_standard/Results/Blasted_genes.list", col.names=c("EC", "GeneID"))
x <- as.data.frame(str_match(blasted_genes$GeneID, "^(.*)\\..$")[,-1])
blastedimputedNodes <- unique(as.character(x[,1]))
sum(as.character(exp_valEC$gene_id)%in%blastedimputedNodes)

#check against arabidopsis fasta
headers <- readLines("Data/databaseHeaders.txt")
headers <- as.data.frame(str_match(headers, "^>(.{9})\\..*$")[,-1])
headers <- unique(as.character(headers[,1]))

sum(as.character(exp_valEC$gene_id)%in%headers)


################################################################################
##                                Using KEGG                                  ##
################################################################################
ara_ec <- getBM(filters = "biotype",
                   values = "protein_coding",
                   attributes = c("ensembl_gene_id",
                                  "ensembl_transcript_id",
                                  "kegg_enzyme"),
                   mart=ensemblPlants)

ara_ec <- ara_ec[sapply(strsplit(as.character(ara_ec$kegg_enzyme), "\\+"), length) == 2,]
ara_ec$kegg_enzyme <- do.call(rbind, strsplit(as.character(ara_ec$kegg_enzyme), "\\+"))[,2]

#delete general reaction
x <- strsplit(as.character(ara_ec$kegg_enzyme), split="\\.", perl=TRUE)
x <- as.data.frame(do.call(rbind, x))
ara_ec <- ara_ec[x$V4 != "-",]

#remove dups
ara_ec <- ara_ec[!duplicated(ara_ec),]

#if gene id will be used:
kegg_ec <- data.frame("EC" = ara_ec$kegg_enzyme,
                        "gene_id" = ara_ec$ensembl_gene_id)
kegg_ec <- kegg_ec[!duplicated(kegg_ec),]

#if transcript id will be used:
exp_valEC_KEGG <- data.frame("transcript_id" = ara_onlyEC$ensembl_transcript_id,
                        "EC" = ec)

write.table(kegg_ec, file="Data/arab_GS_EC-gene_KEGG.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)



