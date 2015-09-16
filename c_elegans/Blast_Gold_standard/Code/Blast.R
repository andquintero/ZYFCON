################################################################################
##                                Dependencies                                ##
################################################################################

setwd("~/master_thesis/s_cerevisiae/Blast_Gold_standard/")
load(".RData")
save.image(".RData")

library(stringr)
################################################################################
##                         Load uniprot reference table                       ##
################################################################################
referencet <- readLines("Data/uniprot-ec%3A-+AND+reviewed%3Ayes.tab")
x <- strsplit(as.character(referencet), split="\t", perl=TRUE)
colnames <- x[[1]]
referencet <- as.data.frame(do.call(rbind, x[2:length(x)]))
colnames(referencet)  <- colnames
rm(x, colnames)

################################################################################
##                               extract EC number                            ##
################################################################################

x <- as.data.frame(str_match(referencet$"Protein names", "^.*EC ([0-9]+.[0-9]*.[0-9]*.[0-9]*).*$")[,-1])
referencet$EC  <- x[,1]

#delete general reaction descriptions
completeannotation <- referencet[sapply(strsplit(as.character(referencet$EC), "\\."), length) == 4,]

################################################################################
##                        create EC number GeneID dicc                        ##
################################################################################
#WARNING dictionary not working, need to use hash for faster processing
#create ec2geneid dictionary 
ecgene_dicc <- vector(mode="list", length=length(completeannotation$EC))
names(ecgene_dicc) <- completeannotation$"Entry name"
x <- as.character(completeannotation$EC)
for(i in 1:length(completeannotation$EC)) {
  ecgene_dicc[[i]] <- x[i]
}
rm(i, x)

################################################################################
##                             read BLAST result                              ##
################################################################################

blastedgenes <- read.table("Results/alluniprot.blasted")
x <- as.character(blastedgenes$V1)
x <- as.data.frame(str_match(x, "^(.*)\\|(.*)\\|(.*)$")[,-1])
blastedgenes$entryname <- x$V3
colnames(blastedgenes) <- c("query_id" , "subject_id" , "percent_identity" , "alignment_length" , "mismatches" , "gap_opens" , "q_start" , "q_end" , "s_start" , "s_end" , "evalue" , "bit_score", "entryname")
rm(x)

################################################################################
##                    filter genes with incomplete EC                         ##
################################################################################

x <- blastedgenes$entryname%in%completeannotation$"Entry name"
blastedgenes <- blastedgenes[x,]

x <- as.character(blastedgenes$subject_id)
x <- as.data.frame(str_match(x, "^(.*)\\|(.*)$")[,-1])
blasted_query_gene <- data.frame("geneID" = x$V1,
                                 "uniprot" = blastedgenes$entryname)

#Write blasted genes and query
write.table(blasted_query_gene, "Results/blasted_query_gene.txt", 
            col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


#run perl script to obtain gene -> ec file

#Write dictionary
dictionary <- data.frame("entry" = completeannotation$"Entry name",
                       "EC" = completeannotation$EC)
write.table(dictionary, "Data/uniprot_EC_dictionary.txt",
            col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


#wite blasted genes and query plus GoldStandards of GO
gogold <- read.table("../Gold_Standard/data/yeast_GS_EC-gene_GO.txt")
ec_gene <- read.table("Results/Blasted_genes.list")

ec_gene_plusgold <- merge(gogold,ec_gene, all = TRUE, sort = FALSE)
ec_gene_plusgold <- ec_gene_plusgold[!duplicated(ec_gene_plusgold),]
write.table(ec_gene_plusgold, "Results/Blasted_genes_plusGS.list", 
            col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

#run script for make network

################################################################################
##                        search GeneID in ec2geneid dicc                     ##
################################################################################
#WARNING this step takes to long, better a hash in perl
#search geneID in diccionary to make a putative EC number annotation
ec <- c()
# for(i in 1:length(blastedgenes$entryname)) {
#   ifelse(is.null(ecgene_dicc[[as.character(blastedgenes$entryname[i])]]), 
#          ec[i] <- "NA", 
#          ec[i] <- ecgene_dicc[[as.character(blastedgenes$entryname[i])]])
# }
x <- as.character(blastedgenes$entryname)
for(i in 1:length(blastedgenes$entryname)) {
  ec[i] <- ecgene_dicc[[x[i]]]
} #not working because blast was against all gene with ec, including incomplete ecs4

#EC number imputation to transcripts IDs
ec_transcript <- data.frame("EC" = ec,
                            "transcriptID" = blastedgenes$subject_id)
ec_transcript <- ec_transcript[!duplicated(ec_transcript),]


