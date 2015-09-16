################################################################################
##                                Dependencies                                ##
################################################################################

source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
library(stringr)

setwd("~/master_thesis/Arabidopsis/genetic_interaction/")
load(".RData")

################################################################################
##                             Read BioGRID table                             ##
################################################################################

ara_biogrid <- read.table("Data/BIOGRID-ORGANISM-3.3.124.tab2/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-3.3.124.tab2.txt", sep="\t")

ara_biogrid <- read.table("Data/BIOGRID-ORGANISM-3.3.124.mitab/BIOGRID-ORGANISM-Arabidopsis_thaliana_Columbia-3.3.124.mitab.txt", sep="\t", fill=TRUE)

