################################################################################
##                                Dependencies                                ##
################################################################################

setwd("~/master_thesis/Arabidopsis/modules_each/")
load(".RData")
save.image(".RData")

install.packages("igraph")

library(stringr)
library(igraph)


################################################################################
##                           read Ecs and maps of KEGG                        ##
################################################################################

#read table of link between pathway and enzyme
ec_map <- read.table("http://rest.kegg.jp/link/pathway/enzyme", col.names = (c("ec", "map")))
ec_map$map <- substring(ec_map$map, 6)
ec_map$ec <- substring(ec_map$ec, 4)

#keep only KEGG maps
x <- substring(ec_map$map, 1, 3) == "map"
ec_map <- ec_map[x,]

#remove incomplete reactions
x <- strsplit(as.character(ec_map$ec), split="\\.", perl=TRUE)
x <- as.data.frame(do.call(rbind, x))
ec_map <- ec_map[x$V4 != "-",];x <- x[x$V4 != "-",]
ec_map <- ec_map[x$V3 != "-",];x <- x[x$V3 != "-",]
ec_map <- ec_map[x$V2 != "-",];x <- x[x$V2 != "-",]

rm(x)

################################################################################
##                           transform to binary matrix                       ##
################################################################################

ctable_ec_map <- table(ec_map) #as table
ctable_ec_map <- as.data.frame.matrix(ctable_ec_map, stringsAsFactors=TRUE) #as data.frame
ctable_ec_map <- data.matrix(ctable_ec_map) #as matrix
ctable_ec_map[ctable_ec_map > 1] <- 1
dim(ctable_ec_map)

################################################################################
##                  find expected number of pathways                          ##
################################################################################

expected_paths <- readLines("http://rest.kegg.jp/list/pathway/ath")
expected_paths <- as.data.frame(str_match(expected_paths[3:length(expected_paths)], "^(.*)\t(.*)$")[,-1])
expected_paths <- sum(substring(expected_paths[,1], 9)%in%substring(colnames(ctable_ec_map), 4))

################################################################################
##                                 read all networks                          ##
################################################################################

#function to read and filter network
read.net <- function(path) {
  #path= path of network
  edges <- read.table(path, header=TRUE)
  if(length(colnames(edges))==3) {
    colnames(edges)  <- c("Node1", "Node2", "weight")
  } else{
    colnames(edges)  <- c("Node1", "Node2", "weight", "num.p") 
  }
  edges <- edges[!duplicated(edges),]
  
  #filter small p.values < 0.05
  #
  edges <- edges[edges$weight >= 0.05,]
  #filter NA
  edges <- edges[complete.cases(edges),]

  #filter nodes autopredicting
  x <- as.data.frame(str_match(edges$Node1, "^(.*);(.*$)")[,-1])
  y <- as.data.frame(str_match(edges$Node2, "^(.*);(.*$)")[,-1])
  x <- as.character(x[,1]) == as.character(y[,1])
  edges <- edges[!x,]
  
  #read gold standard genesIDs and its EC 
  gogold <- read.table("../Gold_Standard/data/arab_GS_EC-gene_GO.txt")
  gold_geneID_ec <- paste(gogold$V2,";",gogold$V1, sep="")
  
  #node attributes
  nodeID <- unique(c(as.character(edges$Node1), as.character(edges$Node2)))
  nodes  <- data.frame("name" = nodeID,
                       "gold" = nodeID%in%gold_geneID_ec)
  
  #make it a graph
  organismic_net <- graph.data.frame(edges, directed=FALSE, vertices=nodes)
  return(organismic_net)
}

# integrated <- read.net("../integration/results/edgelist_weighted.txt");integrated
# exprs <- read.net("../expression/results/edgelist_weighted.txt")
# gomf <- read.net("../GO_similarity_mf/results/edgelist_weighted.txt")
# gobp <- read.net("../GO_similarity_bp/results/edgelist_weighted.txt")
# tf <- read.net("../transcription_factors/results/edgelist_weighted.txt")
# pfam <- read.net("../pfam/results/edgelist_weighted.txt")
# ppi <- read.net("../protein_protein/results/edgelist_weighted.txt")

integrated <- read.net("../integration/results/edgelist_weighted.txt")
exprs <- read.net("../integration/results/edgelist_weighted_exprs.txt")
gomf <- read.net("../integration/results/edgelist_weighted_gomf.txt")
gobp <- read.net("../integration/results/edgelist_weighted_gobp.txt")
tf <- read.net("../integration/results/edgelist_weighted_tf.txt")
pfam <- read.net("../integration/results/edgelist_weighted_pfam.txt")
ppi <- read.net("../integration/results/edgelist_weighted_ppi.txt")

integrated
exprs
gomf
gobp
tf
pfam
ppi


################################################################################
##                               find modularity                              ##
################################################################################



################################################################################
##        FUNCTION  calculate score of community presence in KEGG maps        ##
################################################################################

#function that looks all ec in the kegg maps 
#and compute the proportion of ec in a comunity for each map
proportion_KEGGmaps <- function(commu, maps) {
  #commu = communities data object
  #maps = binary matrix of EC number presence in KEGG maps
  commu_membership <- as.data.frame(membership(commu))
  commu_membership$ec <- as.data.frame(str_match(row.names(commu_membership), "^(.*);(.*$)")[,-1])[,2]
  commu_membership$geneID <- as.data.frame(str_match(row.names(commu_membership), "^(.*);(.*$)")[,-1])[,1]
#   commu_membership$ec <- substring(row.names(commu_membership), 11)
#   commu_membership$geneID <- substring(row.names(commu_membership), 1, 9)
  colnames(commu_membership)  <- c("membership", "ec", "geneID")
  
  map_presence <- matrix(nrow = length(unique(commu_membership$membership)), 
                         ncol = length(colnames(ctable_ec_map)), 
                         #paste("com",c(1:length(unique(commu_membership$membership))), sep="")
                         dimnames = list(c(1:length(unique(commu_membership$membership))), 
                                         colnames(ctable_ec_map)))
  
  # phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
  # q number of common ECs
  # m number of correct EC (map ECs)
  # n totalECs-mapECs
  # k number of EC in community
  pb <- txtProgressBar(min = 0, max = length(unique(commu_membership$membership)), style = 3)
  for(i in 1:length(unique(commu_membership$membership))) {
    #ecs_com <- commu_membership$ec[commu_membership$membership == i]
    ecs_com <- unique(commu_membership$ec[commu_membership$membership == i])
    for (j in 1:length(colnames(ctable_ec_map))) {
      ecs_map <- rownames(ctable_ec_map)[as.numeric(ctable_ec_map[,j]) == 1]
      tot_map <- sum(as.numeric(ctable_ec_map[,j]) == 1)
      #map_presence[i,j] <- phyper(sum(ecs_map%in%ecs_com), tot_map, 
      #                            length(unique(c(ecs_map, ecs_com)))-tot_map, length(ecs_com))
      map_presence[i,j] <- phyper(sum(ecs_map%in%ecs_com)-1, tot_map, 
                                  length(rownames(ctable_ec_map))-tot_map, length(ecs_com))
      
      #map_presence[i,j] <- sum(ecs_map%in%ecs_com)/tot_map
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(map_presence)
}

proportion_score <- function(data, expected) {
  #data= proportion_KEGGmaps output matrix
  #expected = an integer, number of expected pathways in orgamism
  x  <- c()
  for(i in 1:length(rownames(data))) {
    if(sum(data[i,]) > 0) {
      x[i] <- prod(data[i,data[i,] > 0])
    } else {
      x[i]  <- NA
    }
    #ifelse(sum(data[i,]) > 0 , x[i] <- prod(data[i,data[i,] > 0]),NA)
  }
  #y <- sum(na.omit(x))/(length(rownames(data))*abs(length(colnames(data))-length(rownames(data))+1))
  y <- sum(na.omit(x))/(length(rownames(data))*abs(expected-length(rownames(data))+1))
  return(list(row_product = x, proportion_score = y))
  #return(y)
}

npi <- function(graph, maps, expected) {
  #graph = a graph to randomize and calculate proportion score
  #maps = binary matrix of EC number presence in KEGG maps
  #expected = an integer, number of expected pathways in orgamism
  commu_net <- walktrap.community(graph, 
                                  weights = E(graph)$weight, 
                                  steps = 4)
  map_presence  <- proportion_KEGGmaps(commu = commu_net, 
                                       maps = maps)
  
  y <- proportion_score(map_presence, expected)
  y[["communities"]]  <- length(commu_net)
  y[["modularity_score"]]  <- modularity(commu_net)
  return(y)
} 

################################################################################
##              calculate score of community presence in KEGG maps            ##
################################################################################

npiintegrated <- npi(integrated, ctable_ec_map, expected_paths)
npiexprs <- npi(exprs, ctable_ec_map, expected_paths)
npigomf <- npi(gomf, ctable_ec_map, expected_paths)
npigobp <- npi(gobp, ctable_ec_map, expected_paths)
npitf <- npi(tf, ctable_ec_map, expected_paths)
npipfam <- npi(pfam, ctable_ec_map, expected_paths)
npippi <- npi(ppi, ctable_ec_map, expected_paths)

npiintegrated[-1]
npiexprs[-1]
npigomf[-1]
npigobp[-1]
npitf[-1]
npipfam[-1]
npippi[-1]

npiintegrated[[2]]
npiexprs[[2]]
npigomf[[2]]
npigobp[[2]]
npitf[[2]]
npipfam[[2]]
npippi[[2]]

corr_npi <- function(npires, uni){
  #npires= results of a NPI
  #uni = graph of only one data type
  #x <- npires[[2]]*sum(E(uni)$weight >0.8)
  #x <- npires[[2]]*sum(!(V(uni)$gold))
  x <- npires[[2]]*(sum((V(uni)$gold))/sum(!(V(uni)$gold)))
  if(!is.null(E(uni)$num.p)){
    x <- x*mean(E(uni)$num.p)
  }
  return(x)
}

corr_npi(npiintegrated, integrated)
corr_npi(npiexprs, exprs)
corr_npi(npigomf, gomf)
corr_npi(npigobp, gobp)
corr_npi(npitf, tf)
corr_npi(npipfam, pfam)
corr_npi(npippi, ppi)



