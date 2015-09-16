################################################################################
##                                Dependencies                                ##
################################################################################

setwd("~/master_thesis/s_cerevisiae/modules/")
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

expected_paths <- readLines("http://rest.kegg.jp/list/pathway/sce")
expected_paths <- as.data.frame(str_match(expected_paths[3:length(expected_paths)], "^(.*)\t(.*)$")[,-1])
expected_paths <- sum(substring(expected_paths[,1], 9)%in%substring(colnames(ctable_ec_map), 4))



################################################################################
##                        organismic net read edge weights                    ##
################################################################################

#expression
edges <- read.table("../integration/results/edgelist_weighted.txt", header=TRUE)
colnames(edges)  <- c("Node1", "Node2", "weight", "num.p")
edges <- edges[!duplicated(edges),]

#filter small p.values < 0.05
edges <- edges[edges$weight >= 0.05,]

#filter nodes autopredicting
x <- as.data.frame(str_match(edges$Node1, "^(.*);(.*$)")[,-1])
y <- as.data.frame(str_match(edges$Node2, "^(.*);(.*$)")[,-1])

x <- as.character(x[,1]) == as.character(y[,1])
edges <- edges[!x,]


rm(x ,y)

################################################################################
##                           read gold standard geneID and ec                 ##
################################################################################

#read gold standard genesIDs and its EC 
gogold <- read.table("../Gold_Standard/data/yeast_GS_EC-gene_GO.txt")
gold_geneID_ec <- paste(gogold$V2,";",gogold$V1, sep="")

################################################################################
##                              create nnode atributes                        ##
################################################################################

nodeID <- unique(c(as.character(edges$Node1), as.character(edges$Node2)))

nodes  <- data.frame("name" = nodeID,
                     "gold" = nodeID%in%gold_geneID_ec)

################################################################################
##                                make it a graph                             ##
################################################################################

organismic_net <- graph.data.frame(edges, directed=FALSE, vertices=nodes)
organismic_net

################################################################################
##                               find modularity                              ##
################################################################################

commu_organismic_net <- walktrap.community(organismic_net, weights = E(organismic_net)$weight, steps = 4)
#commu_organismic_net <- walktrap.community(organismic_net, weights = E(organismic_net)$weight, steps = 5)

# organismic_net2 <- simplify(organismic_net)
# commu_organismic_net2 <- fastgreedy.community(organismic_net2, weights = E(organismic_net)$weight)


length(commu_organismic_net)
print(commu_organismic_net)
modularity(commu_organismic_net)


commu_membership <- as.data.frame(membership(commu_organismic_net))
commu_membership$ec <- as.data.frame(str_match(row.names(commu_membership), "^(.*);(.*$)")[,-1])[,2]
commu_membership$geneID <- as.data.frame(str_match(row.names(commu_membership), "^(.*);(.*$)")[,-1])[,1]
colnames(commu_membership)  <- c("membership", "ec", "geneID")
#plot(commu_organismic_net, organismic_net)


################################################################################
##        FUNCTION  calculate score of community presence in KEGG maps        ##
################################################################################
#probabilty of finding M ECs 
# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
# q number of common ECs
# m number of correct EC (map ECs)
# n totalECs-mapECs
# k number of EC in community

  

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

randomnet_proportion_score <- function(graph, reps, maps, expected) {
  #graph = a graph to randomize and calculate proportion score
  #reps = and integer to specify number of random networks to generate
  #maps = binary matrix of EC number presence in KEGG maps
  #expected = an integer, number of expected pathways in orgamism
  results=list() 
  pb <- txtProgressBar(min = 0, max = reps, style = 3)
  for(i in 1:reps) {
    random_net  <- rewire.edges(graph, prob = 0.7)
    commu_random_net <- walktrap.community(random_net, 
                                           weights = E(random_net)$weight, 
                                           steps = 4)
    randommap_presence  <- proportion_KEGGmaps(commu = commu_random_net, 
                                               maps = maps)
    x <- paste0("rand", i)
    y <- proportion_score(randommap_presence, expected)
    y[["communities"]]  <- length(commu_random_net)
    y[["modularity_score"]]  <- modularity(commu_random_net)
    results[[x]] <- y
    setTxtProgressBar(pb, i)
  }
  return(results)
  close(pb)
} 

length(commu_organismic_net)
print(commu_organismic_net)
modularity(commu_organismic_net)
################################################################################
##              calculate score of community presence in KEGG maps            ##
################################################################################

map_presence  <- proportion_KEGGmaps(commu = commu_organismic_net, maps = ctable_ec_map)
#map_presence[map_presence < 0.05] <- 0
x <- proportion_score(map_presence, expected_paths) #0.000172619
x[["communities"]]  <- length(commu_organismic_net)
x[["modularity_score"]]  <- modularity(commu_organismic_net)
x


################################################################################
##      RANDOM NETS  calculate score of community presence in KEGG maps       ##
################################################################################

random_2 <- randomnet_proportion_score(graph = organismic_net, reps = 2, maps = ctable_ec_map, expected = expected_paths)
random_2[[1]][[2]]
random_2[[2]][[2]]


random_100.1 <- randomnet_proportion_score(graph = organismic_net, reps = 100, maps = ctable_ec_map, expected = expected_paths)
random_100.2 <- randomnet_proportion_score(graph = organismic_net, reps = 100, maps = ctable_ec_map, expected = expected_paths)
random_100.3 <- randomnet_proportion_score(graph = organismic_net, reps = 100, maps = ctable_ec_map, expected = expected_paths)
random_100.4 <- randomnet_proportion_score(graph = organismic_net, reps = 100, maps = ctable_ec_map, expected = expected_paths)
random_100.5 <- randomnet_proportion_score(graph = organismic_net, reps = 100, maps = ctable_ec_map, expected = expected_paths)
random_100.6 <- randomnet_proportion_score(graph = organismic_net, reps = 100, maps = ctable_ec_map, expected = expected_paths)
random_100.7 <- randomnet_proportion_score(graph = organismic_net, reps = 100, maps = ctable_ec_map, expected = expected_paths)
random_100.8 <- randomnet_proportion_score(graph = organismic_net, reps = 100, maps = ctable_ec_map, expected = expected_paths)
random_100.9 <- randomnet_proportion_score(graph = organismic_net, reps = 100, maps = ctable_ec_map, expected = expected_paths)
random_100.10 <- randomnet_proportion_score(graph = organismic_net, reps = 100, maps = ctable_ec_map, expected = expected_paths)

names(random_100.2) <- paste0("rand", c(101:200))
names(random_100.3) <- paste0("rand", c(201:300))
names(random_100.4) <- paste0("rand", c(301:400))
names(random_100.5) <- paste0("rand", c(401:500))
names(random_100.6) <- paste0("rand", c(501:600))
names(random_100.7) <- paste0("rand", c(601:700))
names(random_100.8) <- paste0("rand", c(701:800))
names(random_100.9) <- paste0("rand", c(801:900))
names(random_100.10) <- paste0("rand", c(901:1000))

#make a unique list of all random networks
all_random <- c(random_100.1, random_100.2, random_100.3, random_100.4, random_100.5, random_100.6, random_100.7, random_100.8, random_100.9, random_100.10)
names(all_random)

y <- c()
for(i in 1:length(all_random)) {
  y[i] <- all_random[[i]][[4]]
}
y
summary(y)


y <- c()
for(i in 1:length(all_random)) {
  y[i] <- length(all_random[[i]][[1]])
}
y
summary(y)
################################################################################
##                           order edges by p.value                           ##
################################################################################

edges_pvalue_order <- edges[order(edges$weight, decreasing = TRUE),]
write.table(edges_pvalue_order, "results/edgelist_weighted.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#identify gold standards in organismic net
x <- edges_pvalue_order$Node1%in%gold_geneID_ec
y <- edges_pvalue_order$Node2%in%gold_geneID_ec

gold <- data.frame("Node1" = ifelse(x, 1, 0),
                   "Node2" = ifelse(y, 1, 0))
rm(x,y)

x <- gold$Node1 == 1 & gold$Node2 == 1
x <- gold$Node2 == 1

edges_gold_pvalue_order <- edges_pvalue_order[x,]
edges_putative_pvalue_order <- edges_pvalue_order[!x,]

write.table(edges_putative_pvalue_order, "results/edgelist_putative_weighted.txt", 
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


################################################################################
##                           annotated net for image                          ##
################################################################################
#create net
organismic_net_annot <- organismic_net

#read table of not combined probabilities
nocomb <- read.table("../integration/results/edgelist_weighted_all.txt", header=TRUE)
edges <- read.table("../integration/results/edgelist_weighted.txt", header=TRUE)
colnames(edges)  <- c("Node1", "Node2", "weight", "num.p")
x <- !duplicated(edges)
nocomb <- nocomb[x,]
edges <- edges[x,]

#filter small p.values < 0.05
x <- edges$weight >= 0.05
nocomb <- nocomb[x,]
edges <- edges[x,]

#filter nodes autopredicting
x <- as.data.frame(str_match(edges$Node1, "^(.*);(.*$)")[,-1])
y <- as.data.frame(str_match(edges$Node2, "^(.*);(.*$)")[,-1])
x <- as.character(x[,1]) == as.character(y[,1])
nocomb <- nocomb[!x,]
edges <- edges[!x,]
rm(x ,y)

organismic_net_annot

E(organismic_net_annot)$exprs <- !is.na(nocomb$expr)
E(organismic_net_annot)$gomf <- !is.na(nocomb$go_mf)
E(organismic_net_annot)$gobp <- !is.na(nocomb$go_bp)
E(organismic_net_annot)$tf <- !is.na(nocomb$tf_och)
E(organismic_net_annot)$pfam <- !is.na(nocomb$pfam_och)
E(organismic_net_annot)$ppi <- !is.na(nocomb$ppi_och)
E(organismic_net_annot)$ppi <- !is.na(nocomb$ppi_och)

E(organismic_net_annot)$p.exprs <- nocomb$expr
E(organismic_net_annot)$p.gomf <- nocomb$go_mf
E(organismic_net_annot)$p.gobp <- nocomb$go_bp
E(organismic_net_annot)$p.tf <- nocomb$tf_och
E(organismic_net_annot)$p.pfam <- nocomb$pfam_och
E(organismic_net_annot)$p.ppi <- nocomb$ppi_och
E(organismic_net_annot)$p.ppi <- nocomb$ppi_och


E(organismic_net_annot)$num.p <- paste0("sum_", edges$num.p)
V(organismic_net_annot)$commu <- paste0("com_", commu_membership$membership)

list.edge.attributes(organismic_net_annot)
list.vertex.attributes(organismic_net_annot)

write.graph(organismic_net_annot,"Results/all.graphml",format="graphml")

################################################################################
##                              assigned vs KO                                ##
################################################################################

KO <- read.table("http://rest.kegg.jp/link/enzyme/sce")
x <- as.data.frame(str_match(KO$V1, "^(.*):(.*$)")[,-1])
y <- as.data.frame(str_match(KO$V2, "^(.*):(.*$)")[,-1])
KO$compare <- paste0(x$V2,";", y$V2)
KO$ec  <- y$V2
KO$id  <- x$V2

y <- as.data.frame(str_match(edges_putative_pvalue_order$Node2, "^(.*);(.*$)")[,-1])
y$only3 <- as.data.frame(str_match(y$V2, "^([0-9]*.[0-9]*.[0-9]*).(.*$)")[,-1])[,1]
edges_putative_pvalue_order$compare3 <- paste0(y$V1,";", y$only3)
y$only2 <- as.data.frame(str_match(y$V2, "^([0-9]*.[0-9]*).(.*$)")[,-1])[,1]
edges_putative_pvalue_order$compare2 <- paste0(y$V1,";", y$only2)
y$only1 <- as.data.frame(str_match(y$V2, "^([0-9]*).(.*$)")[,-1])[,1]
edges_putative_pvalue_order$compare1 <- paste0(y$V1,";", y$only1)


KOcom <- KO[which(as.character(x$V2)%in%unique(as.character(y$V1))),]



KOcom$only3 <- as.data.frame(str_match(KOcom$ec, "^([0-9]*.[0-9]*.[0-9]*).(.*$)")[,-1])[,1]
KOcom$compare3 <- paste0(KOcom$id,";", KOcom$only3)
KOcom$only2 <- as.data.frame(str_match(KOcom$ec, "^([0-9]*.[0-9]*).(.*$)")[,-1])[,1]
KOcom$compare2 <- paste0(KOcom$id,";", KOcom$only2)
KOcom$only1 <- as.data.frame(str_match(KOcom$ec, "^([0-9]*).(.*$)")[,-1])[,1]
KOcom$compare1 <- paste0(KOcom$id,";", KOcom$only1)

x <- unique(as.character(KOcom$compare))%in%unique(as.character(edges_putative_pvalue_order$Node2))
sum(x)/nrow(KOcom) #[1] 0.5897674 complete 

x <- unique(as.character(KOcom$compare3))%in%unique(as.character(edges_putative_pvalue_order$compare3))
sum(x)/nrow(KOcom) #[1] 0.7646512 3

x <- unique(as.character(KOcom$compare2))%in%unique(as.character(edges_putative_pvalue_order$compare2))
sum(x)/nrow(KOcom) #[1] 0.7813953 2

x <- unique(as.character(KOcom$compare1))%in%unique(as.character(edges_putative_pvalue_order$compare1))
sum(x)/nrow(KOcom) #[1] 0.8195349 1
