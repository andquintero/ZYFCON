################################################################################
##                                Dependencies                                ##
################################################################################

setwd("~/master_thesis/s_cerevisiae/expression/")
load(".RData")
save.image(".RData")

source("http://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("Heatplus")

library(stringr)

library(affy)
library(doParallel)
library(minet)
library(bigmemory)
library(infotheo)


################################################################################
##                             read CELL files                                ##
################################################################################

cellfiles <- list.files(recursive = TRUE, pattern = "\\.[cC][eE][lL]\\.gz$|\\.[cC][eE][lL]$")
#273 cell files
myAB <- ReadAffy(filenames = cellfiles)
#rma normalization
eset <- rma(myAB) 

#extrac expression data
exp_mat <- exprs(eset) 
##Revisa dimensiones y los primeros probes almacenados en E:
dim(exp_mat)

################################################################################
##                     Function to format matrix                              ##
################################################################################

#read table of correspondency between probeid and geneID
probe=read.table("./data/probe.txt",sep="",header=TRUE, colClasses = "character")

#Funcion para mapear los probes a nombres de genes
##Argumentos de entrada:
##- a: Matriz de expresion (probes)
##- plataforma: Nombre de la plataforma de microarreglos ("Affymetrix", "BGIYale", "Agilent44K" y "Agilent-012106")
##- probe: archivo que guarda la correspondencia entre probes y genes para cada plataforma 
##Objetos de salida:
##- c: Matriz de expresion (genes)

do.map=function(a,plataforma,probe){
  #Filtra por plataforma y deja solo la columna MSU6 (nombres de fila son probeid)
  probeB=probe[probe$platform==plataforma,names(probe)%in%c("probeid","ORF")]#head(probeB)
  #Convierte la matriz en un data.frame
  probeC=as.data.frame(probeB,row.names=probeB$probeid)#head(probeC) 
  #Combinacion de data.frames
  b=new.cbind(probeC,as.data.frame(a)) #head(b)    
  #Elimina columnas donde estaban los nombres de probes y de genes para acelerar la conversi?n a matriz
  c=b[,-c(1,2)]#head(c)
  #Convierte de data.frame a matrix para poder agregar nombres de fila repetidos
  c=as.matrix(c)#head(c)
  #Convierte los datos de las columnas a valores numericos
  c=apply(c,2,as.numeric)#head(c)
  #Agrega nombres de fila
  rownames(c)<-b[,2]#head(c)  
  #Para genes que aparecen en m?s de un probe, toma el valor m?ximo entre los probes para cada ensayo
  c=aggregate(c, list(rownames(c)), max)#head(c)
  #Agrega nombres de fila
  rownames(c)<-c[,1]#head(c)
  #Elimina columna donde estaban los nombres de fila
  c=c[,-1]#head(c)  
  return(c)
}

#Funci?n para combinar tablas por nombres de filas
new.cbind <- function(...){
  input <- eval(substitute(list(...), env = parent.frame()))  
  names.orig <- NULL
  nrows <- numeric()
  for (i in 1:length(input))
  {
    nrows[i] <- nrow(input[[i]])
    names.orig <- c(names.orig, colnames(input[[i]])) 
  }  
  idx <- (1:length(input))[order(nrows, decreasing=T)]
  x <- NULL
  for (i in 1:length(input))
  {
    x <- c(x, rownames(input[[idx[i]]]))
  }  
  r <- data.frame(row.names=unique(x))
  for (i in 1:length(input))
  {
    r <- cbind(r, data.frame(input[[i]][match(rownames(r), rownames(input[[i]])),]))
  }  
  colnames(r) <- names.orig  
  return(r)
}

################################################################################
##                             format expression matrix                       ##
################################################################################

#read table of correspondency between probeid and geneID
probe=read.table("./data/probe.txt",sep="",header=TRUE, colClasses = "character")

#Replace probes to genes
E <- do.map(exp_mat,"Affymetrix",probe)#El n?mero indica la longitud del nombre del gen

#check E
dim(E)
head(E)
#remove N/A row
E <- E[-1,]

#save expression matrix
write.table(E,"./results/E.txt",sep="\t",col.names=TRUE)


################################################################################
##                          mutual information coefficient                    ##
################################################################################

#Crea matrices de bigmemory
dbs <- E
n=nrow(dbs)
#Crea bigmatrix x (similitudpares) para guardar la similitud por pares
x=filebacked.big.matrix(n*(n-1)/2, 1, type="double", backingfile="similitudpares.bin", descriptorfile="similitudpares.desc") #Por MI  
#Crea bigmatrix x (similitudmatriz) para guardar la similitud en matriz
mm=filebacked.big.matrix(n, n, type="double", backingfile="similitudmatriz.bin", descriptorfile="similitudmatriz.desc")
#Crea bigmatrix z (pares) para guardar las permutaciones
z=filebacked.big.matrix(n*(n-1)/2, 2, type="integer", backingfile="pares.bin", descriptorfile="pares.desc")
z[,1]<-unlist(mclapply(1:n,function(i) rep(i,n-i)))
z[,2]<-unlist(mclapply(2:n,function(i) seq(i,n)))#z[1:10,];is.big.matrix(z);is.filebacked(z);is.separated(z);nrow(z);ncol(z);nrow(x);ncol(x)


#Selecciona procesadores (? crea clusters) para trabajar en paralelo
Sys.getenv(5)
cl <- makeCluster(5)
registerDoParallel(cl)
#stopCluster(cl)#Detiene el cluster. 

#Discretiza base de datos en bins seg?n:
dbd=discretize(t(dbs), disc="globalequalwidth",nbins=floor(1+log2(dim(dbs)[2]))) #Notas: En soybean da un poco sesgada a la izquierda, en tomato con nbins=4 aparecen muchas estimaciones negativas

#Argumentos para cargar el archivo descriptor y luego las bigmatrix
desc<-dget("similitudpares.desc")
x=attach.big.matrix(desc);flush(x)#Para actualizar el filebacked #write.big.matrix(z, "per.txt", row.names = FALSE, col.names = FALSE, sep=",")#describe(z)

#Funcion de informacion mutua modificada
mimodificada=function(X,Y,ensy){
  U <- cbind(Y, X)
  Hx <- .Call("entropyR",X, ensy, 1, 2, DUP = FALSE, PACKAGE = "infotheo")
  Hy <- .Call("entropyR",Y, ensy, 1, 2, DUP = FALSE, PACKAGE = "infotheo")
  Hyx <- .Call("entropyR",U, ensy, 2, 2, DUP = FALSE, PACKAGE = "infotheo")
  res <- Hx + Hy - Hyx
  res
}

#Calcula similitud para cada par de la permutacion
(n=dim(dbd)[2])#N?mero de genes
(ens=dim(dbd)[1])#Numero de ensayos; sqrt(ens);i=1
k=1#Gen de inicio
paso=100#Cada cuantos genes 
while(k<n){
  print(k)
  if((k+paso)>n){paso=1}
  lista<-foreach(i=k:(k+paso), .export=c("mutinformation")) %dopar% {
    unlist(lapply((i+1):n, function(j) {mimodificada(dbd[,i],dbd[,j],ens)}))
  }
  IM=unlist(lista)#
  IM[which(IM<0)]<-0#
  end=sum(seq(n-1,n-(k+paso)))
  start=sum(seq(n-1,n-k))-(n-(k+1))
  x[start:end,1]<-sqrt(1-exp(-2*IM))
  flush(x)
  remove(lista)
  gc()
  k=k+paso
  #print(k)
  #print(end)
  write(k,"conteo.txt")
}

#Crea bigmatrix S (similitudmatriz) para guardar la similitud en matriz
desc<-dget("similitudmatriz.desc")
S=attach.big.matrix(desc)
dim(S);flush(S)

#Pasa los pares de similitud a la matriz completa de similitud
start=1#Posicion de inicio en la matriz de permutacion #k=5321#Ultima fila #combin=mwhich(z,cols=1,vals=list(c(k)),comps=list(c('eq')))
for(i in 1:(n-1)){
  end=start+(n-i-1)#Posicion final
  S[(i:n),i]<-c(1,x[(start:end),1])#Agrega correlacion 1 a la posicion i,i y pega el resto del vector
  start=end+1#Actualiza posicion inicial
  print(i)#Control de progreso
  flush.console()  
}
S[n,n]=1#Agrega correlacion 1 a la ultima celda mm[n-3,n-15]

dim(S)
S[1:10,1:10]

write.big.matrix(S, "results/E.csv")
write.table(rownames(E), "results/E_names.txt", quote = F, row.names = F, col.names = F)
################################################################################
##                           read big similarity matrix                       ##
################################################################################

expression_sim <- read.csv("results/E.csv")
expression_sim <- rbind(rep(0, times=length(colnames(expression_sim))), expression_sim)
expression_sim <- expression_sim + t(expression_sim)
diag(expression_sim) <- 1

sim_geneID <- readLines("results/E_names.txt")

rownames(expression_sim)  <- sim_geneID
colnames(expression_sim)  <- sim_geneID

################################################################################
##                           organismic net read node list                    ##
################################################################################

#read organismic network and extrac gene node names
organismicGOfiltered <- read.table("../Blast_Gold_standard/Results/Gen-Gen.edgelist_GS-GO_filtered.txt", col.names=c("Node1", "Node2"))

x <- as.data.frame(str_match(organismicGOfiltered$Node1, "^(.*);(.*$)")[,-1])
y <- as.data.frame(str_match(organismicGOfiltered$Node2, "^(.*);(.*$)")[,-1])

node1 <- x
node2 <- y
colnames(node1) <- c("GeneID", "EC")
colnames(node2) <- c("GeneID", "EC")

rm(x ,y)

################################################################################
##              assign weight to edges in organismicGOfiltered                ##
################################################################################

x <- as.character(node1$GeneID)
y <- as.character(node2$GeneID)
z <- as.character(rownames(expression_sim))

organismicGOfiltered$expr <- rep(0, times=length(organismicGOfiltered$Node1))
pb <- txtProgressBar(min = 0, max = length(organismicGOfiltered$expr), style = 3)
for(i in 1:length(organismicGOfiltered$expr)) {
  ifelse(x[i]%in%z & y[i]%in%z, 
         organismicGOfiltered$expr[i] <- expression_sim[x[i], y[i]],
         organismicGOfiltered$expr[i] <- NA)
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

sim <- as.matrix(expression_sim[gogold$V2,gogold$V2])
dim(sim)
sim[is.na(sim)] <- 0

pdf("results/gold_sim_image_exprs.pdf", width=4, height=4)
heatmap_2(sim, Rowv=F, Colv=F, do.dendro=c(F, F), col=grey(seq(0, 1, length = 256)), scale="none", legfrac = 0 )
dev.off()

tiff("results/gold_sim_image_exprs.tiff", width=4, height=4, units="in", res = 300)
heatmap_2(sim, Rowv=F, Colv=F, do.dendro=c(F, F), col=grey(seq(0, 1, length = 256)), scale="none", legfrac = 0 )
dev.off()
