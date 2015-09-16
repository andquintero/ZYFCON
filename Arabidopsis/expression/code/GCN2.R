#################################################
# Construcci?n GCN Arabidopsis
#################################################
install.packages("bigmemory")
install.packages("doParallel")
install.packages("bigmemory")

source("http://bioconductor.org/biocLite.R")
biocLite("minet")



#Selecciona directorio
path="~/master_thesis/Arabidopsis/GCN/"
setwd(path)
load("GCN.RData")

#Selecciona los experimentos y las condiciones
scontrol=read.csv("samplescontrol.csv",header=TRUE,colClasses = "character")

#Ubica samples en matriz de expresi?n global
path="C:/Users/BIOESTADISTICA/Desktop/GCN/3. Matrices"
setwd(path)
db=read.table("dbarabidopsis.txt",row.names=1,sep="\t",header=TRUE)
dbs=db[,names(db)%in%paste(scontrol$samples,".CEL",sep="")]
write.csv(names(db)[names(db)%in%paste(scontrol$samples,".CEL",sep="")],"samplescontrolfinal.csv",row.names = FALSE)

#######################
#Matrices de Similitud por Informaci?n Mutua
#######################

library(doParallel)
library(minet)
library(bigmemory)

#Crea matrices de bigmemory
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
Sys.getenv('7')
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


