#################################################
# Construcción GCN Arabidopsis
#################################################

#Selecciona directorio
path="C:/Users/BIOESTADISTICA/Desktop/GCN/"
setwd(path)

#Selecciona los experimentos y las condiciones
scontrol=read.csv("samplescontrol.csv",header=TRUE,colClasses = "character")

#Ubica samples en matriz de expresión global
path="C:/Users/BIOESTADISTICA/Desktop/GCN/3. Matrices"
setwd(path)
db=read.table("dbarabidopsis.txt",row.names=1,sep="\t",header=TRUE)
dbs=db[,names(db)%in%paste(scontrol$samples,".CEL",sep="")]
write.csv(names(db)[names(db)%in%paste(scontrol$samples,".CEL",sep="")],"samplescontrolfinal.csv",row.names = FALSE)

#######################
#Matrices de Similitud por Información Mutua
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


#Selecciona procesadores (ó crea clusters) para trabajar en paralelo
Sys.getenv('NUMBER_OF_PROCESSORS')
cl <- makeCluster(5)
registerDoParallel(cl)
#stopCluster(cl)#Detiene el cluster. 

#Discretiza base de datos en bins según:
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
(n=dim(dbd)[2])#Número de genes
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


#######################
#Umbrales de Similitud por Coeficientes de Agrupamiento
#######################

library(bigmemory);library(igraph);library(robfilter);library(minet);library(ade4)

#Crea bigmatrix para guardar los Ci a cada umbral
C=filebacked.big.matrix(n, 100, type="double", backingfile="ci.bin", descriptorfile="ci.desc")

#Crea bigmatrix para guardar los ki a cada umbral
K=filebacked.big.matrix(n, 100, type="double", backingfile="ki.bin", descriptorfile="ki.desc")

#Crea vector de taos
(ltaos=seq(0.99,0.01,by=-0.01))

#Calcula grado del nodo y CC en cada tao
for(tao in ltaos){
  print(tao)
  ##Construye o cargamatriz de adyacencia
  A=matrix(0,nrow=n,ncol=n)    
  ##Completa matriz de adyacencia
  for(i in 1:n){  
    v.unos=mwhich(S,cols=i,vals=tao,comps='gt')#Vector de unos para el gen i
    v.ceros=(1:n)[((1:n)%in%v.unos)==FALSE]#Vector de ceros para el gen i
    A[v.unos,i]<-1#Agrega unos 
    A[v.ceros,i]<-0#Agrega ceros
  }
  ##Calcula los Ci y ki
  Ag=graph.adjacency(A,mode="undirected")#Matriz de adyacencia
  remove(A);gc()
  Cv=transitivity(Ag,type="local")#Coeficiente de agrupamiento local
  Kv=degree(Ag,loops=FALSE)#Grado de los nodos
  desc<-dget("ki.desc");K=attach.big.matrix(desc)  
  desc<-dget("ci.desc");C=attach.big.matrix(desc)
  K[,round(tao*100,0)]<-Kv;C[,round(tao*100,0)]<-Cv;flush(K);flush(C)
  remove(K);remove(C)
}

##Calcula Cr y Co a los umbrales
desc<-dget("ki.desc");K=attach.big.matrix(desc)  
desc<-dget("ci.desc");C=attach.big.matrix(desc)
Cr=Co=rep(0,100);(ltaos=seq(0.01,0.99,by=0.01))
for(i in round(ltaos*100,0)){  
  gn<-which(K[,i]>=1)
  kn=length(gn)#número de genes con ki>=1
  k1=1/kn*sum(K[gn,i])#Variable para calcular el CC de la red aleatoria
  k2=1/kn*sum(K[gn,i]^2)#Variable para calcular el CC de la red aleatoria
  Co[i]=((k2-k1)^2)/(kn*k1^3)
  if(kn==0){Co[i]=0}
  gn<-which(K[,i]>1)
  kn=length(gn)#número de genes con ki>1
  Cr[i]=1/kn*sum(C[gn,i])#Coeficiente de clustering para la red real
  if(kn==0){Cr[i]=0}
}

#Encuentra el primer maximo local de la diferencia
(ltaos=seq(0.01,0.99,by=0.01))
dif=runmed(abs(Cr-Co),k=23,endrule="constant")[1:100]#dif=abs(Cr-Co)  
plot(ltaos,dif[ltaos*100],t="l",xlab="Umbral",ylab="|C-Co|")
plot(ltaos,abs(Cr-Co)[ltaos*100],t="l",xlab="Umbral",ylab="|C-Co|")
##Aplicando filtro de la mediana
umbral=identify(ltaos,abs(Cr-Co)[ltaos*100],n=1)/100

#######################
#Matrices de Adyacencia
#######################

remove(desc)
tao=umbral
desc<-dget("similitudmatriz.desc")
S=attach.big.matrix(desc)
write.csv(as.matrix(S),"S.csv",row.names = FALSE)
n=dim(S)[1]

#Filtra la matriz de similitud
genef=read.table("gene_names.txt",sep=">",header=FALSE)
genef=as.character(genef[,-1])
genef=substr(genef,start=1,stop=9)
genef2=genef[which(genef%in%rownames(dbs))]
Sf=as.matrix(S)[match(genef2,rownames(dbs)),match(genef2,rownames(dbs))]
write.csv(Sf,"Sfo.csv",row.names = FALSE)

gc()
##Completa matriz de adyacencia
A=filebacked.big.matrix(n, n, type="double", backingfile="ady.bin", descriptorfile="ady.desc")  

for(k in 1:n){  
  v.unos=mwhich(S,cols=k,vals=tao,comps='gt')#Vector de unos para el gen i
  v.ceros=(1:n)[((1:n)%in%v.unos)==FALSE]#Vector de ceros para el gen i
  A[v.unos,k]<-1#Agrega unos 
  A[v.ceros,k]<-0#Agrega ceros
  print(k)
}
flush(A)
write.csv(as.matrix(A),"A.csv",row.names = FALSE)
Af=as.matrix(A)[match(genef2,rownames(dbs)),match(genef2,rownames(dbs))]
write.csv(Af,"Afo.csv",row.names = FALSE)
write.csv(genef2,"gene_names2.csv",row.names = FALSE)
write.csv(rownames(dbs),"gene_names21122.csv",row.names = FALSE)
gc()

sum(Af[upper.tri(Af)])