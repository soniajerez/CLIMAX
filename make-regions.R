
regionalizacion.f <- function(ejecucion,modeS,matrizin,matrizinaux,comp,etapacorte,ngrupos,coordenadas,nest) {

library("MASS")
library("cluster")
library("mclust")

## EXECUTION #1 ######################################################################################################

if(ejecucion == 1) {

# READ INPUT DATA
matrizdatos <<- read.table(matrizin)
matrizdatosaux  <<- read.table(matrizinaux)
pcp <<- matrizdatos
nest <<- dim(pcp)[2]
etapas <<- nest-1
estaciones <<- as.character(1:nest)

print('Datos leidos')

if(modeS == "TRUE"){

# CORRELATION MATRIX BASED PCA

Correlacion <<- cor(pcp,pcp)
print('Matriz de correlaciones calculada')

comp.pcp.cov <<- eigen(Correlacion)
print('EOFs calculadas')

pcs <<- pcp*comp.pcp.cov$vectors
print('PCs obtenidas')

write.table(pcs,"PC.txt")
write.table(comp.pcp.cov$vectors,"EOF.txt")

varianza <<- array(dim=nest)
for(i in 1:nest){
varianza[i] <<- comp.pcp.cov$values[i]/sum(comp.pcp.cov$values)
}

sumavarianza <<- array(dim=nest)
for(i in 1:nest) {
 if(i==1){
  sumavarianza[i] <<- varianza[i]
  }
 if(i!=1){
 j <- i-1
 sumavarianza[i] <<- varianza[i] + sumavarianza[j]
 }
}
write.table(sumavarianza,"varianzaexplicadaacumulada.txt")
write.table(varianza,"varianzaexplicada.txt")

}

} 

## EXECUTION #2 #######################################################################################################

if(ejecucion==2){

if(modeS == "TRUE"){
retencion <<- as.matrix(comp.pcp.cov$vectors[,1:comp])
distcor.comp <<- dist(retencion,"euclidean")

# WARD
wardcor.comp <<- hclust(distcor.comp,"ward")
distancias.ward <<- wardcor.comp$height

netapas <- dim(pcp)[2]-1
inicio <- netapas-40

write.table(distancias.ward,"distancias.txt")

}

} 

## EXECUTION #3 ##########################################################################################################

if(ejecucion==3){

lonlat <<- read.table(coordenadas) 

# WARD SEEDS
corte <<- distancias.ward[etapacorte]
Semillasward <<- matrix(nrow=ngrupos,ncol=comp)


pdf("dendograma.pdf")
plot(wardcor.comp)
rect.hclust(wardcor.comp,h=corte)

GRupos <<- rect.hclust(wardcor.comp,h=corte,border=2:3)

for(i in 1:ngrupos){
 for(j in 1:comp) {
   Semillasward[i,j] <- mean(retencion[GRupos[[i]],j])
  }
 }

# K-MEANS
t.kmeans <<- kmeans(retencion,Semillasward)
retencion <<- as.matrix(comp.pcp.cov$vectors[,1:comp])
grupos.kmeans <- cbind(1:nest,as.data.frame(t.kmeans$cluster))

tamano <<- array(dim=ngrupos)
for(i in 1:ngrupos) {
tamano[i] <- dim(subset(grupos.kmeans,grupos.kmeans[,2]==i))[1]
}

ss <- sort(tamano,TRUE)
orden <- array(dim=length(ss))

tt <- cbind(tamano,1:ngrupos)

for(i in 1:ngrupos){
 cont <- 0
 for(j in 1:dim(tt)[1]){
  if(cont == 0){
  if(ss[i]==tt[j,1]){
  orden[i] <- tt[j,2]
  tt[j,] <- c(0,0)
  cont <- cont+1}
  }
 }
}

ordenado <- array(dim=(dim(grupos.kmeans)[1]))
gruposorden <- grupos.kmeans
for(i in 1:length(ss)){
  for(j in 1:(dim(gruposorden)[1])){
  if(gruposorden[j,2]==orden[i]){
    ordenado[j] <- i }
 } }

grupos.kmeans.ordenado <- cbind(gruposorden[,1],ordenado)

agrupacion.final <<- cbind(lonlat,ordenado)

write.table(agrupacion.final,"regiones.txt",row.names=FALSE,col.names=FALSE)

REgionalizacion <<- vector("list",2)

REgionalizacion[[1]] <<- agrupacion.final 
REgionalizacion[[2]] <<- sort(tamano,decreasing=TRUE) 

names(REgionalizacion) <<- c("agrupacion","tamanogrupos")

agrupacion <- REgionalizacion[[1]]
ngrupos <- length(REgionalizacion[[2]])

valorMedioRegional <- matrix(nrow=dim(matrizdatos)[1],ncol=ngrupos)
for(i in 1:ngrupos){
  elementos <- which(agrupacion[,3]==i)
  valorMedioRegional[,i] <- rowMeans(matrizdatos[,elementos])
  j <- as.character(i)
  write.table(valorMedioRegional[,i],paste("seriemedia-region",j,".txt",sep=""))
      
}

REgionalizacion <<- vector("list",2)

REgionalizacion[[1]] <<- agrupacion.final 
REgionalizacion[[2]] <<- sort(tamano,decreasing=TRUE)  

names(REgionalizacion) <<- c("agrupacion","tamanogrupos")

agrupacion <- REgionalizacion[[1]]
ngrupos <- length(REgionalizacion[[2]])

nn=ngrupos+1
valorMedioRegional <- matrix(nrow=dim(matrizdatosaux)[1],ncol=nn)
for(i in 1:ngrupos){
  elementos <- which(agrupacion[,3]==i)
  valorMedioRegional[,i] <- rowMeans(matrizdatosaux[,elementos])
  j <- as.character(i)
  write.table(valorMedioRegional[,i],paste("seriemedia-aux-region",j,".txt",sep=""))     
}
  valorMedioRegional[,nn] <- rowMeans(matrizdatosaux[,])
  write.table(valorMedioRegional[,nn],paste("seriemedia-aux-allregions.txt"))


}

}




