rm(list = ls(all = TRUE))
library("Matrix")
library("mgcv") #smooth.con
library("dplyr")
library("kableExtra")
library("lattice")# 3D plot
library("rgl") #3D plot
library("scatterplot3d") #3d plot
source("funcoes/funcoes.R")
source("funcoes/funcoesMisturas.R")
#cores<-c("#BAF3DE","#FF9B95","#C9E69E","#FFC29A")
cores<-c("#7fcdbb","#ff5858")
par<-par(pch=19)

g<-2

# # # Cenario 1.1
pii<-c(0.35, 0.65)
nome.plot<-"c1.1_2_"

beta.verd<-list(c(4,4,6),c(2,3,-5))

sigma2.verd<-list(1,1)

arg.grupos<-list(g1=list(beta=beta.verd[[1]],curva=list(f="c1",a=-1,b=1,d=2),sigma2=sigma2.verd[[1]],intercepto=T),
                 g2=list(beta=beta.verd[[2]],curva=list(f="c2",a=-1,b=1,d=4),sigma2=sigma2.verd[[2]],intercepto=T),
                 intervalos=list(c(0,1),c(0,1)))

n=1000
alfas<-c(0.1,0.1)

amostra<-rMisUniPLM(n, pii, p=length(beta.verd[[1]]), arg=arg.grupos)
y=amostra$y
t=amostra$t
X=amostra$X
grupo<- apply(amostra$clu,2,FUN=function(x)which(as.logical(x))) # numero do grupo de cada observação

k=3

par(mfrow=c(1,3), mar=c(2,2,1,1))

s3d <- scatterplot3d(z=y,x=X[,3],y=X[,2], 
                     xlab = "x2",
                     ylab = "x1",
                     type = "p", 
                     color = cores[grupo],
                     angle = 55, scale.y = 0.7, 
                     pch = c(17,20)[grupo], cex.symbols = 1.3,
                     main = paste0("Gráfico de dispersão 3D - C",k) )
s3d$plane3d(beta.verd[[1]][1],x.coef = beta.verd[[1]][3], y.coef = beta.verd[[1]][2], col=cores[1], lwd=2)

s3d$plane3d(beta.verd[[2]][1],x.coef = beta.verd[[2]][3], y.coef = beta.verd[[2]][2], col=cores[2], lwd=2)

legend("topright", legend = c("Grupo 1","Grupo 2"), fill = cores[1:2])


s3d <- scatterplot3d(z=y,x=X[,3],y=t, 
                     xlab = "x2",
                     ylab = "t",
                     type = "p", 
                     color = cores[grupo],
                     angle = 55, scale.y = 0.7, 
                     pch = c(17,20)[grupo], cex.symbols = 1.3,
                     main = paste0("Gráfico de dispersão 3D - C",k) )
s3d$plane3d(mean(y),x.coef =0, y.coef = 0, lwd=2)

legend("topright", legend = c("Grupo 1","Grupo 2"), fill = cores[1:2])


s3d <- scatterplot3d(z=y,x=X[,2],y=t, 
                     xlab = "x1",
                     ylab = "t",
                     type = "p", 
                     color = cores[grupo],
                     angle = 55, scale.y = 0.7, 
                     pch = c(17,20)[grupo], cex.symbols = 1.3,
                     main = paste0("Gráfico de dispersão 3D - C",k) )
s3d$plane3d(mean(y),x.coef =0, y.coef = 0, lwd=2)

legend("topright", legend = c("Grupo 1","Grupo 2"), fill = cores[1:2])




jpeg(file=paste0("Resultados/", "teste","graficos3d.jpg"), width = 800, height = 800, quality = 100, pointsize = 30)

s3d

dev.off()
jpeg(file=paste0("Resultados/", nome.plot,"graficos3d.jpg"), width = 800, height = 800, quality = 100, pointsize = 30)


print(cloud(y~X[,2]*X[,3], col=cores[grupo], pch=c(8,19)[grupo], cex=0.5), split=c(1,1,2,2), more=T)
print(cloud(y~t*X[,3], col=cores[grupo], pch=c(8,19)[grupo],cex=0.5), split=c(1,2,2,2), more=T)
print(cloud(y~t*X[,2], col=cores[grupo], pch=c(8,19)[grupo],cex=0.5), split=c(2,1,2,2), more=T)

dev.off()

n<-length(y)
km<-kmeans(y,centers=g, iter.max=100, nstart=50)
p0<-km$size/sum(km$size)

# Mudar o label, o grupo que possui o maior valor de centro é o primeiro
# Segue em ordem decrescente
o.g<-order(km$centers, decreasing = T)
p0<-p0[o.g]
rotulos_reordenados<-sapply(km$cluster,FUN=function(y)which(o.g==y))


#--------------------------- SVM----------------------------
rotulos_reordenados <- as.factor(rotulos_reordenados)

# 3. Treinar SVM com rótulos do k-means
# opções de kernel: pode avaliar qual tem um comportamento melhor ou se são semelhantes #para o nosso cenário
#svm_model <- svm(rotulos_reordenados ~ y,kernel = "radial", probability = TRUE,fitted=T)
#svm_model <- svm(rotulos_reordenados ~ y, kernel = "radial", probability = TRUE, fitted=T)
svm_model <- svm(rotulos_reordenados ~ y, kernel = "polynomial",degree=5, probability = TRUE, fitted=TRUE)

label_svm <- svm_model$fitted

grupo<- apply(amostra$clu,2,FUN=function(x)which(as.logical(x))) # numero do grupo de cada observação


tabela<-table(grupo, rotulos_reordenados)
acuracia<-sum(diag(tabela))/sum(tabela)
acuracia


tabela2<-table(grupo, label_svm)
acuracia2<-sum(diag(tabela2))/sum(tabela2)
acuracia2

p0.2<-sapply(1:g,FUN=function(x)mean(label_svm==x))

#------------------------------------ Cenario 2 mal separados

nome.plot<-"c2.1_2_"

beta.verd<-list(c(3,4,6),c(2,3,-5))
sigma2.verd<-list(1,1)


arg.grupos<-list(g1=list(beta=beta.verd[[1]],curva=list(f="c1",a=-1,b=1,d=2),sigma2=sigma2.verd[[1]],intercepto=T),
                 g2=list(beta=beta.verd[[2]],curva=list(f="c2",a=-1,b=1,d=4),sigma2=sigma2.verd[[2]],intercepto=T),
                 intervalos=list(c(0,1),c(0,1)))

n=1000
alfas<-c(0.1,0.1)
amostra<-rMisUniPLM(n, pii, p=length(beta.verd[[1]]), arg=arg.grupos)

y=amostra$y
t=amostra$t
X=amostra$X

jpeg(file=paste0("Resultados/", nome.plot,"graficos3d.jpg"), width = 800, height = 800, quality = 100, pointsize = 30)

print(cloud(y~X[,2]*X[,3], col=cores[grupo], pch=c(8,19)[grupo], cex=0.5), split=c(1,1,2,2), more=T)
print(cloud(y~t*X[,3], col=cores[grupo], pch=c(8,19)[grupo], cex=0.5), split=c(2,1,2,2), more=T)
print(cloud(y~t*X[,2], col=cores[grupo], pch=c(8,19)[grupo], cex=0.5), split=c(1,2,2,2), more=T)

dev.off()
prop.gerada<-apply(amostra$clu,1,mean)


n<-length(y)
km<-kmeans(y,centers=g, iter.max=100, nstart=50)
p0<-km$size/sum(km$size)

# Mudar o label, o grupo que possui o maior valor de centro é o primeiro
# Segue em ordem decrescente
o.g<-order(km$centers, decreasing = T)
p0<-p0[o.g]
rotulos_reordenados<-sapply(km$cluster,FUN=function(y)which(o.g==y))


#--------------------------- SVM----------------------------
rotulos_reordenados <- as.factor(rotulos_reordenados)

# 3. Treinar SVM com rótulos do k-means
# opções de kernel: pode avaliar qual tem um comportamento melhor ou se são semelhantes #para o nosso cenário
#svm_model <- svm(rotulos_reordenados ~ y,kernel = "radial", probability = TRUE,fitted=T)
svm_model <- svm(rotulos_reordenados ~ y, kernel = "polynomial",degree=5, probability = TRUE, fitted=T)
#svm_model <- svm(label_kmeans ~ ., data = df, kernel = "polynomial",degree=3, probability = #TRUE)

label_svm <- svm_model$fitted
p0.2<-sapply(1:g,FUN=function(x)mean(label_svm==x))

grupo<-apply(amostra$clu,2,FUN=function(x)which(as.logical(x))) # numero do grupo de cada observação


tabela<-table(grupo, rotulos_reordenados)
acuracia<-sum(diag(tabela))/sum(tabela)
acuracia


tabela2<-table(grupo, label_svm)
acuracia2<-sum(diag(tabela2))/sum(tabela2)
acuracia2


p0.2<-sapply(1:g,FUN=function(x)mean(label_svm==x))

p0.2


#-------------------- Pouco separados-------------------
# Cenario 3.2
pii<-c(0.35,0.65)

beta.verd<-list(c(4,4,6),c(2,3,-5))
nome.amostra<-"Amostras/c3.2_2_"
nome.plot<-"c3.2_2_"


#------------------------------
sigma2.verd<-list(1,1)


arg.grupos<-list(g1=list(beta=beta.verd[[1]],curva=list(f="c1",a=-1,b=1,d=2),sigma2=sigma2.verd[[1]],intercepto=T),
                 g2=list(beta=beta.verd[[2]],curva=list(f="c2",a=-1,b=1,d=4),sigma2=sigma2.verd[[2]],intercepto=T),
                 intervalos=list(c(0,1),c(0,1)))


n=1000
alfas<-c(0.1,0.1)
amostra<-rMisUniPLM(n, pii, p=length(beta.verd[[1]]), arg=arg.grupos)

y=amostra$y
t=amostra$t
X=amostra$X

jpeg(file=paste0("Resultados/", nome.plot,"graficos3d.jpg"), width = 800, height = 800, quality = 100, pointsize = 30)

print(cloud(y~X[,2]*X[,3], col=cores[grupo], pch=c(8,19)[grupo], cex=0.5), split=c(1,1,2,2), more=T)
print(cloud(y~t*X[,3], col=cores[grupo], pch=c(8,19)[grupo], cex=0.5), split=c(2,1,2,2), more=T)
print(cloud(y~t*X[,2], col=cores[grupo], pch=c(8,19)[grupo], cex=0.5), split=c(1,2,2,2), more=T)

dev.off()


# 
# 
# # Plot
# plot3d(
#   x=t, y=X[,2], z=y, 
#   col = cores[grupo], 
#   type = 's', 
#   radius = .1,
#   xlab="t", ylab="X2", zlab="Y")
# 
# plot3d(
#   x=y, y=t, z=X[,3], 
#   col = cores[grupo], 
#   type = 's', 
#   radius = .1,
#   xlab="y", ylab="t", zlab="X3")
# 
# plot3d(
#   x=X[,2], y=X[,3], z=y, 
#   col = cores[grupo], 
#   type = 's', 
#   radius = .1,
#   xlab="X2", ylab="X3", zlab="Y")
# 
