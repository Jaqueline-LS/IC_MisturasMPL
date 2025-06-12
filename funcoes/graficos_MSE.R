
# Gráfico MSE dos Betas
calculaMSE<-function(j,p)
{
  mse.betas<-matrix(0,ncol=length(sizes), nrow=500)
  curva<-arg.grupos[[j]]$curva
  f<-curva$f
  d<-curva$d
  curva.f<-get(f)
  for(i in seq_along(sizes))
  {
    arquivo<-paste0(nome.amostra,sizes[i],".RDS")
    amostras<-readRDS(arquivo)
    
    betas<-lapply(1:M,FUN=function(x){amostras[[x]][[2]]$theta$b})
    
    betas1<-sapply(1:M,FUN=function(x){betas[[x]][[j]][p,]})

    mse.betas[,i]<-((beta.verd[[j]][p]-betas1)^2)
  }
  return(mse.betas)
}
mseBetas<-sapply(1:3,FUN=function(x)sapply(1:2, FUN=calculaMSE,p=x,simplify = F), simplify = F)
jpeg(file=paste0("Resultados/", nome.plot,"mse.b.jpg"), width = 1500, height = 800, quality = 100, pointsize = 20)
par(mfrow=c(1,3), mar=c(3,2,2,1))
for(i in 1:3)
{
  x<-seq_along(sizes)
  titulo<-c(expression(paste("Gráfico do MSE de ", beta[0])), 
            expression(paste("Gráfico do MSE de ", beta[1])),
            expression(paste("Gráfico do MSE de ", beta[2])))
  plot(x, colMeans(mseBetas[[i]][[1]]), type = "b", pch=19, col=cores[1], ylim=c(0,2), 
       ylab="MSE", xlab = "", xaxt = "n",
       main=titulo[i], lwd=2)
  axis(1, at = x, labels = paste0("n = ",sizes)) # Personaliza o eixo X com intervalos de 3
  sd<-apply(mseBetas[[i]][[1]], 2,sd)
  arrows(x, y0=colMeans(mseBetas[[i]][[1]])-sd, x, y1 = colMeans(mseBetas[[i]][[1]])+sd,
         col = cores[1], lwd = 2, code=3, angle = 90, length = 0.10)
  points(x+0.1, colMeans(mseBetas[[i]][[2]]), pch=19, col=cores[2])
  lines(x+0.1, colMeans(mseBetas[[i]][[2]]), col=cores[2], lwd=2)
  sd<-apply(mseBetas[[i]][[2]], 2,sd)
  arrows(x+0.1, y0=colMeans(mseBetas[[i]][[2]])-sd, x+0.1, y1 = colMeans(mseBetas[[i]][[2]])+sd,
         col = cores[2], lwd = 2, code=3, angle = 90, length = 0.10)
  legend("topright", legend = c("grupo 1", "grupo 2"), lty = 1,lwd=3, col = cores[1:2], cex=0.8, bty="n")
  
}
dev.off()
#----------------------- curvas------------------------
ASE.n<-function(j)
{
  ASE.k<-matrix(0,ncol=length(sizes), nrow=500)
  curva<-arg.grupos[[j]]$curva
  f<-curva$f
  d<-curva$d
  curva.f<-get(f)
  for(i in seq_along(sizes))
  {
    arquivo<-paste0(nome.amostra,sizes[i],".RDS")
    amostras<-readRDS(arquivo)
    gammas<-lapply(1:500,FUN=function(x){amostras[[x]][[2]]$theta$a})
    gammas1<-sapply(1:500,FUN=function(x){gammas[[x]][[j]]})
   
    aux.g.t.hat<-function(x)
    {
      N<-amostras[[x]][[2]]$N 
      return((N%*%gammas1[,x] - curva.f(amostras[[x]]$amostra$t,d=d))^2)
    }
    
    g.t.hat<-sapply(1:500,FUN=aux.g.t.hat)
    ASE.k[,i]<-apply(g.t.hat, 2,mean)
  }
  return(ASE.k)
}



ASE.j<-sapply(1:2, FUN=ASE.n, simplify = F)

#---------------------------------------Segunda imagem----------------

jpeg(file=paste0("Resultados/", nome.plot,"mse.outros.jpg"), width = 1500, height = 800, quality = 100, pointsize = 20)
par(mfrow=c(1,3), mar=c(3,2,2,1))
#------------- Grafico do ASE das curvas----------------------------
x<-seq_along(sizes)
titulo<-c(expression(paste("Gráfico do ASE das curvas")))
plot(x, colMeans(ASE.j[[1]]), type = "b", pch=19, col=cores[1], ylim=c(0,5), 
ylab="ASE", xlab = "", xaxt = "n",
main=titulo, lwd=2)
axis(1, at = x, labels = paste0("n = ",sizes)) # Personaliza o eixo X com intervalos de 3
sd<-apply(ASE.j[[1]], 2,sd)
arrows(x, y0=c(colMeans(ASE.j[[1]])-sd), x, y1 = c(colMeans(ASE.j[[1]])+sd),
           col = cores[1], lwd = 2, code=3, angle = 90, length = 0.10)
points(x+0.05, colMeans(ASE.j[[2]]), pch=19, col=cores[2])
lines(x+0.05, colMeans(ASE.j[[2]]), pch=19, col=cores[2], lwd=2)
sd<-apply(ASE.j[[2]], 2,sd)
arrows(x+0.05, y0=c(colMeans(ASE.j[[2]])-sd), x+0.05, y1 = c(colMeans(ASE.j[[2]])+sd),
         col = cores[2], lwd = 2, code=3, angle = 90,length = 0.10)
  
legend("topright", legend = c("grupo 1", "grupo 2"), lty = 1,lwd=3, col = cores[1:2], cex=0.8, bty="n")
  



# MSE proporções

MSE.p<-function(j)
{
  mse.k<-matrix(0,ncol=length(sizes), nrow=500)
  
  for(i in seq_along(sizes))
  {
    arquivo<-paste0(nome.amostra,sizes[i],".RDS")
    amostras<-readRDS(arquivo)
    p<-lapply(1:M,FUN=function(x){amostras[[x]][[2]]$theta$p})
    mse.k[,i]<-sapply(1:500,FUN = function(x){(pii[j] - p[[x]][j])^2})
  }
  return(mse.k)
}

mse.j<-sapply(1:2, FUN=MSE.p, simplify = F)

x<-seq_along(sizes)
titulo<-c(expression(paste("Gráfico do MSE ",p[j])))
plot(x, colMeans(mse.j[[1]]), type = "b", pch=19, col=cores[1], ylim=c(0,0.1), 
     ylab="MSE", xlab = "", xaxt = "n",
     main= titulo, lwd=2)
axis(1, at = x, labels = paste0("n = ",sizes)) # Personaliza o eixo X com intervalos de 3
sd<-apply(mse.j[[1]], 2,sd)
arrows(x, y0=colMeans(mse.j[[1]])-sd, x, y1 = colMeans(mse.j[[1]])+sd,
       col = cores[1], lwd = 2, code=3, angle = 90, length = 0.10)
points(x+0.05, colMeans(mse.j[[2]]), pch=19, col=cores[2])
lines(x+0.05, colMeans(mse.j[[2]]), col=cores[2], lwd=2)
sd<-apply(mse.j[[2]], 2,sd)
arrows(x+0.05, y0=colMeans(mse.j[[2]])-sd, x+0.05, y1 = colMeans(mse.j[[2]])+sd,
       col = cores[2], lwd = 2, code=3, angle = 90, length = 0.10)
legend("topright", legend = c("grupo 1", "grupo 2"), lty = 1,lwd=3, col = cores[1:2], cex=0.8, bty="n")



# MSE sigma2

MSE.sd<-function(j)
{
  mse.k<-matrix(0,ncol=length(sizes), nrow=500)
  
  for(i in seq_along(sizes))
  {
    arquivo<-paste0(nome.amostra,sizes[i],".RDS")
    amostras<-readRDS(arquivo)
    sigmas2<-sapply(1:M,FUN=function(x){amostras[[x]][[2]]$theta$sd2})
    mse.k[,i]<-sapply(1:500,FUN = function(x){(as.numeric(sigma2.verd[[j]]) - sigmas2[j,x])^2})
  }
  return(mse.k)
}

mse.j<-sapply(1:2, FUN=MSE.sd, simplify = F)

x<-seq_along(sizes)
titulo<-c(expression(paste("Gráfico do MSE ", sigma[j]^2)))
plot(x, colMeans(mse.j[[1]]), type = "b", pch=19, col=cores[1], ylim=c(0,4), 
     ylab="MSE", xlab = "", xaxt = "n",
     main= titulo, lwd=2)
axis(1, at = x, labels = paste0("n = ",sizes)) # Personaliza o eixo X com intervalos de 3
sd<-apply(mse.j[[1]], 2,sd)
arrows(x, y0=colMeans(mse.j[[1]])-sd, x, y1 = colMeans(mse.j[[1]])+sd,
       col = cores[1], lwd = 2, code=3, angle = 90, length = 0.10)
points(x+0.05, colMeans(mse.j[[2]]), pch=19, col=cores[2])
lines(x+0.05, colMeans(mse.j[[2]]), col=cores[2], lwd=2)
sd<-apply(mse.j[[2]], 2,sd)
arrows(x+0.05, y0=colMeans(mse.j[[2]])-sd, x+0.05, y1 = colMeans(mse.j[[2]])+sd,
       col = cores[2], lwd = 2, code=3, angle = 90, length = 0.10)
legend("topright", legend = c("grupo 1", "grupo 2"), lty = 1,lwd=3, col = cores[1:2], cex=0.8, bty="n")

dev.off()




