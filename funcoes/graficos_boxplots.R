
jpeg(file=paste0("Resultados/", nome.plot,"boxplots-outros.jpg"), width = 800, height = 800, quality = 100, pointsize = 20)

par(mfrow=c(g,g), mar=c(1,2,1.5,1.5))
titulo<-c(expression(paste(sigma[1]^2)), 
          expression(paste(sigma[2]^2)),
          expression(paste(sigma[3]^2)))
recupera.sigmas<-function(i,j)
{
  arquivo<-paste0(nome.amostra,sizes[i],".RDS")
  amostras<-readRDS(arquivo)
  sigmas2<-sapply(1:M,FUN=function(x){amostras[[x]][[2]]$theta$sd2})
  return(sigmas2[j,])
}

recupera.proporcoes<-function(i,j)
{
  arquivo<-paste0(nome.amostra,sizes[i],".RDS")
  amostras<-readRDS(arquivo)
  p<-lapply(1:M,FUN=function(x){amostras[[x]][[2]]$theta$p[j]})
  return(unlist(p))
}




l<-lapply(1:g,FUN=function(j)lapply(seq_along(sizes),FUN=recupera.sigmas,j))
for(j in 1:g)
{
  l.j<-l[[j]]
  boxplot(l.j,boxwex=0.2,main=titulo[j], pch=1,at=c(1,1.25,1.50,1.75,2), col=cores,xaxt="n", ylim=c(min(unlist(l.j)),max(unlist(l.j))))
  abline(h=sigma2.verd[j], lty=2, lwd=2)
  if(j==1)legend("topright",legend = paste0("n = ",sizes), col=cores, bty="n", fill=cores, cex=0.8)
}


titulo<-c(expression(paste(p[1])), 
          expression(paste(p[2])),
          expression(paste(p[3])))
l2<-lapply(1:(g-1),FUN=function(j)lapply(seq_along(sizes),FUN=recupera.proporcoes,j))

for(j in 1:(g-1))
{
  l.j<-l2[[j]]
  boxplot(l.j,boxwex=0.2,main=titulo[j], pch=1,at=c(1,1.25,1.50,1.75,2), col=cores,xaxt="n", ylim=c(min(unlist(l.j)),max(unlist(l.j))))
  abline(h=pii[j], lty=2, lwd=2)
  #legend("topright",legend = paste0("n = ",sizes), col=cores, bty="n", fill=cores, cex=0.5)
}


dev.off()

# ------------------------------- segunda imagem---------
jpeg(file=paste0("Resultados/", nome.plot,"boxplots-betas.jpg"), width = 800, height = 1200, quality = 100, pointsize = 20)
par(mfrow=c(3,2), mar=c(1.5,2,1.5,1.5))

titulo<-c(expression(paste(beta[10])), 
          expression(paste(beta[11])),
          expression(paste(beta[12])),
          expression(paste(beta[20])), 
          expression(paste(beta[21])),
          expression(paste(beta[22])))
recupera.betas<-function(i,j,z)
{
  arquivo<-paste0(nome.amostra,sizes[i],".RDS")
  amostras<-readRDS(arquivo)
  betas<-lapply(1:M,FUN=function(x){amostras[[x]][[2]]$theta$b[[j]][z,]})
  return(unlist(betas))
}

for(x in 1:3)
{
  l3<-lapply(1:g,FUN=function(j)lapply(seq_along(sizes),FUN=recupera.betas,j=j,z=x))
  for(j in 1:g)
  {
    l.p<-l3[[j]]
    boxplot(l.p,boxwex=0.2,main=titulo[((j-1)*3)+x], pch=1,at=c(1,1.25,1.50,1.75,2), col=cores,xaxt="n", ylim=c(min(unlist(l.p)),max(unlist(l.p))))
    abline(h=beta.verd[[j]][x], lty=2, lwd=2)
    if(x==1 && j==1) legend("topright",legend = paste0("n = ",sizes), col=cores, bty="n", fill=cores, cex=0.8)
  }
}


dev.off()
