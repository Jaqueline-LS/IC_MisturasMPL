gera.amostras<-function(n)
{
  teptetam<-vector("list",M)
  for (i in 1:M)
  {
    aux=0
    while (aux==0){
      amostra<-rMisUniPLM(n, pii,p=length(beta.verd[[1]]), arg=arg.grupos) 
     # theta<-try(EMisUniPLM(g=g,alfas=alfas,y=amostra$y, t=amostra$t, X=amostra$X),TRUE)
      theta<-try(EMisUniPLM(g=g,alfas=alfas,y=amostra$y, t=amostra$t, X=amostra$X, clu = amostra$clu),TRUE)
      aux1<-sum(class(theta) != "try-error")
      if(aux1!=0)
      {
        MI<-MI.MisUniPLM(y=amostra$y,X=amostra$X,alfas=alfas,theta)
        teste<-try(solve(MI), silent = TRUE)
        aux2<-sum(class(teste) != "try-error")
        aux<- aux1 & aux2
      }else{
        aux<-0
      }
    }
    cov.theta<-teste
    p<-ncol(amostra$X)
    n<-nrow(amostra$X)
    q<-ncol(theta$N)
    g<-length(theta$theta$p)+1
    
    sd.p<-sqrt(diag(cov.theta)[1:(g-1)])
    sd.b<-sqrt(diag(cov.theta)[g:((g-1)+g*p)])
    
    #sd.a<-sqrt(diag(cov.theta)[(g+g*p):((g-1)+g*p+g*q)])
    
    sd.sigma2<-sqrt(diag(cov.theta)[(g+g*p+g*q):((g-1)+g*p+g*q+g)])
    
    #erro.padrao<-list(sd.p=sd.p,sd.b=sd.b, sd.a=sd.a,sd.sigma2=sd.sigma2, ginv=aux)
    
    erro.padrao<-list(sd.p=sd.p,sd.b=sd.b, sd.sigma2=sd.sigma2)
    
    teptetam[[i]]<-list(t=amostra$t, theta, erro.padrao, y=amostra$y,X=amostra$X,alfas=alfas)
    
  }
  saveRDS(teptetam, file=paste0(nome.amostra,n,".RDS"))
}


tabela.n<-function(n)
{
  arquivo<-paste0(nome.amostra,n,".RDS")
  amostras<-readRDS(arquivo)
  
  # Recuperando as proporções estimadas
  p<-lapply(1:M,FUN=function(x){amostras[[x]][[2]]$theta$p})
  
  p1<-sapply(1:M,FUN=function(x){p[[x]][[1]]})
  p1.medias<-mean(p1)
  
  sd.p<-sd(p1)
  
  # Betas
  betas<-lapply(1:M,FUN=function(x){amostras[[x]][[2]]$theta$b})
  
  betas1<-sapply(1:M,FUN=function(x){betas[[x]][[1]]})
  b1.medios<-apply(betas1, MARGIN = 1, FUN=mean)
  
  betas2<-sapply(1:M,FUN=function(x){betas[[x]][[2]]})
  b2.medios<-apply(betas2, MARGIN = 1, FUN=mean)
  
  sd.b1<-apply(betas1, MARGIN = 1, FUN=sd)
  sd.b2<-apply(betas2, MARGIN = 1, FUN=sd)
  
  # Sigmas
  sigmas2<-sapply(1:M,FUN=function(x){amostras[[x]][[2]]$theta$sd2})
  sigmas2.medios<-apply(sigmas2, 1, mean)
  
  sd.s2<-apply(sigmas2, 1, sd)
  
  sd.Emp<-lapply(1:M,FUN=function(x){amostras[[x]][[3]]})
  sd.Emp<-sapply(1:M, FUN=function(x){unlist(sd.Emp[[x]])})
  sd.Emp<-apply(sd.Emp, MARGIN = 1,mean)
  
  tabela<-data.frame(
    estimativa=c(p1.medias,b1.medios,b2.medios,sigmas2.medios), 
    ep=c(sd.p,sd.b1,sd.b2, sd.s2), ep.MI=sd.Emp)
  
  colnames(tabela)<-c("$\\hat{\\theta}$", "sd", "sd.emp")
  
  acuracias<-sapply(1:M,FUN=function(x){amostras[[x]][[2]]$acuracia})
  acuracias2<-sapply(1:M,FUN=function(x){amostras[[x]][[2]]$acuracia2})
  cat("A acuracia do k-means no agrupamento inicial foi: ", mean(acuracias)*100)
  cat("A acuracia depois do SVM foi: ", mean(acuracias2)*100)

  
  plot(x="",xlim = c(-1,1), ylim=c(-6,6), 
       main=paste0("n = ",n), 
       cex.main=1, xlab="t", ylab="f(t)", cex.axis=1)
  
  for(j in 1:M)
  {
    N<-amostras[[j]][[2]]$N
    t<-amostras[[j]][[1]]
    y_hat<-as.numeric(N%*%amostras[[j]][[2]]$theta$a[[1]])
    dados<-data.frame(t,y_hat)
    dados<-dados[order(dados$t, decreasing=FALSE),]
    lines(dados$t,dados$y_hat, col=cores[1])
  }
  curve(c1(x,d=2),from=-1,to=1, lwd=1,add=T)
  
  
  for(j in 1:M)
  {
    N<-amostras[[j]][[2]]$N
    t<-amostras[[j]][[1]]
    y_hat<-as.numeric(N%*%amostras[[j]][[2]]$theta$a[[2]])
    dados<-data.frame(t,y_hat)
    dados<-dados[order(dados$t, decreasing=FALSE),]
    lines(dados$t,dados$y_hat, col=cores[2])
  }
  curve(c2(x,d=4),from=-1,to=1, lwd=1,add=T)
  
  legend("topright", legend=c("Grupo 1","Grupo 2"), bty = "n",lwd = 3, col=cores[1:2], cex=1)
  
  return(tabela)
}
