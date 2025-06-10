gera.amostras<-function(n)
{
  teptetam<-vector("list",M)
  for (i in 1:M)
  {
    aux=0
    while (aux==0){
      amostra<-rMisUniPLM(n, pii,p=length(beta.verd[[1]]), arg=arg.grupos) 
      theta<-try(EMisUniPLM(g=g,alfas=alfas,y=amostra$y, t=amostra$t, X=amostra$X),TRUE)
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
    
    erro.padrao<-list(sd.p=sd.p,sd.b=sd.b, sd.sigma2=sd.sigma2)
    
    teptetam[[i]]<-list(amostra=amostra,saida.EM=theta,erro.padrao=erro.padrao,alfas=alfas)
    
  }
  saveRDS(teptetam, file=paste0(nome.amostra,n,".RDS"))
}



gera.amostras.SVM<-function(n)
{
  teptetam<-vector("list",M)
  for (i in 1:M)
  {
    aux=0
    while (aux==0){
      amostra<-rMisUniPLM(n, pii,p=length(beta.verd[[1]]), arg=arg.grupos) 
      theta<-try(EMisUniPLM.SVM(g=g,alfas=alfas,y=amostra$y, t=amostra$t, X=amostra$X, clu = amostra$clu),TRUE)
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
    
    teptetam[[i]]<-list(t=amostra$t, theta, erro.padrao, y=amostra$y,X=amostra$X,alfas=alfas, clu=amostra$clu)
    
  }
  saveRDS(teptetam, file=paste0(nome.amostra,n,".RDS"))
}




tabela.n<-function(n) # Função para ler os arquivos e pegar os valores estimados para as M replicas 
{
  arquivo<-paste0(nome.amostra,n,".RDS")
  amostras<-readRDS(arquivo)
  
  # Recuperando as proporções estimadas
  p<-lapply(1:M,FUN=function(x){amostras[[x]]$saida.EM$theta$p})
  
  p1<-sapply(1:M,FUN=function(x){p[[x]][[1]]})
  p1.medias<-mean(p1)
  
  sd.p<-sd(p1)
  
  # Betas
  betas<-lapply(1:M,FUN=function(x){amostras[[x]]$saida.EM$theta$b})
  
  betas1<-sapply(1:M,FUN=function(x){betas[[x]][[1]]})
  b1.medios<-apply(betas1, MARGIN = 1, FUN=mean)
  
  betas2<-sapply(1:M,FUN=function(x){betas[[x]][[2]]})
  b2.medios<-apply(betas2, MARGIN = 1, FUN=mean)
  
  sd.b1<-apply(betas1, MARGIN = 1, FUN=sd)
  sd.b2<-apply(betas2, MARGIN = 1, FUN=sd)
  
  # Sigmas
  sigmas2<-sapply(1:M,FUN=function(x){amostras[[x]]$saida.EM$theta$sd2})
  sigmas2.medios<-apply(sigmas2, 1, mean)
  
  sd.s2<-apply(sigmas2, 1, sd)
  
  sd.Emp<-lapply(1:M,FUN=function(x){amostras[[x]]$erro.padrao})
  sd.Emp<-sapply(1:M, FUN=function(x){unlist(sd.Emp[[x]])})
  sd.Emp<-apply(sd.Emp, MARGIN = 1,mean)
  
  tabela<-data.frame(
    estimativa=c(p1.medias,b1.medios,b2.medios,sigmas2.medios), 
    ep=c(sd.p,sd.b1,sd.b2, sd.s2), ep.MI=sd.Emp)
  
  colnames(tabela)<-c("$\\hat{\\theta}$", "sd", "sd.emp")

  # acuracias<-sapply(1:M,FUN=function(x){amostras[[x]][[2]]$acuracia})
  # acuracias2<-sapply(1:M,FUN=function(x){amostras[[x]][[2]]$acuracia2})
  # boxplot(acuracias,acuracias2, col=cores[1], main="Acurácia do agrupamento inicial das das 500 réplicas n=2000",xaxt="n")
  # axis(1,at=c(1,2) ,labels = c("K-Means","SVM") )# Personaliza o eixo X com intervalos de 3
  # 
  
  
  
  plot(x="",xlim = c(-1,1), ylim=c(-6,6), 
       main=paste0("n = ",n), 
       cex.main=1, xlab="t", ylab="f(t)", cex.axis=1)
  
  for(j in 1:M)
  {
    N<-amostras[[j]]$saida.EM$N
    t<-amostras[[j]]$amostra$t
    y_hat<-as.numeric(N%*%amostras[[j]]$saida.EM$theta$a[[1]])
    dados<-data.frame(t,y_hat)
    dados<-dados[order(dados$t, decreasing=FALSE),]
    lines(dados$t,dados$y_hat, col=cores[1])
  }
  curva<-get(arg.grupos[[1]]$curva$f)
  d<-arg.grupos[[1]]$curva$d
  a<-arg.grupos[[1]]$curva$a
  b<-arg.grupos[[1]]$curva$b
  curve(curva(x,d),a,b, lwd=1,add=T)
  
  
  for(j in 1:M)
  {
    N<-amostras[[j]]$saida.EM$N
    t<-amostras[[j]]$amostra$t
    y_hat<-as.numeric(N%*%amostras[[j]][[2]]$theta$a[[2]])
    dados<-data.frame(t,y_hat)
    dados<-dados[order(dados$t, decreasing=FALSE),]
    lines(dados$t,dados$y_hat, col=cores[2])
  }
  curva<-get(arg.grupos[[2]]$curva$f)
  d<-arg.grupos[[2]]$curva$d
  a<-arg.grupos[[2]]$curva$a
  b<-arg.grupos[[2]]$curva$b
  curve(curva(x,d),a,b, lwd=1,add=T)
  
  legend("topright", legend=c("Grupo 1","Grupo 2"), bty = "n",lwd = 3, col=cores[1:2], cex=1)
  
  return(tabela)
}
