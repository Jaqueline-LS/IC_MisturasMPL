rm(list = ls(all = TRUE))
library("Matrix")
library("mgcv")
library("dplyr")
library("kableExtra") 
source("funcoes/funcoes.R")
source("funcoes/funcoesMisturas.R")
source("funcoes/GerarAmostras_geral.R") 
cores<-c("#BAF3DE","#C9E69E","#FF9B95","#FFC29A")
cores2<-c("#1d91c0","#ff5858")
par<-par(pch=19)



gerarAmostra<-function(n, beta, curva, sigma2)
{
  x1<-runif(n,0,1)
  x2<-runif(n,0,1)
  t<-runif(n,-1,1)
  erro<-rnorm(n,mean=0,sd=sqrt(sigma2))
  X<-cbind(rep(1,n),x1,x2)
  f<-get(curva)
  y<-X%*%beta+f(t, d=2)+erro
  dados<-list(X=X,t=t,y=y)
  return(dados)
}

gera.amostras1<-function(n)
{
  teptetam<-vector("list",M)
  for (i in 1:M)
  {
    aux=0
    while (aux==0){
      amostra<-gerarAmostra(n, beta.verd, curva="c1", sigma2=2)
      theta<-try(EMisUniPLM(g=g,alfas=alfas,y=amostra$y, t=amostra$t, X=amostra$X),TRUE)
      t<-amostra$t
      d<-data.frame(t=t)
      ZZ<-smoothCon(s(t,bs="ps",m=c(2,3)),data=d,knots=NULL,absorb.cons=T) 
      # absorb.cons = T ( restrição de idenficabilidade )
      # bs="ps",m=c(2,3) P-splines, B-splines cúbicos com penalização de ordem 2
      N<-ZZ[[1]]$X
      K<-(ZZ[[1]]$S)[[1]]
      theta2<-try(estimar.semiParametrico(y=amostra$y,X=amostra$X,N,K,alfa=alfas,precision=1e-5),TRUE)
      
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
    if(g==1)
    {
      sd.p<-NULL
    }else{
      sd.p<-sqrt(diag(cov.theta)[1:(g-1)])
    }
    sd.b<-sqrt(diag(cov.theta)[g:((g-1)+g*p)])
    
    #sd.a<-sqrt(diag(cov.theta)[(g+g*p):((g-1)+g*p+g*q)])
    
    sd.sigma2<-sqrt(diag(cov.theta)[(g+g*p+g*q):((g-1)+g*p+g*q+g)])
    
    #erro.padrao<-list(sd.p=sd.p,sd.b=sd.b, sd.a=sd.a,sd.sigma2=sd.sigma2, ginv=aux)
    
    erro.padrao<-list(sd.p=sd.p,sd.b=sd.b, sd.sigma2=sd.sigma2)
    erro.padrao2<-MI.semiParametrico(X=amostra$X,N,alfa=alfas,sigma2.hat=theta2$sigma2,K)

    teptetam[[i]]<-list(t=amostra$t,t1=theta,t2=theta2, erro.padrao, erro.padrao2, y=amostra$y,X=amostra$X,alfas=alfas)
    
  }
  saveRDS(teptetam, file=paste0(nome.amostra,n,".RDS"))
}

estimar.semiParametrico<-function(y,X,N,K,alfa,precision=1e-5)
{
  n<-dim(X)[1]
  # valor inicial
  b.beta<- solve ( t(X)%*%X , t(X)%*%(y) )
  
  sigma2.hat<-as.numeric((t(y-(X%*%b.beta)) %*% (y-(X%*%b.beta))) / n)
  
  a.gama<-solve( t(N)%*%N + (alfa*sigma2.hat*K) , t(N)%*%(y-(X%*%b.beta)) )
  
  theta0<-matrix(c(b.beta, a.gama, sigma2.hat))
  # convergência
  criterio<-precision+1
  while(criterio>precision)
  {
    b.beta<-solve( t(X)%*%X , t(X)%*%(y-(N%*%a.gama) ))
    sigma2.hat<-as.numeric(( t( y - (X%*%b.beta) - (N%*%a.gama) ) %*% ( y - (X%*%b.beta) - (N%*%a.gama)) ) / n)
    a.gama<-solve(t(N)%*%N + ((alfa*sigma2.hat)*K), t(N)%*%(y-(X%*%b.beta)))
    theta<-matrix(c(b.beta, a.gama, sigma2.hat))
    criterio<-sqrt(sum((theta-theta0)^2)) #norma l2
    theta0<-theta
  }
  return(list(b=b.beta,a=a.gama,sigma2=sigma2.hat ))
}



#---------Função que monta a MI esperada, retorna os erros das estimativas
MI.semiParametrico<-function(X,N,alfa,sigma2.hat,K)
{
  n<-dim(X)[1]
  p<-dim(X)[2]
  q<-dim(N)[2]
  MI<-matrix(0, nrow=(p+q+1), ncol=(p+q+1)) 
  MI[1:p,1:p]<-t(X)%*%X
  MI[1:p,(p+1):(p+q)]<-t(X)%*%N
  MI[(p+1):(p+q),1:p]<- t(N)%*%X
  MI[(p+1):(p+q),(p+1):(p+q)]<- t(N)%*%N + alfa*sigma2.hat*K
  MI[p+q+1, p+q+1]<-n/(2*sigma2.hat)
  MI<-(1/sigma2.hat)*MI
  teste<-try(solve(MI), silent = T)
  aux<-sum(class(teste)=="try-error")
  if(aux==0)
  {
    cov.theta<-teste
  }
  else
  {
    #cov.theta<-chol2inv(chol(MI))
    cov.theta<-MASS::ginv(MI)
  }
  var.b<-diag(cov.theta[1:p,1:p])
  sd.b<-sqrt(var.b)
  sd.sigma2<-sqrt(cov.theta[p+q+1, p+q+1])
  cov.y.hat<-N%*%cov.theta[(p+1):(p+q),(p+1):(p+q)]%*%t(N)
  ep.curva<-sqrt(diag(cov.y.hat))
  #return(list(sdEmp.b=sd.b, sdEmp.sigma2=sd.sigma2, ep.curva=ep.curva))
  return(list(sdEmp.b=sd.b, sdEmp.sigma2=sd.sigma2))
  
}


tabela.n<-function(n) # Função para ler os arquivos e pegar os valores estimados para as M replicas 
{
  arquivo<-paste0(nome.amostra,n,".RDS")
  amostras<-readRDS(arquivo)
  
  
  # Betas EMis
  betas<-lapply(1:M,FUN=function(x){amostras[[x]]$t1$theta$b})
  
  betas1<-sapply(1:M,FUN=function(x){betas[[x]][[1]]})
  b1.medios<-apply(betas1, MARGIN = 1, FUN=mean)
  
  
  sd.b1<-apply(betas1, MARGIN = 1, FUN=sd)
  
  # Betas funçao do semiparametrico
  betas<-lapply(1:M,FUN=function(x){amostras[[x]]$t2$b})
  
  betas1<-sapply(1:M,FUN=function(x){betas[[x]]})
  b2.medios<-apply(betas1, MARGIN = 1, FUN=mean)
  
  
  sd.b2<-apply(betas1, MARGIN = 1, FUN=sd)
  
  # Sigmas
  sigmas2<-sapply(1:M,FUN=function(x){amostras[[x]][[2]]$theta$sd2})
  sigmas2.medios<-mean(sigmas2)
  sd.s2<-sd(sigmas2)
  
  sigmas2.2<-sapply(1:M,FUN=function(x){amostras[[x]]$t2$sigma2})
  sigmas2.2.medios<-mean(sigmas2.2)
  sd.s2.2<-sd(sigmas2.2)
  
  
  
  sd.Emp<-lapply(1:M,FUN=function(x){amostras[[x]][[4]]})
  sd.Emp<-sapply(1:M, FUN=function(x){unlist(sd.Emp[[x]])})
  sd.Emp<-apply(sd.Emp, MARGIN = 1,mean)
  
  
  sd.Emp.2<-lapply(1:M,FUN=function(x){amostras[[x]][[5]]})
  sd.Emp.2<-sapply(1:M, FUN=function(x){unlist(sd.Emp.2[[x]])})
  sd.Emp.2<-apply(sd.Emp.2, MARGIN = 1,mean)
  
  tabela<-data.frame(
      estimativa1=c(b1.medios,sigmas2.medios),c(b2.medios,sigmas2.2.medios), 
      ep=c(sd.b1, sd.s2), c(sd.b2, sd.s2.2), ep.MI=c(sd.Emp),c=(ep.MI2=sd.Emp.2))
  
  
  
  colnames(tabela)<-c("$\\hat{\\theta}$","$\\hat{\\theta_{2}}$", "sd", "sd_2","sd.emp", "sd.emp2")
  
  
  plot(x="",xlim = c(-1,1), ylim=c(-6,6), 
       main=paste0("n = ",n), 
       cex.main=1, xlab="t", ylab="f(t)", cex.axis=1)
  
  for(j in 1:M)
  {
    N<-amostras[[j]][[2]]$N
    t<-amostras[[j]][[1]]
    y_hat<-as.numeric(N%*%amostras[[j]]$t1$theta$a[[1]])
    dados<-data.frame(t,y_hat)
    dados<-dados[order(dados$t, decreasing=FALSE),]
    lines(dados$t,dados$y_hat, col=cores[1], lwd=2)
  }
  curve(c1(x,d=2),from=-1,to=1, lwd=1,add=T)
  
  for(j in 1:M)
  {
    N<-amostras[[j]][[2]]$N
    t<-amostras[[j]][[1]]
    y_hat<-as.numeric(N%*%amostras[[j]]$t2$a)
    dados<-data.frame(t,y_hat)
    dados<-dados[order(dados$t, decreasing=FALSE),]
    lines(dados$t,dados$y_hat, col=cores[3], lty=2, lwd=2)
  }
  curve(c1(x,d=2),from=-1,to=1, lwd=1,add=T)
  
  legend("topright", legend=c("Algoritmo de misturas", "Algoritmo PLM"), bty = "n",lwd = 3, col=cores[c(1,3)],lty=c(1,2), cex=1)
  
  return(tabela)
}



#------------------- Geração das amostras pelos 2 métodos--------------------

n=500
beta.verd<-c(8,4,6)
sigma2<-2


alfas<-1
g=1
M=20
nome.amostra<-"comparacao.amostra1"

#gera.amostras1(1000)
tabela.n(1000)
