########### Cross-validation - K-folds ( Hastie and Tibishirani )

nfolds<-function(n,folds=5){
  s = rep(floor(n/folds),folds)
  s[0:(n%%folds)] = s[0:(n%%folds)] + 1
  
  S = list()
  x = 1:n
  for(j in 1:folds){
    S[[j]] = sample(x,size=s[j])
    x = x[which(!x%in%S[[j]])]
  }
  return(S)
}

CV.NaoParametrico<-function(lambda,y,N,t0,ndx,bdeg,S){
  n=length(y)
  folds=length(S)
  CV = numeric(folds)
  for(i in 1:folds){
    Si=S[[i]]
    Ni=bspline(t0[-Si],ndx,bdeg)
    Dk=diff(diag(ncol(Ni)),differences=2)               
    K=t(Dk)%*%Dk 
    a<-solve(t(Ni)%*%Ni+lambda*K, t(Ni)%*%y[-Si])
    Yest<-N[Si,]%*%a
    CV[i] = sum(as.vector(y[Si]-Yest)^2)
  }
  #print(sum(CV)/n)
  return(sum(CV)/n)
}

CV.SemiParametrico<-function(alfa,y,t0,X,N,K,ndx,bdeg,S)
{
  n<-dim(X)[1]
  folds=length(S)
  CV = numeric(folds)
  for(i in 1:folds)
  {
    Si=S[[i]]
    Ni=bspline(t0[-Si],ndx,bdeg)
    Dk=diff(diag(ncol(Ni)),differences=2)               
    K=t(Dk)%*%Dk 
    theta<-estimar.semiParametrico(y=y[-Si],X=X[-Si,],N=Ni,K=K,alfa=alfa)
    Yest<-N[Si,]%*%theta$a+X[Si,]%*%theta$b
    CV[i] = sum(as.vector(y[Si]-Yest)^2)
  }
  return(sum(CV)/n)
}

CV.SemiParametricoSMN<-function(alfa,y,t0,X,N,K,ndx,bdeg,S,family)
{
  n<-dim(X)[1]
  folds=length(S)
  CV = numeric(folds)
  for(i in 1:folds)
  {
    Si=S[[i]]
    Ni=bspline(t0[-Si],ndx,bdeg)
    Ni.2<-bspline(t0[Si],ndx,bdeg)
    Dk=diff(diag(ncol(Ni)),differences=2)               
    K=t(Dk)%*%Dk
    theta<-EMSMNsp(y=y[-Si],X=X[-Si,],N=Ni,K=K,alpha=alfa,family=family)
    Yest<-Ni.2%*%theta$gamma+X[Si,]%*%theta$beta
    CV[i] = mean(as.vector(y[Si]-Yest)^2)
  }
  return(mean(CV))
}

CV.APL<-function(alfas,y,dados,X,ndx,S) #Testar a função
{
  n<-dim(X)[1]
  p<-ncol(X)
  folds=length(S)
  CV = numeric(folds)
  for(i in 1:folds)
  {
    Si=S[[i]]
    Curvas.treino<-vector("list",k)
    Curvas.teste<-vector("list",k)
    for(j in 1:k)
    {
      d<-data.frame(t=dados[-Si,p+j])
      ZZ<-smoothCon(s(t,bs="cr",k=ndx[j]),data=d,knots=NULL, absorb.cons=T)
      N<-ZZ[[1]]$X
      K<-(ZZ[[1]]$S)[[1]] 
      Curvas.treino[[j]]<-list(N=N,K=K,ndx=ndx[j])
      
      d<-data.frame(t=dados[Si,p+j])
      ZZ2<-smoothCon(s(t,bs="cr",k=ndx[j]),data=d,knots=NULL, absorb.cons=T)
      N2<-ZZ2[[1]]$X
      Curvas.teste[[j]]<-list(N=N2)
    }
    
    theta<-EAPL(y=y[-Si],X=X[-Si,], Curvas.treino, alfas=alfas)
    cat("alfas:",alfas)
    Curvas.treino<-vector("list",k)
    Curvas.teste<-vector("list",k)
    for(i in 1:k)
    {
      d<-data.frame(t=dados[-Si,i])
      ZZ<-smoothCon(s(t,bs="cr",k=ndx[i]),data=d,knots=NULL, absorb.cons=T)
      N<-ZZ[[1]]$X
      K<-(ZZ[[1]]$S)[[1]] 
      Curvas.treino[[i]]<-list(N=N,K=K,ndx=ndx[i])
      
      d<-data.frame(t=dados[Si,i])
      ZZ2<-smoothCon(s(t,bs="cr",k=ndx1),data=d,knots=NULL, absorb.cons=T)
      N2<-ZZ2[[1]]$X
      Curvas.teste[[i]]<-list(N=N2)
    }
    theta<-EAPL(y=y[-Si],X=X[-Si,], Curvas.treino, alfas=alfas)
    vetor<-1:k
    aux<-function(j)
    {
      Curvas.teste[[j]]$N%*%theta$a[[j]]
    }
    N.gamma<-sapply(vetor,FUN=aux)
    sum.N.Gamma<-apply(N.gamma,1,sum)
    Yest<-X[Si,]%*%theta$b+sum.N.Gamma
    CV[i] = mean(as.vector(y[Si]-Yest)^2)
  }
  return(mean(CV))
}

