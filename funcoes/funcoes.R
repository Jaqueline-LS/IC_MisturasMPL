library("splines")
# Função que monta a matriz N
bspline<-function(x, ndx=15, bdeg=4)
{
  d<-bdeg-1
  nos1<-seq(min(x),max(x), length=ndx+2)
  delta<-nos1[2]-nos1[1]
  nos<-c(min(x)-delta*(d:1),nos1,max(x)+delta*(1:d))
  N<-splineDesign(knots = nos, x, bdeg, 0*x)
  return(N)
}


# Estimação de gama, beta e sigma 2 de um modelo semiparamétrico
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
  

  #y.hat<-N%*%a.gama #colocar fora
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
  return(list(sdEmp.b=sd.b, sdEmp.sigma2=sd.sigma2, ep.curva=ep.curva))
}

miEspNsp <-function(theta,y,X,N,K){
  # teta = [beta;sigma2;f;alpha];
  n=length(y)
  p=ncol(X)
  q0=ncol(K)
  beta0=theta[1:p]
  sigma2=as.numeric(theta[p+1])
  gama=theta[(p+2):(length(theta)-1)]
  alpha=as.numeric(theta[length(theta)])
  ################### I1(beta,gamma,sigma2) ######################
  mil=p+q0+1
  I1=matrix(0,mil,mil)
  I1[1:p,1:p]=1/sigma2*t(X)%*%X   # I_beta,beta
  I1[mil,mil]= n/(2*(sigma2^2))   # I_sigma2,sigma2
  I1[(p+1):(p+q0),(p+1):(p+q0)]=1/sigma2*t(N)%*%N+alpha*K   # I_gamma,gamma
  I1[1:p,(p+1):(p+q0)]= 1/sigma2*t(X)%*%N                # I_beta_gamma
  I1[(p+1):(p+q0),1:p]= 1/sigma2*t(N)%*%X                # I_gama_beta
  return(I1)
}


# Estimação dos gamas, beta e sigma 2 de um modelo aditivo parcialmente linear
EAPL<-function(y,X,Curvas,alfas,precision=1e-5)   
{
  n<-dim(X)[1]
  k<-length(Curvas)
  
  # valor inicial
  b.beta<- solve ( t(X)%*%X , t(X)%*%(y) )
  sigma2.hat<-as.numeric((t(y-(X%*%b.beta)) %*% (y-(X%*%b.beta))) / n)
  gammas<- vector("list", length = k)
  for(i in 1:k)
  {
    gammas[[i]]<-solve(t(Curvas[[i]]$N)%*%Curvas[[i]]$N + sigma2.hat*alfas[i]*Curvas[[i]]$K)%*%t(Curvas[[i]]$N)%*%(y-X%*%b.beta)
  }
  
  theta0<-matrix(c(b.beta,unlist(gammas),sigma2.hat))
  # convergência
  criterio<-precision+1
  while(criterio>precision)
  {
    vetor<-1:k
    aux<-function(j)
    {
      Curvas[[j]]$N%*%gammas[[j]]
    }
    N.gamma<-sapply(vetor,FUN=aux)
    sum.N.Gamma<-apply(N.gamma,1,sum)
    b.beta<-solve(t(X)%*%X)%*%t(X)%*%(y-sum.N.Gamma)
    sigma2.hat<-as.numeric(( t( y - (X%*%b.beta) - sum.N.Gamma) %*% ( y - (X%*%b.beta) - sum.N.Gamma) ) / n)
    if(k>1)
    {
      for(i in 1:k)
      {
        vetor<-1:k
        vetor<-vetor[-i]
        aux<-function(j)
        {
          Curvas[[j]]$N%*%gammas[[j]]
        }
        outras.curvas<-sapply(vetor,FUN=aux)
        sum.N.Gamma<-apply(outras.curvas,1,sum)
        gammas[[i]]<-solve(t(Curvas[[i]]$N)%*%Curvas[[i]]$N + sigma2.hat*alfas[i]*Curvas[[i]]$K)%*%t(Curvas[[i]]$N)%*%(y-X%*%b.beta-sum.N.Gamma)
      }
      outras.curvas<-sapply(vetor,FUN=aux)
      sum.N.Gamma<-apply(outras.curvas,1,sum)
      gammas[[i]]<-solve(t(Curvas[[i]]$N)%*%Curvas[[i]]$N + sigma2.hat*alfas[i]*Curvas[[i]]$K)%*%t(Curvas[[i]]$N)%*%(y-X%*%b.beta-sum.N.Gamma)
    }
    if(k==1)
    {
      gammas[[1]]<-solve(t(Curvas[[1]]$N)%*%Curvas[[1]]$N + sigma2.hat*alfas[1]*Curvas[[1]]$K)%*%t(Curvas[[1]]$N)%*%(y-X%*%b.beta) 
    }
    theta<-matrix(c(b.beta, unlist(gammas), sigma2.hat))
    criterio<-sqrt(sum((theta-theta0)^2)) #norma l2
    theta0<-theta
  }
  return(list(b=b.beta,a=gammas,sigma2=sigma2.hat))
}

MI.espAPL<-function(theta,y,X,Curvas,alfas)
{
  beta<-theta$b
  sigma2<-theta$sigma2
  p<-ncol(X)
  n<-dim(X)[1]
  q<-sapply(theta$a,FUN=length)
  MI<-matrix(0, nrow = p+sum(q)+1, ncol=p+sum(q)+1)
  
  #----------------------
  # qsi<-c(theta$b,unlist(theta$a))
  # K<-matrix(0,ncol=p+sum(q),nrow = p+sum(q))
  # aux<-0 # auxiliar
  # for(i in 1:length(q))
  # {
  #   K[(aux+p+1):(aux+q[i]+p),(aux+p+1):(aux+q[i]+p)]<- (alfas[i]/2)*Curvas[[i]]$K
  #   aux<-aux+q[i]
  # }
  # t(qsi)%*%K%*%qsi
  #-----------------------
  #bloco dos betas(linha)
  
  MI[1:p,1:p]<-t(X)%*%X
  posicao<-0
  for(i in 1:length(q))
  {
    MI[1:p,(p+posicao+1):(p+posicao+q[i])]<-t(X)%*%Curvas[[i]]$N
    posicao<-posicao+q[i]
  }
  
  # bloco (linha) dos gammas
  posicao<-0
  for(i in 1:length(q)) # i é coluna
  {
    
    MI[(p+posicao+1):(p+posicao+q[i]),1:p]<-t(Curvas[[i]]$N)%*%X
    posicao2<-0 # posicao nas colunas
    for(j in 1:length(q))
    {
      if(i==j)
      {
        MI[(p+posicao+1):(p+posicao+q[i]),(p+posicao2+1):(p+posicao2+q[j])]<-t(Curvas[[i]]$N)%*%Curvas[[i]]$N+alfas[i]*sigma2*Curvas[[i]]$K
        posicao2<-posicao2+q[j]
      }
      else
      {
         MI[(p+posicao+1):(p+posicao+q[i]),(p+posicao2+1):(p+posicao2+q[j])]<-t(Curvas[[i]]$N)%*%Curvas[[j]]$N
         posicao2<-posicao2+q[j]
      }
    }
    posicao<-posicao+q[i]
    
    
  }
  

  # linha do sigma 2
  MI[p+sum(q)+1,p+sum(q)+1]<-n/(2*sigma2)
  MI<-MI/sigma2
  return(MI)
  
}

BIC.APL<-function(alfas,y,X,Curvas)
{
  n<-nrow(X)
  p<-ncol(X)
  theta<-EAPL(y,X,Curvas,alfas=alfas)
  gammas<-theta$a
  vetor<-1:k
  aux<-function(j)
  {
    Curvas[[j]]$N%*%gammas[[j]]
  }
  N.gamma<-sapply(vetor,FUN=aux)
  sum.N.Gamma<-apply(N.gamma,1,sum)
  
  aux2<-function(j)
  {
    (alfas[[j]]/2)*t(gammas[[j]])%*%Curvas[[j]]$K%*%gammas[[j]]
  }
  penalty.j<-sapply(vetor,FUN=aux2)
  mu<-X%*%theta$b + sum.N.Gamma 
  
  lp <- (-n/2)*log(theta$sigma2) - ((1/(2*theta$sigma2)) * t(y-mu)%*%(y-mu)) - sum(penalty.j)
  
  aux3<-function(j)
  {
    sum(diag(Curvas[[j]]$N%*%solve(t(Curvas[[j]]$N)%*%Curvas[[j]]$N+alfas[j]*Curvas[[j]]$K*theta$sigma2)%*%t(Curvas[[j]]$N)))
  }
  
  # J<-cbind(X,Curvas[[1]]$N,Curvas[[2]]$N, Curvas[[3]]$N)
  # H<-J%*%(t(J))
  df.curvas<-sapply(vetor,FUN=aux3)
  num_param<-p+sum(df.curvas)+1
  bic<- -2*lp + num_param*log(n)
}

