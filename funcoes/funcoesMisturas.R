# Curvas que respeitam a restrição

c1<- function(x,d=2) d*cos(2*pi*x)
c2<- function(x,d=2) d*sin(2*pi*x)*exp(-0.5*x^2) #curva
c3<- function(x,d=3) d*sin(3*pi*x)*exp(0.5*x^(4)) #curva


#-------- Curvas-------------
# c1<- function(x,d=1) d*exp(sin(pi*x))
# c2<- function(x,d=1) d*cos(pi*x)

# #-------- Curvas-------------
# c1<- function(x,d=0) 2*cos(pi*x)+d
# c2<- function(x,d=0) 3*cos(pi*x)+d

# curve(c1,-1,1, ylim=c(-5,5))
# curve(c2,-1,1, add=T)
# curve(c3,-1,1, add=T)
# 


#--------------------
amostraMU<-function(n,pii,arg.grupos)
{
  z<-rmultinom(n,1,pii)
  indice<-sapply(1:length(pii),FUN=function(i)which(z[i,]==1))
  y<-numeric(n)
  t<-numeric(n)
  beta<-arg.grupos[[1]]$beta
  p<-length(beta)
  X<-matrix(0,n,p)
  g<-length(indice)

  for(i in 1:g)
  {
    intercepto = arg.grupos[[i]]$intercepto
    l<-indice[[i]]
    ni<-length(l)
    beta<-arg.grupos[[i]]$beta
    p<-length(beta)
    
    if(!intercepto) # sem intercepto (itercepto = F)
    {
      X[l,]<-sapply(1:p, FUN=function(p)X[l,p]<-runif(ni,min=0, max = 1))
    }else{
      X[l,1]<-rep(1,ni)
      vec<-1:p
      vec<-vec[-1]
      X[l,vec]<-sapply(vec, FUN=function(k)X[l,k]<-runif(ni,min=0, max = 1))
    }
    curva<-arg.grupos[[i]]$curva
    a<-curva$a
    b<-curva$b
    t[l]<-runif(ni,a, b)
  }
  return(list(X=X,t=t,z=z))
}

rMisUniPLM.mu<-function(X,t,z, arg.grupos)
{
  indice<-sapply(1:length(pii),FUN=function(i)which(z[i,]==1))
  y<-numeric(n)
  erro<-numeric(n)
  g<-length(indice)
  
  for(i in 1:g)
  {
    l<-indice[[i]] # obs. do grupo i
    n<-length(l)
    beta<-arg.grupos[[i]]$beta

    curva<-arg.grupos[[i]]$curva
    f<-curva$f
    sigma2<-arg.grupos[[i]]$sigma2
    d<-curva$d
    erro[l]<-rnorm(n, mean = 0, sd=sqrt(sigma2))
    curva.f<-get(f)
    #constrained.f<-curva.f(t[l],d=d)-sum(curva.f(t[l],d=d))
    # y[l]<-X[l,]%*%beta + constrained.f + erro[l]
    y[l]<-X[l,]%*%beta + curva.f(t[l],d=d) + erro[l]
    
    
  }
  dados<-list(y=y,X=X,t=t,clu=z)
  return(dados)
  
}



#-------------------------
verificarGeração<-function(amostra, arg.grupos)
{
  g<-length(arg.grupos)-1
  y<-amostra$y
  t<-amostra$t
  X<-amostra$X
  p<-ncol(X)
  
  
  hist(y,freq=F, breaks=15)
  
  grupo<- apply(amostra$clu,2,FUN=function(x)which(as.logical(x))) # numero do grupo de cada observação
  
  #plot(amostra$t, amostra$y, pch=19, cex=0.5,col=cores[grupo])
  
  # for (i in 1:p)
  # {
  #   plot(amostra$X[,2],amostra$y, cex=0.5,col=cores[grupo])
  # }
  
  apply(amostra$clu,1, mean) # proporções geradas 
  
  #--------Plot da variável t pela y menos a parte linear-------
  g.i<-vector("list",g)
  
  for(i in 1:g)
  {
    g.i[[i]]<-which(as.logical(amostra$clu[i,])) # indice das obs. do grupo i
  }
  matrizX<-lapply(g.i,FUN=function(x)X[x,])
  X.grand<-as.matrix(bdiag(matrizX)) # matriz bloco diagonal
  
  beta.grand<-unlist(beta.verd)
  partelin<-X.grand%*%beta.grand # Efeito das variáveis lineares em Y
  
  ordem<-unlist(g.i) #indice em que as observações estão em X.g
  
  plot(amostra$t[ordem], amostra$y[ordem]-partelin, pch=19, cex=0.5,col=cores[grupo[ordem]])
  c<-vector("list",g)
  for(i in 1:g)
  {
    c[[i]]<-get(arg.grupos[[i]]$curva$f)
    d<-arg.grupos[[i]]$curva$d
    a<-arg.grupos[[i]]$curva$a
    b<-arg.grupos[[i]]$curva$b
    curve(c[[i]](x,d),a,b,add=T, lwd=3, lty=2)
  }
}

#-------------------------

rPLM<-function(n,p,beta,curva,sigma2, intercepto=FALSE, intervalos) #curva é list "c1", a, b                                                                           #intervalos é list, cada objeto é um vetor do intervalo de uma variável explicativa 
{
  a<-curva$a
  b<-curva$b
  f<-curva$f
  d<-curva$d
  X<-matrix(0,nrow=n, ncol=p)
  if(!intercepto) # sem intercepto (itercepto = F)
  {
    X<-sapply(1:dim(X)[2], FUN=function(p)X[,p]<-runif(n,min=intervalos[[p]][1], max = intervalos[[p]][2]))
  }else{
    X[,1]<-rep(1,n)
    vec<-(1:dim(X)[2])
    for(i in vec[-1])
    {
      X[,i]<-runif(n,min=intervalos[[i-1]][1], max = intervalos[[i-1]][2])
    }
   
  }
  
  t<-runif(n,a, b)
  curva.f<-get(f)
  erro<-rnorm(n, 0, sd=sqrt(sigma2))
  #constrained.f<-curva.f(t,d=d)-sum(curva.f(t,d=d))
  #y<-X%*%beta + constrained.f + erro
  y<-X%*%beta + curva.f(t,d=d)  + erro
  
  
  modelo<-list(y=y,X=X,t=t)
  return(modelo)
}


rMisUniPLM<- function(n, pii, p, arg) 
{
  # n: tamanho da amostra gerada
  # pii: proporções da mistura
  # family: Cada grupo é um modelo parcialmente linear (PLM)
  # arg: Uma lista com os argumentos de cada grupo,intercepto,betas, curva, intervalo da curva
  # As curvas disponíveis são c1, c2.
  y<-vector(mode = "numeric", length = n)
  X<-matrix(NA, nrow = n, ncol=p)
  t<-vector(mode="numeric", length = n)
  z<-rmultinom(n,1,pii)
  g <- length(pii)
  for(i in 1:n) 
  {
    j<-which(z[,i]==T)
    modelo<-rPLM(n=1,p=p,beta=arg[[j]]$beta,curva = arg[[j]]$curva,
                 sigma2=arg[[j]]$sigma2, intercepto = arg[[j]]$intercepto,
                 intervalos = arg$intervalos)
    y[i] <-modelo$y 
    X[i,]<-modelo$X
    t[i]<-modelo$t
  }
  dados<-list(y=y,X=X,t=t,clu=z)
  return(dados)
}



log.vero.MisUniPLM<-function(y,p,g,mu,sigma2)
{
  p[g]<-1-sum(p)
  s<-sapply(1:g,FUN=function(j)p[j]*dnorm(y, mean = mu[[j]], sd=sqrt(sigma2[j])))
  l<-sum(log(apply(s, 1, sum)))
  return(l)
}

EMisUniPLM<-function(g,alfas,y,t,X,iter.max=100, n.start=50)
{
  # Agrupamento inicial K-means
  n<-length(y)
  km<-kmeans(y,centers=g, iter.max=iter.max, nstart=n.start, algorithm = c("Hartigan-Wong"))
  p0<-km$size/sum(km$size)
  
  # Mudar o label, o grupo que possui o maior valor de centro é o primeiro
  # Segue em ordem decrescente
  o.g<-order(km$centers, decreasing = T)
  p0<-p0[o.g]
  rotulos_reordenados<-sapply(km$cluster,FUN=function(y)which(o.g==y))

  
  #------- Construção da matriz N que irá definir as curvas-----
  d<-data.frame(t=t)
  ZZ<-smoothCon(s(t,bs="ps",m=c(2,3)),data=d,knots=NULL,absorb.cons=T) 
  # absorb.cons = T ( restrição de idenficabilidade )
  # bs="ps",m=c(2,3) P-splines, B-splines cúbicos com penalização de ordem 2
  N<-ZZ[[1]]$X
  K<-(ZZ[[1]]$S)[[1]]
  
  
  # Outra opção para montar a N, não tem restrição de identificabilidade

  # N<-bspline(t, ndx=g*15)
  # q<-ncol(N)
  # D<-diag(q)
  # D<-diff(D, differences = 2)
  # K<-t(D)%*%D
  
  #----------------------- Valores iniciais-----------------------
  # Valor inicial dos parâmetros de acordo com o agrupamento inicial do K-means
  beta0<-vector("list",g)
  gamma0<-vector("list",g)
  sigma2.0<-vector("numeric",g)
  mu0<-vector("list",g)

  Z<-list(g)
  for(j in 1:g)
  {
    # matriz onde a diagonal diz se a observação pertence ao grupo j
    
    Z[[j]]<-diag(as.numeric(rotulos_reordenados==j)) 
    #Z[[j]]<-diag(as.numeric(km$cluster==j)) 
  
  }
  
  
  for(j in 1:g)
  {
    H.j<-Z[[j]]
  
    beta0[[j]]<-solve(t(X)%*%H.j%*%X)%*%t(X)%*%H.j%*%y

    sigma2.0[[j]]<-as.numeric(t(y-X%*%beta0[[j]])%*%H.j%*%(y-X%*%beta0[[j]]))/sum(Z[[j]])

    gamma0[[j]]<-solve(t(N)%*%H.j%*%N+alfas[j]*sigma2.0[[j]]*K)%*%t(N)%*%H.j%*%(y-X%*%beta0[[j]])

    sigma2.0[j]<-as.numeric(t(y-X%*%beta0[[j]]-N%*%gamma0[[j]])%*%H.j%*%(y-X%*%beta0[[j]]-N%*%gamma0[[j]]))/sum(Z[[j]])

    gamma0[[j]]<-solve(t(N)%*%H.j%*%N+alfas[j]*sigma2.0[j]*K)%*%t(N)%*%H.j%*%(y-X%*%beta0[[j]])

    mu0[[j]]<-as.vector(X%*%beta0[[j]]+N%*%gamma0[[j]])

  }
  
  theta0<-list(b=beta0,a=gamma0,sd2=sigma2.0,p=p0[-g])
  l0<-log.vero.MisUniPLM(y,p0[-g],g,mu=mu0,sigma2=sigma2.0)

  criterio<-1E5
  cont<-0

  p<-vector("numeric",g-1)
  beta<-vector("list",g)
  gamma<-vector("list",g)
  sigma2<-vector("numeric",g)
  mu<-vector("list",g)

  
  while(criterio>=1E-5 && cont<5000)
  {
    beta0<-theta0$b
    gamma0<-theta0$a
    sigma2.0<-theta0$sd2
    p0<-theta0$p

    #----------- Etapa E -----------

    z<-matrix(0,nrow = n, ncol=g)
    for(j in 1:g)
    {
      # proporção do grupo j vezes o valor da densidade para cada observação
      # esse vetor preenche a coluna da matriz z
      ifelse(j==g,z[,j]<-(1-sum(p0))*dnorm(y,mean = mu0[[j]], sd = sqrt(sigma2.0[j])),z[,j]<-p0[j]*dnorm(y,mean = mu0[[j]], sd = sqrt(sigma2.0[j])))
    }
    sum.pj.dnorm<-apply(z, 1, FUN=sum) # um vetor, onde cada valor é soma das colunas de z
    
    for(j in 1:g)
    {
      # z_ij chapéu
      z[,j]<-z[,j]/sum.pj.dnorm
    }

    #--------------Etapa M-------------------
    # Reescrevendo o algoritmo para estimar apenas g-1 proporções a outra é obtida pelo complementar
    if(g==2)
    {
      p<-mean(z[,1])
    }else{
      p<-apply(z[,-g], 2, mean) # media de cada coluna de z
    }
   

    for(j in 1:g)
    {
      H.j<-diag(z[,j]) # matriz diagonal com os valores de z_ij chapéu

      beta[[j]] <- solve(t(X)%*%H.j%*%X) %*%t(X)%*%H.j%*% (y-N%*%gamma0[[j]])

      gamma[[j]] <- solve(t(N)%*%H.j%*%N + alfas[j]*sigma2.0[[j]]*K ) %*%t(N)%*%H.j%*% (y-X%*%beta[[j]])

      sigma2[[j]] <- as.numeric(t(y-X%*%beta[[j]]-N%*%gamma[[j]]) %*%H.j%*% (y-X%*%beta[[j]]-N%*%gamma[[j]]) ) / sum(z[,j])

      mu[[j]] <- as.vector(X%*%beta[[j]]+N%*%gamma[[j]])
      
    }


    l<-log.vero.MisUniPLM(y,p,g,mu=mu,sigma2=sigma2)
    criterio<-abs(l/l0-1)

    l0<-l

    theta<-list(b=beta,a=gamma,sd2=sigma2,p=p)

    theta0<-theta

    cont<-cont+1

  }

  return(list(theta=theta,z=z,N=N, K=K))

}




MI.MisUniPLM<-function(y,X,alfas,theta) # aproximação (pela soma do produto das scores de uma observação)
{
  N<-theta$N
  K<-theta$K
  z<-theta$z
  
  params<-theta$theta
  g<-length(params$p)+1
  n<-nrow(X)
  q<-ncol(N)
  p<-sapply(params$b,FUN=length)
  MI<-matrix(0, nrow = (g-1)+sum(p)+(g*q)+g, ncol=(g-1)+sum(p)+(g*q)+g)
  
  S.p<-vector("list",g-1)
  S.b<-vector("list",g)
  S.a<-vector("list",g)
  S.sigma<-vector("list",g)
  
  for(i in 1:n)
  {
    for(j in 1:g)
    {
      alfa<-alfas[j]
      sigma2<-params$sd2[[j]]
      beta<-params$b[[j]]
      gama<-params$a[[j]]
      if(j!=g)
      {
        pii<-params$p[[j]]
      }
      p.g<-1-sum(params$p)
      z.j<-z[,j]
      z.g<-z[,g]
      sigma<-sqrt(sigma2)
      
      # S_p, p1,...,p_{g-1}
      if(j!=g)
      {
        S.p[[j]]<- z.j[i]/pii - z.g[i]/p.g
      }
      S.b[[j]] <- (1/sigma2) * z[i,j] %*% (y[i] -t(X[i,])%*%beta-t(N[i,])%*%gama) %*% X[i,]
    
      S.a[[j]] <- (z[i,j]/sigma2) * (y[i] -t(X[i,])%*%beta-t(N[i,])%*%gama) %*% N[i,] - t((alfa/n)*(K%*%gama))
    
      S.sigma[[j]] <- (-z[i,j]/(2*sigma2)) + ((z[i,j]/(2*(sigma2^2)))*(y[i] - t(X[i,])%*%beta-t(N[i,])%*%gama)^2)
      
    }
    S.i<-c(unlist(S.p),unlist(S.b), unlist(S.a), unlist(S.sigma))
    S<-S.i%*%t(S.i)
    MI<-MI+S
  }
  return(MI)
  
}

BICMisUniPLM<-function(alfas,g,y,t,X1)
{
  X<-X1
  theta<-EMisUniPLM(g=g,alfas=alfas,y=y,t=t,X=X,iter.max=100, n.start=50)
  N<-theta$N
  K<-theta$K
  z<-theta$z
  params<-theta$theta
  pii<-params$p
  n<-nrow(X)
  q<-ncol(N)
  p<-ncol(X)
  beta<-params$b
  sigma2<-as.numeric(params$sd2)
  gamma<-params$a
  mu<-vector("list",g)
  dfN<-numeric(g)
  
  for(j in 1:g)
  {
    H.j<-diag(z[,j])
    mu[[j]] <- as.vector(X%*%beta[[j]]+N%*%gamma[[j]])
    
    # Degrees of freedom curva j
    
    auxH1<-solve(t(N)%*%H.j%*%N + alfas[j]*sigma2[j]*K)%*%t(N)%*%H.j
    auxH2<-X%*%solve(t(X)%*%H.j%*%X)%*%t(X)%*%H.j
    auxH3<-(diag(q)-auxH1%*%auxH2%*%N)
    auxH3<-(auxH3+t(auxH3))/2
    auxH3<-solve(auxH3)
    H<-N%*%auxH3%*%(auxH1-auxH1%*%auxH2)
    dfN[j]<-sum(diag(H))
    
  } 
  
  
  l<-log.vero.MisUniPLM(y,pii,g,mu=mu,sigma2=sigma2)
  
  dfN<-sum(dfN)
  
  BIC<-(-2*l)+((g-1)*p+g+g+dfN)*log(n)
  return(BIC)
}




