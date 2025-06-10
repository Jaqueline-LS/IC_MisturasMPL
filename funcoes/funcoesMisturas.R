library("e1071")     # para SVM


# Curvas que respeitam a restrição
c1<- function(x,d=2) d*cos(2*pi*x)
c2<- function(x,d=2) d*sin(2*pi*x)*exp(-0.5*x^2) #curva
c3<- function(x,d=3) d*sin(3*pi*x)*exp(0.5*x^(4)) #curva


# curve(c1,-1,1, ylim=c(-5,5))
# curve(c2,-1,1, add=T)
# curve(c3,-1,1, add=T)



#-------------------------

rPLM<-function(n,p,beta,curva,sigma2, intercepto=FALSE, intervalos) #curva é uma lista "c1", a, b                                                                           #intervalos é list, cada objeto é um vetor do intervalo de uma variável explicativa 
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
  y<-X%*%beta + curva.f(t,d=d)  + erro
  
  
  modelo<-list(y=y,X=X,t=t)
  return(modelo)
}


rMisUniPLM<- function(n, pii, p, arg) 
{
  # n: tamanho da amostra gerada
  # pii: proporções da mistura
  # arg: Uma lista com os argumentos de cada grupo,intercepto,betas, curva, intervalo da curva
  y<-vector(mode = "numeric", length = n)
  X<-matrix(NA, nrow = n, ncol=p)
  t<-vector(mode="numeric", length = n)
  z<-rmultinom(n,1,pii)
  g <- length(pii)
  for(i in 1:n) 
  {
    if(g==1)
    {
      j<-1
    }else{
      j<-which(z[,i]==T)
    }
    modelo<-rPLM(n=1,p=p,beta=arg[[j]]$beta,curva = arg[[j]]$curva,
                 sigma2=arg[[j]]$sigma2, intercepto = arg[[j]]$intercepto,
                 intervalos = arg$intervalos) 
    # intervalos 
    # # é o intervalo da uniforme que será gerada para X, nas simulações deixei todos iguais a (0,1)
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
# 
# #-----------------Uma amostra de teste do Cenario 2-----------------
# pii<-c(0.35, 0.65)
# beta.verd<-list(c(3,4,6),c(2,3,-5))
# sigma2.verd<-list(2,1)
# arg.grupos<-list(g1=list(beta=beta.verd[[1]],curva=list(f="c1",a=-1,b=1,d=2),sigma2=sigma2.verd[[1]],intercepto=T),
#                  g2=list(beta=beta.verd[[2]],curva=list(f="c2",a=-1,b=1,d=4),sigma2=sigma2.verd[[2]],intercepto=T),
#                  intervalos=list(c(0,1),c(0,1)))
# n=300
# alfas<-c(0.1,0.1)
# amostra<-rMisUniPLM(n, pii, p=length(beta.verd[[1]]), arg=arg.grupos)
# g=2
# alfas=alfas
# y=amostra$y
# t=amostra$t
# X=amostra$X
# clu=amostra$clu
# iter.max=100
# n.start=50
# #-----------------------------------------------


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
  
  
  #------- Construção das matrizes N e K para as curvas-----
  d<-data.frame(t=t)
  ZZ<-smoothCon(s(t,bs="ps",m=c(2,3)),data=d,knots=NULL,absorb.cons=T) 
  # absorb.cons = T ( restrição de idenficabilidade )
  # bs="ps",m=c(2,3) P-splines, B-splines cúbicos com penalização de ordem 2
  N<-ZZ[[1]]$X
  K<-(ZZ[[1]]$S)[[1]]
  
  
  # Outra opção para montar N e K, não tem restrição de identificabilidade

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

    sigma2.0[j]<-as.numeric(t(y-X%*%beta0[[j]])%*%H.j%*%(y-X%*%beta0[[j]]))/sum(Z[[j]])

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

      sigma2[j] <- as.numeric(t(y-X%*%beta[[j]]-N%*%gamma[[j]]) %*%H.j%*% (y-X%*%beta[[j]]-N%*%gamma[[j]]) ) / sum(z[,j])

      mu[[j]] <- as.vector(X%*%beta[[j]]+N%*%gamma[[j]])
      
    }


    l<-log.vero.MisUniPLM(y,p,g,mu=mu,sigma2=sigma2)
    criterio<-abs(l/l0-1)

    l0<-l

    theta<-list(b=beta,a=gamma,sd2=sigma2,p=p)

    theta0<-theta
    
    mu0<-mu

    cont<-cont+1

  }
  return(list(theta=theta,z=z,N=N, K=K))
}

# Função com SVM
EMisUniPLM.SVM<-function(g,alfas,y,t,X,iter.max=100, n.start=50, clu)
# clu são os grupos verdadeiros
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
  
  
  #--------------------------- SVM----------------------------
  rotulos_reordenados <- as.factor(rotulos_reordenados)
  svm_model <- svm(rotulos_reordenados ~ y, kernel = "polynomial",degree=3, probability = TRUE, fitted=T)
  label_svm <- svm_model$fitted
  p0<-sapply(1:g,FUN=function(x)mean(label_svm==x))
  
  # Agrupamento verdadeiro
  grupo<-apply(clu,2,FUN=function(x)which(as.logical(x))) # numero do grupo de cada observação
  
  
  # Acuracia com o k-means
  
  tabela<-table(grupo, rotulos_reordenados)
  acuracia<-sum(diag(tabela))/sum(tabela)
  acuracia
  
  # Acuracia depois o SVM
  
  tabela2<-table(grupo, label_svm)
  acuracia2<-sum(diag(tabela2))/sum(tabela2)
  acuracia2
  
  # # Mudar os rotulos de acordo com a proporção menor não funciona depois no SVM
  # # O SVM que separou mal. O que ele definiu como grupo 1 é realmente o grupo 1
  # o.g<-order(p0, decreasing = F)
  # p0<-p0[o.g]
  # label_svm.reordenado<-sapply(label_svm,FUN=function(y)which(o.g==y))
  # # Acurácia depois do reordenamento
  # tabela2<-table(grupo, label_svm.reordenado)
  # acuracia2<-sum(diag(tabela2))/sum(tabela2)
  # acuracia2
  
  
  rotulos_reordenados<-label_svm
  
  
  
  
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
  
  return(list(theta=theta,z=z,N=N, K=K, acuracia=acuracia, acuracia2=acuracia2))
  #return(list(theta=theta,z=z,N=N, K=K))
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



#-------------------------
verificarGeração<-function(amostra, arg.grupos, tipo)
{
  g<-length(arg.grupos)-1
  y<-amostra$y
  t<-amostra$t
  X<-amostra$X
  p<-ncol(X)
  
  
  hist(y,freq=F, breaks=15, main=tipo)
  
  grupo<- apply(amostra$clu,2,FUN=function(x)which(as.logical(x))) # numero do grupo de cada observação
  
  #plot(amostra$t, amostra$y, pch=19, cex=0.5,col=cores[grupo])
  # 
  # for (i in 1:p)
  # {
  #   plot(amostra$X[,2],amostra$y, cex=0.5,col=cores[grupo])
  # }
  # 
  # apply(amostra$clu,1, mean) # proporções geradas 
  
  # #--------Plot da variável t pela y menos a parte linear-------
  # g.i<-vector("list",g)
  # 
  # for(i in 1:g)
  # {
  #   g.i[[i]]<-which(as.logical(amostra$clu[i,])) # indice das obs. do grupo i
  # }
  # matrizX<-lapply(g.i,FUN=function(x)X[x,])
  # X.grand<-as.matrix(bdiag(matrizX)) # matriz bloco diagonal
  # 
  # beta.grand<-unlist(beta.verd)
  # partelin<-X.grand%*%beta.grand # Efeito das variáveis lineares em Y
  # 
  # ordem<-unlist(g.i) #indice em que as observações estão em X.g
  # 
  # plot(amostra$t[ordem], amostra$y[ordem]-partelin, pch=19, cex=0.5,col=cores[grupo[ordem]])
  # c<-vector("list",g)
  # for(i in 1:g)
  # {
  #   c[[i]]<-get(arg.grupos[[i]]$curva$f)
  #   d<-arg.grupos[[i]]$curva$d
  #   a<-arg.grupos[[i]]$curva$a
  #   b<-arg.grupos[[i]]$curva$b
  #   curve(c[[i]](x,d),a,b,add=T, lwd=3, lty=2)
  # }
}

VerificarEstimacao<-function(amostra, arg.grupos, theta)
{
  
  N<-theta$N
  z<-theta$z
  g<-length(arg.grupos)-1
  y<-amostra$y
  t<-amostra$t
  X<-amostra$X
  p<-ncol(X)
  
  
  grupo<-apply(amostra$clu,2,FUN=function(x)which(as.logical(x))) # numero do grupo de cada observação
  
  # Pega os parâmetros estimados
  betas=theta$theta$b
  gamma=theta$theta$a
  
  for (i in 2:p) # y vs. x
  {
    plot(amostra$X[,i],amostra$y, cex=0.5,col=cores[grupo],
         main=paste0("y vs. x - ",tipo),ylab="y", xlab="x_i")
    abline(a=betas[[1]][1],b=betas[[1]][i], lwd=2, lty=2)
    abline(a=betas[[2]][1],b=betas[[2]][i], lwd=2, lty=2)
    abline(a=beta.verd[[1]][1],b=beta.verd[[1]][i], lwd=3)
    abline(a=beta.verd[[2]][1],b=beta.verd[[2]][i], lwd=3)
    legend("topright",legend = c("estimada","verdadeira"),lty = c(2,1), lwd = c(2,3), cex=0.7)
  }
  
  
  plot(amostra$t, amostra$y, pch=19, cex=0.5,col=cores[grupo], 
       main=paste0("y vs. t - ",tipo))
  for(i in 1:g)
  {
    y_hat<-as.numeric(N%*%gamma[[i]])
    dados<-data.frame(t,y_hat)
    dados<-dados[order(dados$t, decreasing=FALSE),]
    lines(dados$t,dados$y_hat, col=cores[i], lwd=3)
    
    curva<-get(arg.grupos[[i]]$curva$f)
    d<-arg.grupos[[i]]$curva$d
    a<-arg.grupos[[i]]$curva$a
    b<-arg.grupos[[i]]$curva$b
    curve(curva(x,d),a,b,add=T, lwd=3, lty=2)
    
  }
  legend("topright",legend = c("verdadeira","estimada"),lty = c(2,1), lwd = c(2,3), cex=0.4)
  
  
  
  #--------Plot da variável t pela y menos a parte linear-----------------
  # ----------------Ou plot da X menos a curva e as demais lineares-------
  # Residuo não paramétrico
  g.i<-vector("list",g)
  for(i in 1:g)
  {
    g.i[[i]]<-which(as.logical(amostra$clu[i,])) # indice das obs. do grupo i
  }

  matrizX<-lapply(g.i,FUN=function(x)X[x,]) # Pega as observações que pertecem ao grupo
  X.grand<-as.matrix(bdiag(matrizX)) # matriz bloco diagonal

  beta.grand<-unlist(betas) # X%*%beta.chapeu
  partelin<-X.grand%*%beta.grand # Efeito das variáveis lineares em Y

  matrizN<-lapply(g.i,FUN=function(x)N[x,]) # Pega as observações que pertecem ao grupo
  N.grand<-as.matrix(bdiag(matrizN)) # matriz bloco diagonal
  a.grand<-unlist(gamma)
  curva<-N.grand%*%a.grand # Efeito das variáveis lineares em Y

  ordem<-unlist(g.i) #indice em que as observações estão em X.g e N.g

  # for (i in 2:p) # y-n_i*gamma vs. x
  # {
  #   plot(amostra$X[ordem,i],amostra$y[ordem]-curva, cex=0.5,col=cores[grupo[ordem]])
  #   abline(a=b[[1]][1],b=b[[1]][i], lwd=3)
  #   abline(a=b[[2]][1],b=b[[2]][i], lwd=3)
  #   abline(a=4,b=4, lwd=3, lty=2)
  #   abline(a=b[[2]][1],b=b[[2]][i], lwd=3)
  # }
  #


  plot(amostra$t[ordem], amostra$y[ordem]-partelin, pch=19, cex=0.5,col=cores[grupo[ordem]])
  for(i in 1:g)
  {
    y_hat<-as.numeric(N%*%gamma[[i]])
    dados<-data.frame(t,y_hat)
    dados<-dados[order(dados$t, decreasing=FALSE),]
    lines(dados$t,dados$y_hat, col=cores[i], lwd=3)

    curva<-get(arg.grupos[[i]]$curva$f)
    d<-arg.grupos[[i]]$curva$d
    a<-arg.grupos[[i]]$curva$a
    b<-arg.grupos[[i]]$curva$b
    curve(curva(x,d),a,b,add=T, lwd=3, lty=2)

  }
  legend("topright",legend = c("estimada","verdadeira"),lty = c(2,1), lwd = c(2,3))
  
}



