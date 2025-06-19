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



# Define o número de réplicas e o número de grupos

M<-500
g<-1


# Define os parametros do cenario
# # Cenario 1
pii<-c(0.35, 0.65)
pii<-1
beta.verd<-list(c(8,4,6),c(2,3,-5))
beta.verd<-list(c(8,4,6))


sigma2.verd<-1


#arg.grupos<-list(g1=list(beta=beta.verd[[1]],curva=list(f="c1",a=-1,b=1,d=2),sigma2=sigma2.verd[[1]],intercepto=T),
                # g2=list(beta=beta.verd[[2]],curva=list(f="c2",a=-1,b=1,d=4),sigma2=sigma2.verd[[2]],intercepto=T),
                # intervalos=list(c(0,1),c(0,1)))

arg.grupos<-list(g1=list(beta=beta.verd[[1]],curva=list(f="c1",a=-1,b=1,d=2),sigma2=sigma2.verd[[1]],intercepto=T),
                 intervalos=list(c(0,1),c(0,1)))


# Para gerar uma amostra de teste com o agrupamento do k-mean

n=1000
alfas<-c(0.1)
amostra<-rMisUniPLM(n, pii, p=length(beta.verd[[1]]), arg=arg.grupos)
t<-amostra$t
y=amostra$y
X=amostra$X
d<-data.frame(t=t)
ZZ<-smoothCon(s(t,bs="ps",m=c(2,3)),data=d,knots=NULL,absorb.cons=T) 
# absorb.cons = T ( restrição de idenficabilidade )
# bs="ps",m=c(2,3) P-splines, B-splines cúbicos com penalização de ordem 2
N<-ZZ[[1]]$X
K<-(ZZ[[1]]$S)[[1]]
theta<-estimar.semiParametrico(y=amostra$y,X=X,N=N,K=K, alfa=alfas)

a=theta$a
b=theta$b
sigma2<-theta$sigma2


par(mfrow=c(1,1))

#-------------Envelope simulado com o resíduo quantíco semi parametrico------------
# Resíduo quantilico
residuo.quantilico<-qnorm(pnorm(y,mean=N%*%a + X%*%b,sd=sqrt(sigma2)))
y.sim<-matrix(0,n,100)
epsilon<-matrix(0,n,100)
e1<-numeric(n)
e2<-numeric(n)

for(i in 1:100)
{
  y.sim[,i]<-rnorm(n,mean=N%*%a + X%*%b,sd=sqrt(sigma2))
  theta<-estimar.semiParametrico(y=y.sim[,i],X=X,N=N,K=K, alfa=alfas)
  a.sim=theta$a
  b.sim=theta$b
  sigma2.sim<-theta$sigma2
  epsilon[,i]<-qnorm(pnorm(y.sim[,i],mean=N%*%a.sim + X%*%b.sim,sd=sqrt(sigma2.sim)))
  epsilon[,i]<-sort(epsilon[,i])
}
for(i in 1:n){
  eo <- sort(epsilon[i,])
  e1[i] <- (eo[2]+eo[3])/2
  e2[i] <- (eo[97]+eo[98])/2 }

med <- apply(epsilon,1,mean)
faixa <- range(residuo.quantilico,e1,e2)
#
par(pty="s")
pp=qqnorm(residuo.quantilico,xlab="Percentis da N(0,1)",
          ylab="Residuo Quantilico", ylim=faixa, pch=16)
par(new=T)
qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
par(new=T)
qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1)
par(new=T)
qqnorm(med,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=2)

alfa<-0.05
ind <- numeric(n)
ind <- sort(residuo.quantilico)<sort(e1) | sort(residuo.quantilico)>sort(e2)
prop.fora<-sum(ind)/n

prop.fora
