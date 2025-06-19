rm(list = ls(all = TRUE))
library("Matrix")
library("mgcv")
library("dplyr")
library("kableExtra") 
source("funcoes/funcoes.R")
source("funcoes/funcoesMisturas.R")
source("funcoes/GerarAmostras_geral.R") 
cores<-c("#c2a5cf","#BAF3DE","#C9E69E","#FF9B95","#FFC29A")
#cores<-c("#BAF3DE","#C9E69E","#FF9B95","#FFC29A")

par<-par(pch=19)

g<-2

pii<-c(0.35, 0.65)

beta.verd<-list(c(5),c(-1))
sigma2.verd<-list(0.09,0.04)
c4<-function(t, d=1)d*exp(2*t - 1) + 1
c5<-function(t, d=1)d*sin(2*pi*t)

arg.grupos<-list(g1=list(beta=beta.verd[[1]],curva=list(f="c4",a=0,b=1,d=1),sigma2=sigma2.verd[[1]],intercepto=F),
                 g2=list(beta=beta.verd[[2]],curva=list(f="c5",a=0,b=1,d=1),sigma2=sigma2.verd[[2]],intercepto=F),
                 intervalos=list(c(0,1),c(0,1)))


n=300
alfas<-c(300,100)
# Tem que colocar a função para montar a N sem a restrição, a curva que eles utilizam não respeitam essa restrição
amostra<-rMisUniPLM(n, pii, p=length(beta.verd[[1]]), arg=arg.grupos)
theta<-try(EMisUniPLM(g=2,alfas=alfas,y=amostra$y, t=amostra$t, X=amostra$X),TRUE)

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


plot(amostra$X[,1],amostra$y, cex=0.5,col=cores[grupo],
       main=paste0("y vs. x - ",tipo),ylab="y", xlab="x_i")
  abline(a=betas[[1]][1],b=betas[[1]][1], lwd=2, lty=2)
  abline(a=betas[[2]][1],b=betas[[2]][1], lwd=2, lty=2)
  abline(a=beta.verd[[1]][1],b=beta.verd[[1]][1], lwd=3)
  abline(a=beta.verd[[2]][1],b=beta.verd[[2]][1], lwd=3)
  legend("topright",legend = c("estimada","verdadeira"),lty = c(2,1), lwd = c(2,3), cex=0.7)


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
