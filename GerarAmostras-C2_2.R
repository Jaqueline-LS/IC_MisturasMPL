rm(list = ls(all = TRUE))
library("Matrix")
library("mgcv")
library("dplyr")
library("kableExtra")
source("funcoes/funcoes.R")
source("funcoes/funcoesMisturas.R")
cores<-c("#BAF3DE","#FF9B95","#C9E69E","#FFC29A")
par<-par(pch=19)


M<-500
g<-2

#Cenario 2
# pii<-c(0.35, 0.65)
# 
# beta.verd<-list(c(3,4,6),c(2,3,-5))
# nome.amostra<-"Amostras/c2_2_"
# nome.plot<-"c2_2_"

# # # Cenario 2.1
# 
# pii<-c(0.47, 0.53)
# 
# beta.verd<-list(c(3,4,6),c(2,3,-5))
# nome.amostra<-"Amostras/c2.1_2_"
# nome.plot<-"c2.1_2_"


# #Cenario 3.1
# beta.verd<-list(c(4,4,6),c(2,3,-5))
# nome.amostra<-"Amostras/c2.2_2_"
# nome.plot<-"c2.2_2_"


# Cenario 3.2
pii<-c(0.35,0.65)

beta.verd<-list(c(4,4,6),c(2,3,-5))
nome.amostra<-"Amostras/c3.2_2_"
nome.plot<-"c3.2_2_"


#------------------------------
sigma2.verd<-list(1,1)


arg.grupos<-list(g1=list(beta=beta.verd[[1]],curva=list(f="c1",a=-1,b=1,d=2),sigma2=sigma2.verd[[1]],intercepto=T),
                 g2=list(beta=beta.verd[[2]],curva=list(f="c2",a=-1,b=1,d=4),sigma2=sigma2.verd[[2]],intercepto=T),
                 intervalos=list(c(0,1),c(0,1)))

n=3000

amostra<-rMisUniPLM(n, pii, p=length(beta.verd[[1]]), arg=arg.grupos) 


verificarGeração(amostra,arg.grupos)

#---------------------------------------------------
alfas<-c(0.1,0.1)

sizes<-c(300,500,1000,2000)


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
    
    #erro.padrao<-list(sd.p=sd.p,sd.b=sd.b, sd.a=sd.a,sd.sigma2=sd.sigma2, ginv=aux)
    
    erro.padrao<-list(sd.p=sd.p,sd.b=sd.b, sd.sigma2=sd.sigma2)
    
    teptetam[[i]]<-list(t=amostra$t, theta, erro.padrao, y=amostra$y,X=amostra$X,alfas=alfas,theta)
    
  }
  saveRDS(teptetam, file=paste0(nome.amostra,n,".RDS"))
}
sapply(sizes,FUN = gera.amostras)

source("graficos_MSE.R")

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
  
  legend("topright", legend=c("grupo 1","grupo 2"), bty = "n",lwd = 3, col=cores[1:2], cex=0.5 )
  
  return(tabela)
}

jpeg(file=paste0("Resultados/", nome.plot,"curvas.jpg"), width = 1500, height = 800, quality = 100, pointsize = 20)

par(mfrow=c(2,2), mar=c(2,2,2,2))

tabela.300<-tabela.n(300)
tabela.500<-tabela.n(500)
tabela.1000<-tabela.n(1000)
tabela.2000<-tabela.n(2000)

dev.off()


tabela.inicial<-data.frame(theta=c(pii[-g],unlist(beta.verd),unlist(sigma2.verd)))
colnames(tabela.inicial)<-c("$\\theta$")
tabela.final<-cbind.data.frame(tabela.inicial,tabela.300,tabela.500,tabela.1000, tabela.2000)
row.names(tabela.final)<-paste0("$\\",c("p_1","beta_{10}","beta_{11}", "beta_{12}","beta_{20}","beta_{21}","beta_{22}", "sigma^2_1","sigma^2_2"),"$")


format_scientific_latex_aux<-function(x) 
{
  formatted <- formatC(x, format = "e", digits = 2)
  formatted <- gsub("e([+-])", " \\\\times 10^{\\1", formatted)
  formatted <- gsub("([0-9])$", "\\1}$", formatted)
  formatted<-paste0("$", formatted)
  return(formatted)
}
format_scientific_latex <- function(x) {
  ifelse(abs(x)<1e-2, format_scientific_latex_aux(x),round(x, digits = 3))
}

# Aplicar a formatação a todos os valores numéricos da tabela
tabela.final[]<- lapply(tabela.final, format_scientific_latex)

tabela_latex <- knitr::kable(tabela.final, caption = paste0("C1 - Bem separados"), format = "latex", escape = FALSE, booktabs=T) %>%
  add_header_above(c(" " = 2, "n=300" = 3,"n=500" = 3,"n=1000" = 3,"n=2000"=3)) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

tabela_latex



