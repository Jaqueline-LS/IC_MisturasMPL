rm(list = ls(all = TRUE))
library("Matrix")
library("mgcv")
library("dplyr")
library("kableExtra")
library("scatterplot3d") 
source("funcoes/funcoes.R")
source("funcoes/funcoesMisturas.R")
source("funcoes/GerarAmostras_geral.R")

cores<-c("#c2a5cf","#BAF3DE","#C9E69E","#FF9B95","#FFC29A")

par<-par(pch=19)


M<-500
g<-2


pii<-c(0.35,0.65)

beta.verd<-list(c(4,4,6),c(2,3,-5))
nome.amostra<-"Amostras.corrigidas/c2_"
nome.plot<-"c2_"


#------------------------------
sigma2.verd<-list(1,1)


arg.grupos<-list(g1=list(beta=beta.verd[[1]],curva=list(f="c1",a=-1,b=1,d=2),sigma2=sigma2.verd[[1]],intercepto=T),
                 g2=list(beta=beta.verd[[2]],curva=list(f="c2",a=-1,b=1,d=4),sigma2=sigma2.verd[[2]],intercepto=T),
                 intervalos=list(c(0,1),c(0,1)))


n=300
alfas<-c(0.1,0.1)
amostra<-rMisUniPLM(n, pii, p=length(beta.verd[[1]]), arg=arg.grupos)
theta<-try(EMisUniPLM(g=2,alfas=alfas,y=amostra$y, t=amostra$t, X=amostra$X),TRUE)

verificarGeração(amostra, arg.grupos, tipo="Cenário 2")

y=amostra$y
t=amostra$t
X=amostra$X
grupo<- apply(amostra$clu,2,FUN=function(x)which(as.logical(x))) # numero do grupo de cada observação

k=2

par(mfrow=c(1,3), mar=c(2,2,1,1))

s3d <- scatterplot3d(z=y,x=X[,3],y=X[,2], 
                     xlab = "x3",
                     ylab = "x2",
                     type = "p", 
                     color = cores[grupo],
                     angle = 55, scale.y = 0.7, 
                     pch = c(17,15)[grupo], cex.symbols = 1.3,
                     main = paste0("Gráfico de dispersão 3D - C",k) )
s3d$plane3d(beta.verd[[1]][1],x.coef = beta.verd[[1]][3], y.coef = beta.verd[[1]][2], col=cores[1], lwd=2)

s3d$plane3d(beta.verd[[2]][1],x.coef = beta.verd[[2]][3], y.coef = beta.verd[[2]][2], col=cores[2], lwd=2)

legend("topright", legend = c("Grupo 1","Grupo 2"), fill = cores[1:2], cex=0.8)


s3d <- scatterplot3d(z=y,x=X[,3],y=t, 
                     xlab = "x3",
                     ylab = "t",
                     type = "p", 
                     color = cores[grupo],
                     angle = 55, scale.y = 0.7, 
                     pch = c(17,15)[grupo], cex.symbols = 1.3,
                     main = paste0("Gráfico de dispersão 3D - C",k) )
s3d$plane3d(mean(y),x.coef =0, y.coef = 0, lwd=2)



s3d <- scatterplot3d(z=y,x=X[,2],y=t, 
                     xlab = "x2",
                     ylab = "t",
                     type = "p", 
                     color = cores[grupo],
                     angle = 55, scale.y = 0.7, 
                     pch = c(17,15)[grupo], cex.symbols = 1.3,
                     main = paste0("Gráfico de dispersão 3D - C",k) )
s3d$plane3d(mean(y),x.coef =0, y.coef = 0, lwd=2)

VerificarEstimacao(amostra, arg.grupos, theta, tipo="Cenário 2")



# Replicações------------------------------------------------
alfas<-c(0.1,0.1)
sizes<-c(300,500,1000,2000)

# Esse comando gera as M amostras de tamanho 300 o arquivo é salvo na pasta Amostras 
# gera.amostras(300) 

# Para mais detalhes da geração das replicas das amostras confeir o arquivo "GerarAmostras_geral.R"


# Esse comando vai gerar as M amostras para cada tamanho de amostra definido em sizes
# Os arquivos são salvos na pasta Amostras
# sapply(sizes,FUN = gera.amostras) 


# Se quiser gerar os gráficos do MSE
source("funcoes/graficos_MSE.R")


# Para plotar as curvas e contruir as tabelas

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

tabela_latex <- knitr::kable(tabela.final, caption = paste0("C3 - Pouco separados"), format = "latex", escape = FALSE, booktabs=T) %>%
  add_header_above(c(" " = 2, "n=300" = 3,"n=500" = 3,"n=1000" = 3,"n=2000"=3)) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

tabela_latex


# Para gerar os Box-plots

source("funcoes/graficos_boxplots.R")
