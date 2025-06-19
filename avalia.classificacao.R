jpeg(file=paste0("Resultados/", nome.plot,"boxplots-acuracias.jpg"), width = 800, height = 800, quality = 100, pointsize = 20)

acuracia.n<-function(n) # Função para ler os arquivos e pegar os valores estimados para as M replicas 
{
  arquivo<-paste0(nome.amostra,n,".RDS")
  amostras<-readRDS(arquivo)
  acuracia.aux<-function(i){
    grupo.verd<-apply(amostras[[i]]$amostra$clu,2,FUN=function(x)which(as.logical(x)))
    grupo.em<-apply(amostras[[i]]$saida.EM$z,1,FUN=function(x)which.max(x))
    tabela<-table(grupo.verd, grupo.em)
    acuracia<-sum(diag(tabela))/sum(tabela)
    return(acuracia)
  }
  acuracias<-sapply(1:M, FUN=acuracia.aux)
}
l<-lapply(sizes,FUN=acuracia.n)
boxplot(l,boxwex=0.7,main="Acurácias", pch=1, col=cores,xaxt="n", ylim=c(0,1))
axis(1, at = seq_along(sizes), labels = paste0("n = ",sizes)) # Personaliza o eixo X com intervalos de 3

dev.off()

