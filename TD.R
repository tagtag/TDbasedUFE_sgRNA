TD <- function(Z0,ysel_SPT,log=F,c=1,L1=10,L2=10,L3=10,L4=4,ntop=10){
  if (min(Z0+c)<=0 & log) {print ("please increse the value of c");return()}
  if (log){
    Z0 <- apply(log10(Z0+c),seq_along(dim(Z0))[-1],scale)
    } else {
    Z0 <- apply(Z0,seq_along(dim(Z0))[-1],scale)
    }
  if (!(length(dim(Z0)) %in% 3:4)) {print("Tensor must be 3-way or 4-way");return()}
      require(rTensor)
      if (length(dim(Z0))==3)
      {
      HOSVD0 <- hosvd(as.tensor(Z0),c(L1,L2,L3))
      } else {
        HOSVD0 <- hosvd(as.tensor(Z0),c(L1,L2,L3,L4))  
      }
      library(readxl)
      es_gene <- read_excel("msb145216-sup-0004-datasets4.xlsx",sheet=2)
      es_gene <- es_gene[es_gene[,2]>=3,]
      es_gene_2 <- read_excel("msb145216-sup-0001-datasets1.xlsx")
      index1 <- rownames(ysel_SPT) %in% unlist(es_gene[,1])
      index2 <- rownames(ysel_SPT) %in% unlist(es_gene_2[,5])
      require(caTools)
      AUC <- NULL
      for (i in seq_len(L1))
      {
      AUC <- c(AUC,colAUC(c(HOSVD0$U[[1]][index1,i],HOSVD0$U[[1]][index2,i]),
                          c(rep(1,sum(index1)),rep(2,sum(index2)))))
      }
      l1<- which.max(AUC)
      if (length(dim(Z0))==4) {
      index <- which(abs(HOSVD0$Z@data[l1,,,])>=abs(HOSVD0$Z@data[l1,,,])[rank(-abs(HOSVD0$Z@data[l1,,,]))==ntop], arr.ind = T)
      } else if (length(dim(Z0))==3) {
      index <- which(abs(HOSVD0$Z@data[l1,,])>=abs(HOSVD0$Z@data[l1,,])[rank(-abs(HOSVD0$Z@data[l1,,]))==ntop], arr.ind = T)
      }
      index <- cbind(l1,index)
      index <- index[order(-abs(HOSVD0$Z@data[index])),]
      return(list(HOSVD0=HOSVD0,index=index))
}