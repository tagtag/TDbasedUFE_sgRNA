source("TD.R")
#yusa
x <- read.csv("yusa_raw_v10.tab.gz",sep="\t",header=T)
ysel<- split(x[,1],x[,2])
ysel_SPT <- ysel[lapply(ysel,length)>=5]
ysel_SPT_5 <- lapply(ysel_SPT,function(x){x[1:5]})
ysel_SPT_5_df <- t(data.frame(ysel_SPT_5))
rownames(ysel_SPT_5_df) <- names(ysel_SPT_5)

Z0<- array(NA,c(dim(ysel_SPT_5_df)[1],5,dim(x)[2]-2))
for (i in  seq_len(dim(ysel_SPT_5_df)[1]))
{
  cat(i," ")
  Z0[i,,] <-data.matrix(x[match(ysel_SPT_5_df[i,],x[,1]),-c(1:2)]) 
}
TD_result <-TD(Z0,ysel_SPT_5_df,L1=10,L2=5,ntop=10)
TD_result <-TD(Z0,ysel_SPT_5_df,log=T,L1=10,L2=5,ntop=10)


#TKOv1
x <- read.csv("tko_counts.txt.gz",sep="\t",header=T)
ysel<- split(x[,1],x[,2])
ysel_SPT <- ysel[lapply(ysel,length)>=4]
ysel_SPT_4 <- lapply(ysel_SPT,function(x){x[1:4]})
ysel_SPT_4_df <- t(data.frame(ysel_SPT_4))
rownames(ysel_SPT_4_df) <- names(ysel_SPT_4)

Z0<- array(NA,c(dim(ysel_SPT_4_df)[1],4,dim(x)[2]-2))
for (i in  seq_len(dim(ysel_SPT_4_df)[1]))
{
  cat(i," ")
  Z0[i,,] <-data.matrix(x[match(ysel_SPT_4_df[i,],x[,1]),-c(1:2)]) 
}
TD_result <-TD(Z0,ysel_SPT_4_df,L1=20,L2=4,ntop=10)
TD_result <-TD(Z0,ysel_SPT_4_df,log=T,L1=20,L2=4,ntop=10)

#Whitehead
x <- read.delim("~/RESEARCH/CRISPER/Wang2015_2017_merged_counts.txt.gz")
ysel<- split(x[,1],unlist(lapply(strsplit(x[,1],"_"),"[",1)))
ysel_SPT <- ysel[lapply(ysel,length)==10]
ysel_SPT_df <- t(data.frame(ysel_SPT))
rownames(ysel_SPT_df) <- names(ysel_SPT)
rownames(ysel_SPT_df) <- gsub("sg","",rownames(ysel_SPT_df))
Z0<- array(NA,c(dim(ysel_SPT_df)[1],10,dim(x[,-1])[2]))
for (i in  seq_len(dim(ysel_SPT_df)[1]))
{
  cat(i," ")
  Z0[i,,] <-data.matrix(x[match(ysel_SPT_df[i,],x[,1]),-1]) 
}
TD_result <-TD(Z0,ysel_SPT_df,L1=20,L2=5,ntop=10)
TD_result <-TD(Z0,ysel_SPT_df,log=T,L1=20,L2=5,ntop=10)

#GeCKO
x <- read.delim("~/RESEARCH/CRISPER/Achilles_raw_GeckoV2.tab.gz")
ysel<- split(x[,1],x[,2])
ysel_SPT <- ysel[lapply(ysel,length)>=4]
ysel_SPT_4 <- lapply(ysel_SPT,function(x){x[1:4]})
ysel_SPT_4_df <- t(data.frame(ysel_SPT_4))
rownames(ysel_SPT_4_df) <- names(ysel_SPT_4)
cell_names <- unlist(lapply(strsplit(toupper(colnames(x)),"_",fixed=T,),"[",1))
SPT <- split(seq_len(dim(x)[2]),cell_names)
index <- lapply(SPT,length)>=4
SPT <- SPT[index]
SPT <- lapply(SPT,function(x){x[1:4]})
SPT <- t(data.frame(SPT))
Z0<- array(NA,c(dim(ysel_SPT_4_df)[1],4,dim(SPT)[1]*4))
for (i in  seq_len(dim(ysel_SPT_4_df)[1]))
{
  cat(i," ")
  Z0[i,,] <-data.matrix(x[match(ysel_SPT_4_df[i,],x[,1]),SPT]) 
}
dim(Z0) <- c(18862,4,29,4)
TD_result <-TD(Z0,ysel_SPT_4_df,L1=10,L2=4,ntop=10)
TD_result <-TD(Z0,ysel_SPT_4_df,log=T,L1=10,L2=4,ntop=10)

#Avana
x <- read.csv("Avana_sgrna_raw_readcounts_matched.csv.gz",sep=",",header=T)
y <- read.csv("Avana_sgrnamapping.csv",sep=",",header = T)
ysel<- split(y[,1],y[,3])
ysel_SPT <- ysel[lapply(ysel,length)>=4]
ysel_SPT_4 <- lapply(ysel_SPT,function(x){x[1:4]})
ysel_SPT_4_df <- t(data.frame(ysel_SPT_4))
rownames(ysel_SPT_4_df) <- names(ysel_SPT_4)

Z <- array(NA,c(dim(ysel_SPT_4_df)[1],dim(x)[2],4))
for (i in  seq_len(dim(Z)[1]))
{
  cat(i," ")
  Z[i,,] <- t(data.matrix(x[match(ysel_SPT_4_df[i,],x[,1]),]))
}
Z <- Z[,-1,]

cell_names <- unlist(lapply(strsplit(toupper(colnames(x)),"311CAS9",fixed=T,),"[",1))
SPT <- split(seq_len(dim(Z)[2]),cell_names)
index <- lapply(SPT,length)>=2
SPT <- SPT[index]
SPT <- lapply(SPT,function(x){x[1:2]})
SPT <- t(data.frame(SPT))
Z0 <- array(NA,c(dim(Z)[1],dim(SPT)[1],2,4))
for (i in  seq_len(dim(SPT)[1]))
{
  cat(i," ")
  Z0[,i,,] <- Z[,SPT[i,],]
}
TD_result <-TD(Z0,ysel_SPT_4_df,L1=10,L2=10,L3=2,ntop=10)
TD_result <-TD(Z0,ysel_SPT_4_df,log=T,L1=10,L2=10,L3=2,ntop=10)

