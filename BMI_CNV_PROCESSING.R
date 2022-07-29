BMI_CNV_PROCESSING=function(snpdata,snpmap,wesdata,wesmap,nsample){

ncol=nsample
lrr=snpdata
log2Rchr=wesdata

##### combine wes and snp
map_all=rbind(snpmap,wesmap)
map_all=map_all[order(map_all[,3]),]


#####platform specific robust scaling
lrr_snp_norm=matrix(data=NA,nrow=dim(lrr)[1],ncol=ncol)
log2R_norm=matrix(data=NA,nrow=dim(log2Rchr)[1],ncol=ncol)
for (i in 1:82) {
  dd=lrr[,i]
  lrr_snp_norm[,i]=(dd-median(dd))/(quantile(dd)[4]-quantile(dd)[2])
  dd1=log2Rchr[,i]
  log2R_norm[,i]=(dd1-median(dd1))/(quantile(dd1)[4]-quantile(dd1)[2])
}


#####combine lrr and log2R
LRR=matrix(data=NA,nrow=dim(map_all)[1],ncol=82)
wes_index=which(map_all$platform=="WES")
snp_index=which(map_all$platform=="SNP")
LRR[snp_index,]=lrr_snp_norm
LRR[wes_index,]=log2R_norm
LRR=t(LRR)

####generate cluster matrix

library(mixtools)
cluster=matrix(data=NA,nrow = dim(LRR)[1],ncol = dim(LRR)[2])
wan=vector()
for(j in 1:dim(LRR)[2]){
  x=LRR[,j]  
  wait1 <- normalmixEM(x, lambda = 1/3, mu = c(-2, 0, 1.5), sigma =1,maxit = 4000)
  iter=length(wait1$all.loglik)
  wan[j]=ifelse(iter==4001,1,0)
  if(wan[j]==1){
    index=rep(2,dim(LRR)[1])
    index[which(x< -1.2)]=1
    index[which(x> 1)]=3
  }else{
    prob=wait1$posterior
    index=vector()
    for(i in 1:dim(LRR)[1]){
      index[i]=which(prob[i,]==max(prob[i,]))
    }
    if(length(wait1$mu)==3 & max(wait1$mu)<0){
      mu1=wait1$mu[1]
      mu2=wait1$mu[3]
      if(mu1> -1.2){
        index=rep(2,dim(LRR)[1])
      }
      if(mu2> -1.2){
        index[index==3]=2
      }else{
        index=rep(1,dim(LRR)[1])
      }
    }else if(length(wait1$mu)==3 & min(wait1$mu)>0){
      mu1=wait1$mu[1]
      mu2=wait1$mu[3]
      if(mu1<1){
        index[index==1]=2
      }else{
        index=rep(3,dim(LRR)[1])
      }
      if(mu2<1){
        index=rep(2,dim(LRR)[1])
      }
    }else if(length(wait1$mu)==3){
      mu1=wait1$mu[1]
      mu2=wait1$mu[3]
      if(mu1>-1.2){
        index[index==1]=2
      }
      if(mu2<0.9){
        index[index==3]=2
      }
    }
    if(length(wait1$mu)==2){
      mu1=wait1$mu[1]
      mu2=wait1$mu[2]
      if(mu1>-1.2){
        index[index==1]=2
      }
      if(mu2<0.9){
        index[index==3]=2
      }
    }}
  cluster[,j]=index
}

return(list(LRR=LRR,cluster=cluster))
}







