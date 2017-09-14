   
rm(list=ls())
##install.packages(c('nloptr','parallel','data.table','ggplot2'))
library(parallel)
library(nloptr)
library(ggplot2)
library(data.table)
getwd()
source('./npbin.r')

#### ### We will analyze a subset of one of the sample dataset for illustration purpose.

thrc <- 5     #### minimum total coverage allowed
ncores <- 1   #### the number of cores to be used, can ONLY be 1 if run on Windows.  

dt = read.table('data_atac.txt',header=T)
colnames(dt)

dt.ct = data.table(dt)
dt.ct = dt.ct[m>=thrc,]  ### remove those with low total coverage

dt.ct = dt.ct[1:2000,]  ### for illustration purpose, keep only the 2000 points, remove this line will end up with analyzing the whole data set. It could be slow if only one core is used.
dt.ct[,phat:=xm/m] 
n = nrow(dt.ct)
#### NPBin
nbk <- 11  ### number of breaks
k <- 4     ### order ot splines
breaks <- seq(0,1,length.out=nbk)
pi.init <- hist(dt.ct[,phat],breaks=seq(0,1,length.out=nbk+k-3),plot=F)$density ### initialized the weights using the histogram of phat
### estimate the overall model
emBinBspltemp <- emBinBspl(dt.ct[,xm],dt.ct[,m],breaks=breaks,k=k,pi.init=pi.init,ncores=ncores,err.max=1e-3,iter.max=200)  
### estimate the null model
modtemp <- estNull(dt.ct[,xm],dt.ct[,m],emBinBspltemp,init=NULL,iter.max=200,ncores=ncores,ub=rep(log(1e4),2),err.max=1e-4) 
dt.ct[,fnp:=modtemp[['f']]][,f0np:=modtemp[['f0']]][,locfdrnp:=modtemp[['locfdr']]][,fdrnp:=locfdr2FDR(locfdrnp)][,ranknp:= rank(locfdrnp,ties.method='max')] 

names(modtemp)

modtemp$coef.null  ### null parameters of the NPBin model

## Empirical Bayes test using phat
pct0 <- 0.45         
ebbetahat <- ebBeta(dt.ct[,xm],dt.ct[,m],dt.ct[,phat],breaks=breaks,k=k,pi.init=hist(dt.ct[,phat],breaks=seq(0,1,length.out=nbk+k-3),plot=F)$density,pct0=pct0,init=NULL,iter.max=200,err.max=1e-4,ncores=ncores)
 dt.ct[,fhat:=ebbetahat[['f']]][,f0hat:=ebbetahat[['f0']]][,locfdrhat:=ebbetahat[['locfdr']]][,fdrhat:=locfdr2FDR(locfdrhat)][,rankhat:= rank(locfdrhat,ties.method='max')]

names(ebbetahat)

ebbetahat$coef.null ### null parameters of EBE 

## Binomial test
p.bin <- sapply(1:n, function(y) binom.test(dt.ct[y,xm],dt.ct[y,m])$p.value)
dt.ct[,pvbin:=p.bin][,fdrbin:=p.adjust(pvbin,method='BH')][,rankbin:= rank(pvbin,ties.method='max')]





### Evaluate the results using motifs

dt.ct[,sum(potential_TP)]  ### number of potential TP defined based on motif
dt.ct[,sum(potential_FP)]  ### number of potential FP defined based on motif

### find the number of TP and FP in top ranked SNPs
dt.ct[,tpnp:=rank2nhit(ranknp,potential_TP)][,fpnp:=rank2nhit(ranknp,potential_FP)]
dt.ct[,tphat:=rank2nhit(rankhat,potential_TP)][,fphat:=rank2nhit(rankhat,potential_FP)]
dt.ct[,tpbin:=rank2nhit(rankbin,potential_TP)][,fpbin:=rank2nhit(rankbin,potential_FP)]



### plot the accuracy measure as in the main paper. 
### We presented a zoom-in version in the main paper to the top 20%, 
### because usually there are not many ALI SNPs.
### Note that the default of the demo only select a subset of the data for illustration purpose.
### Thus the figure may not an exact replica of the one in the paper.
### To reproduce the results in the paper, please use the whole dataset

cbfpalette <- c( "#D55E00", "#0072B2","#CC79A7",  "#009E73", "#E69F00", "#56B4E9", "#F0E442")
plotidac <- c(' NPB','EBE','Binom')

tfzoomin.gg <-  ggplot(dt.ct) +
  geom_line(aes(x=sort(ranknp)/n,y=log2(tpnp[order(ranknp)]+1)-log2(fpnp[order(ranknp)]+1),fill=plotidac[1],colour=plotidac[1],linetype=plotidac[1]),size=1)+
  geom_line(aes(x=sort(rankhat)/n,y=log2(tphat[order(rankhat)]+1)-log2(fphat[order(rankhat)]+1),fill=plotidac[2],colour=plotidac[2],linetype=plotidac[2]),size=1)+
  geom_line(aes(x=sort(rankbin)/n,y=log2(tpbin[order(rankbin)]+1)-log2(fpbin[order(rankbin)]+1),fill=plotidac[3],colour=plotidac[3],linetype=plotidac[3]),size=1)+
  theme(legend.position = c(0.7, 0.25),legend.title=element_blank(),legend.text = element_text(size = 12, face = "bold"),
        legend.background = element_rect(colour = 'black', size = 0.2, linetype='solid'),legend.box = "horizontal")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(axis.text.x = element_text(size = 10, angle = 0, hjust = 0, colour = "black"))+
  xlab("Proportion of selected SNPs")+ylab("Log2(TP/FP)")+ ##xlim(0,0.2)+ylim(2.5,4.5)+
  scale_fill_manual(values = cbfpalette)+
  scale_colour_manual(values = cbfpalette)

print(tfzoomin.gg)




