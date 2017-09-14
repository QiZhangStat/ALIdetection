
library(parallel)
library(nloptr)
###library(VGAM)


betabinom.ab <- function (x, size, shape1, shape2, log = FALSE, Inf.shape = exp(20), 
    limit.prob = 0.5) { ### This function was taken from VGAM, so that the code does not require VGAM package anymore
    Bigg <- Inf.shape
    Bigg2 <- Inf.shape
    if (!is.logical(log.arg <- log) || length(log) != 1) 
        stop("bad input for argument 'log'")
    rm(log)
    LLL <- max(length(x), length(size), length(shape1), length(shape2))
    if (length(x) != LLL) 
        x <- rep_len(x, LLL)
    if (length(size) != LLL) 
        size <- rep_len(size, LLL)
    if (length(shape1) != LLL) 
        shape1 <- rep_len(shape1, LLL)
    if (length(shape2) != LLL) 
        shape2 <- rep_len(shape2, LLL)
    is.infinite.shape1 <- is.infinite(shape1)
    is.infinite.shape2 <- is.infinite(shape2)
    ans <- x
    ans[TRUE] <- log(0)
    ans[is.na(x)] <- NA
    ans[is.nan(x)] <- NaN
    ok0 <- !is.na(shape1) & !is.na(shape2) & !is.na(x) & !is.na(size)
    okk <- (round(x) == x) & (x >= 0) & (x <= size) & !is.infinite.shape1 & 
        !is.infinite.shape2 & ok0
    if (any(okk)) {
        ans[okk] <- lchoose(size[okk], x[okk]) + lbeta(shape1[okk] + 
            x[okk], shape2[okk] + size[okk] - x[okk]) - lbeta(shape1[okk], 
            shape2[okk])
        endpt1 <- (x == size) & ((shape1 < 1/Bigg) | (shape2 < 
            1/Bigg)) & ok0
        if (any(endpt1)) {
            ans[endpt1] <- lgamma(size[endpt1] + shape1[endpt1]) + 
                lgamma(shape1[endpt1] + shape2[endpt1]) - (lgamma(size[endpt1] + 
                shape1[endpt1] + shape2[endpt1]) + lgamma(shape1[endpt1]))
        }
        endpt2 <- (x == 0) & ((shape1 < 1/Bigg) | (shape2 < 1/Bigg)) & 
            ok0
        if (any(endpt2)) {
            ans[endpt2] <- lgamma(size[endpt2] + shape2[endpt2]) + 
                lgamma(shape1[endpt2] + shape2[endpt2]) - (lgamma(size[endpt2] + 
                shape1[endpt2] + shape2[endpt2]) + lgamma(shape2[endpt2]))
        }
        endpt3 <- ((Bigg < shape1) | (Bigg < shape2)) & ok0
        if (any(endpt3)) {
            ans[endpt3] <- lchoose(size[endpt3], x[endpt3]) + 
                lbeta(shape1[endpt3] + x[endpt3], shape2[endpt3] + 
                  size[endpt3] - x[endpt3]) - lbeta(shape1[endpt3], 
                shape2[endpt3])
        }
    }
    if (!log.arg) {
        ans <- exp(ans)
    }
    ok1 <- !is.infinite.shape1 & is.infinite.shape2
    ok2 <- is.infinite.shape1 & !is.infinite.shape2
    ok3 <- Bigg2 < shape1 & Bigg2 < shape2
    ok4 <- is.infinite.shape1 & is.infinite.shape2
    if (any(ok3)) {
        prob1 <- shape1[ok3]/(shape1[ok3] + shape2[ok3])
        ans[ok3] <- dbinom(x = x[ok3], size = size[ok3], prob = prob1, 
            log = log.arg)
        if (any(ok4)) {
            ans[ok4] <- dbinom(x = x[ok4], size = size[ok4], 
                prob = limit.prob, log = log.arg)
        }
    }
    if (any(ok1)) 
        ans[ok1] <- dbinom(x = x[ok1], size = size[ok1], prob = 0, 
            log = log.arg)
    if (any(ok2)) 
        ans[ok2] <- dbinom(x = x[ok2], size = size[ok2], prob = 1, 
            log = log.arg)
    ans[shape1 < 0] <- NaN
    ans[shape2 < 0] <- NaN
    ans
}





dbetabinom.vec <- function(x,m,rate1,rate2,log=TRUE){ ## vector of beta binomial density
  n <- length(x)
  lvec <- sapply(1:n, function(ii) dbetabinom.ab(x[ii],size=m[ii],shape1=rate1,shape2=rate2,log=log))
  out <- unlist(lvec)
  out
}


bsplfun <- function(xrange=c(0,1),breaks=seq(xrange[1],xrange[2],length.out=100),k=4,ncores=1){ # create B-spline 
          nbreaks <- length(breaks)
          nbasis <- nbreaks+k-2
          step.break <- mean(abs(diff(breaks)))
          breaks.extend <- c(breaks[1]-step.break*((k-1):1),breaks,breaks[nbreaks]+step.break*(1:(k-1)))
          coef0 <- matrix(1,1,1)
          bspl <- lapply(1:length(breaks.extend),function(ii) list(breaks=breaks.extend[ii:(ii+1)],coef=coef0))
          for(kk in 2:k){
                    bspl.kk <- mclapply(1:(length(bspl)-1), function(ii) bsplfun.updt(ii,bspl),mc.cores=ncores,mc.preschedule=FALSE)
                            bspl <- bspl.kk
                  }
           out <- bspl[1:nbasis]
           for(ii in 1:(k-1)){
              outii <- out[[ii]]
              breaksii <- outii[['breaks']][(k+1-ii):(k+1)]
              coefii <- outii[['coef']][,(k-ii+1):k]
              if(ii==1){coefii <- matrix(coefii,ncol=1)}
              out[[ii]] <- list(breaks=breaksii,coef=coefii)
            }
           for(ii in 1:(k-1)){
              outii <- out[[nbasis+1-ii]]
              breaksii <- outii[['breaks']][1:(ii+1)]
              coefii <- outii[['coef']][,1:ii]
              if(ii==1){coefii <- matrix(coefii,ncol=1)}
              out[[nbasis+1-ii]] <- list(breaks=breaksii,coef=coefii)
            }
           for(ii in 1:nbasis){
              outii <- out[[ii]]
              breaksii <- outii[['breaks']]
              coefii <- outii[['coef']]
              intii <- sum(sapply(1:(length(breaksii)-1), function(jj) sum((breaksii[jj+1]^(1:k)-breaksii[jj]^(1:k))*coefii[,jj]/(1:k))))
              coefii <- coefii/intii
              out[[ii]] <- list(breaks=breaksii,coef=coefii)
            }

           out
        }

bsplfun.updt <- function(ii,bspl){ 
      bs1 <- bspl[[ii]]
          breaks1 <- bs1[['breaks']]
          coef1 <- bs1[['coef']]
          bs2 <- bspl[[ii+1]]
          breaks2 <- bs2[['breaks']]
          coef2 <- bs2[['coef']]
          kk <- length(breaks1)       ## we assume that length(breaks1) and length(breaks2) are the same
          breaks <- c(breaks1[1],breaks2) ## we assume that breaks1[-1] == breaks2[1:(kk-1)]
          coef <- rbind(0,cbind(coef1,0))/(breaks1[kk]-breaks1[1]) - rbind(cbind(coef1,0),0)*breaks1[1]/(breaks1[kk]-breaks1[1]) -rbind(0,cbind(0,coef2))/(breaks2[kk]-breaks2[1]) + rbind(cbind(0,coef2),0)*breaks2[kk]/(breaks2[kk]-breaks2[1])
          out <- list(breaks=breaks,coef=coef)
          out
    }

iBiBsFun <- function(bk,x,m,cf){ ### integration over one spline basis fnction
  k <- length(cf)  
  out <- sapply(1:k, function(jj) cf[jj]*exp(lgamma(m+1)-lgamma(x+1)+lgamma(x+jj)-lgamma(m+jj+1))*(pbeta(bk[2],shape1=jj+x,shape2=m-x+1)-pbeta(bk[1],shape1=jj+x,shape2=m-x+1)))
  out
}


intBinBspl <- function(bs,x,m){ ## integration of Binomial + B spline 
  breaks <- bs[['breaks']]
  cf <- bs[['coef']]
  out <- sapply(1:(length(breaks)-1), function(ii) sum(iBiBsFun(breaks[ii:(ii+1)],x,m,cf[,ii])))
  sum(out)
}

getDesignMtx <- function(bs,x,m,ncores=1){ # get the design matrix c_ij
nbasis <- length(bs)
ndt <- length(x)
out <- mclapply(1:ndt, function(ii) sapply(bs, function(jj) intBinBspl(jj,x[ii],m[ii])),mc.cores=ncores,mc.preschedule=FALSE)
out <- t(as.matrix(as.data.frame(out)))
dimnames(out) <- NULL
out
}

emBinBspl <- function(x,m,breaks=seq(0,1,length.out=101),k=4,pi.init=rep(1,length(breaks)+k-4)/(length(breaks)+k-4),ncores=1,err.max=0.000001,iter.max=200){  ### NPBin model. We remove the basis that only covers bin at the boundary
  n <- length(x)
  bspl <- bsplfun(range(breaks),breaks=breaks,k=k,ncores=ncores) ## create the break points and coefficients of the B spline functions
  bspl <- bspl[-c(1,length(bspl))]
  nb <- length(bspl)
  dmtx <- getDesignMtx(bspl,x,m ,ncores=ncores) ### get the design matrix
  err <- Inf
  iter <- 0
  ll.init <- -Inf
  ll.all <- ll.init
  err.all <- err
  while(err>err.max&iter<iter.max){
    dtot <-  dmtx%*%matrix(pi.init,length(pi.init),1) 
    post <- dmtx*(matrix(pi.init,n,nb,byrow=TRUE))/matrix(dtot,n,nb,byrow=FALSE)
    pi <- colSums(post)
    pi <- pi/sum(pi)
    ll <- sum(apply(dmtx,1, function(ii) log(sum(ii*pi))))
    err <- ll-ll.init
    err <- max(err, err/abs(ll))
    ll.all <- c(ll.all,ll)
    err.all <- c(err.all,err)
    pi.init <- pi
    ll.init <- ll
    iter <- iter+1
  } 
  out <- list(pi=as.numeric(pi.init),post=post,bspl=bspl,dmtx=dmtx,f=as.numeric(dmtx%*%pi.init),ll.all=ll.all,err.all=err.all, convergence=list(err=err,err.max=err.max,converged=(err<err.max),niter=iter,ll=ll.init),controls=list(k=k,nbasis=nb,breaks=breaks))
  out
}



evBspl <- function(p,bspl,ncores=1){ ## evaluate B-spline function
nb <- length(bspl)
k <- dim(bspl[[1]][['coef']])[1]
 evBspl.sg <- function(p,bsplsg){
   coef <- bsplsg[['coef']]
   bk <- bsplsg[['breaks']]
   id <- which(p<bk)
  if(length(id)==0|id[1]==1){
    out <- 0
  }else{
    kk <- id[1]-1
    out <- sum(sapply(1:k, function(jj) coef[jj,kk]*p^(jj-1)))
  }
  out
 }
out <- mclapply(p, function(pp) sapply(bspl, function(ii) evBspl.sg(pp,ii)),mc.cores=ncores,mc.preschedule=FALSE)
out <- t(as.matrix(as.data.frame(out)))
dimnames(out) <- NULL
out
}





evBsplDrv <- function(p,bspl,ncores=1){  ## evaluate the derivative of B-spline
nb <- length(bspl)
k <- dim(bspl[[1]][['coef']])[1]
 evBsplDrv.sg <- function(p,bsplsg){
   coef <- bsplsg[['coef']]
   bk <- bsplsg[['breaks']]
   id <- which(p<bk)
  if(length(id)==0|id[1]==1){
    out <- 0
  }else{
    kk <- id[1]-1
    out <- sum(sapply(2:k, function(jj) (jj-1)*coef[jj,kk]*p^(jj-2)))
  }
  out
 }
out <- mclapply(p, function(pp) sapply(bspl, function(ii) evBsplDrv.sg(pp,ii)),mc.cores=ncores,mc.preschedule=FALSE)
out <- t(as.matrix(as.data.frame(out)))
dimnames(out) <- NULL
out
}





estNull1  <- function(mod,pseq=(1:9999)/1e4,ncores=1){  ## estimate the null model, part 1
    pseq <- (1:9999)/10000
     bspl <- mod[['bspl']]
     pi <- mod[['pi']]
     drvbspl <- evBsplDrv(pseq,bspl,ncores=ncores)
     df <- apply(drvbspl,1, function(ii) sum(ii*pi))
     evbspl <- evBspl(pseq,bspl,ncores=ncores)
     fseq <- apply(evbspl,1, function(ii) sum(ii*pi))
out <- list(fseq=fseq,df=df)
out
}

estNull2 <- function(x,m,mod,prep,init=NULL,iter.max=200,err.max=1e-6,algorithm='NLOPT_GN_DIRECT_L',pseq=(1:9999)/1e4,lb=c(0,0),ub=rep(log(1e4),2)){   ## estimate the null model, part 2
  if(is.null(init)){
  phat <- x/m
  m1p <- mean(phat)
    n <- length(x)
  ss <- n*var(phat)
  mi <- sum(1/m)
  m2p <- (ss-(m1p-m1p^2)*mi)/(n-mi)
    sc <- m1p^2/m2p-1
    init <- c(m1p*sc,(1-m1p)*sc)
}
    fseq <- prep[['fseq']]
    df <- prep[['df']]
   ell <- function(ipt){
        shape1 <- exp(ipt[1])
        shape2 <- exp(ipt[2])
        df0 <- (shape1-1-(shape1+shape2-2)*pseq)*dbeta(pseq,shape1-1,shape2-1,log=FALSE)*exp(lbeta(shape1-1,shape2-1)-lbeta(shape1,shape2))
        f0seq <- dbeta(pseq,shape1,shape2,log=FALSE)
        ovec <-  (fseq*df0-df*f0seq)^2/fseq^3 ##  
        ovec[is.na(ovec)] <- 0
        mean(ovec)
    }
    optout <- nloptr(log(init),ell,lb=lb,ub=ub,opts=list(algorithm=algorithm,maxeval=iter.max,ftol_rel=err.max,xtol_rel=sqrt(err.max)))
    coef.opt <- exp(optout[['solution']])
    f <- mod[['f']]
    f0 <- dbetabinom.vec(x,m,coef.opt[1],coef.opt[2],log=FALSE)
    pi0 <- min(1,1/quantile(f0/f,probs=0.975))
    maxi <- c(x[which.min(f/f0)],m[which.min(f/f0)])
    coef <- list(shape1=coef.opt[1],shape2=coef.opt[2],pi0=pi0)
    out <- mod
    out[['coef.null']] <- coef
    out[['pi0']] <- pi0
    out[['f0']] <- f0
    out[['locfdr']] <- pi0*f0/f
    out[['convergence.null']] <- list(opt.out=optout,p.maxlr=maxi)
    out
}

estNull  <- function(x,m,mod,init=NULL,ncores=1,iter.max=200,err.max=1e-6,algorithm='NLOPT_GN_DIRECT_L',pseq=(1:9999)/1e4,lb=c(0,0),ub=rep(log(1e4),2)){  ## estimate the null model
  prep <- estNull1(mod,pseq,ncores=ncores)
output <- estNull2(x,m,mod,prep,init=init,iter.max=iter.max,err.max=err.max,algorithm=algorithm,pseq=pseq,lb=lb,ub=ub)
  output
}


locfdr2FDR <- function(locfdr){  ### convert locfdr to FDR
n <- length(locfdr)
out <- sapply(1:n, function(ii) mean(locfdr[locfdr<=locfdr[ii]]))
out
}

rank2nhit <- function(r,id){  ### convert ranking to the number of discoveries
out <- sapply(r, function(y) sum((r<=y)&id))
out
}


#####################################################################################################################################################################################################################################################
##################### functions for methods to be compared with

 
betaTrim_mle <- function(x,m,p,pct0=0.25,init=c(1,1,0.5),iter.max=200,err.max=1e-6,lb.opt=c(0,0,0),ub.opt=c(log(1e4),log(1e4),1)){ ## estimate the null model following Efron's approach.
    if(is.null(init)){
       m1p <- mean(p)
       m2p <- var(p)
       sc <- m1p*(1-m1p)/m2p-1
       init <- c(log(m1p*sc),log((1-m1p)*sc),0.9)
    }
   beta_bulk <- function(input){
     rate1 <- exp(input[1])
     rate2 <- exp(input[2])
     pi0 <- input[3]
        if((rate1>0)&(rate2>0)&(pi0>0)&(pi0<1)){
           n <- length(p)
           ub <- quantile(p,probs=0.5+pct0)
           lb <- quantile(p,probs=0.5-pct0)
           ix <- (p<=ub)&(p>=lb)
           dbt <- dbeta(p,shape1=rate1,shape2=rate2,log=TRUE)
           dbm <- log(1-pi0*(pbeta(ub,shape1=rate1,shape2=rate2,log=FALSE)-pbeta(lb,shape1=rate1,shape2=rate2,log=FALSE)))
           betavec <- ix*(dbt+log(pi0))+(1-ix)*dbm
           betavec[is.na(betavec)] <- 0
           -sum(betavec)
         }else{Inf}
    }
    optout <- nloptr(init,beta_bulk,lb=lb.opt,ub=ub.opt,opts=list(algorithm='NLOPT_GN_DIRECT_L',maxeval=iter.max,ftot_abs=err.max,ftol_rel=err.max))
   coefout <-  optout[['solution']]
   rate1 <- exp(coefout[1])
   rate2 <- exp(coefout[2])
   pi0 <- coefout[3]
   f0 <- dbetabinom.vec(x,m,rate1=rate1,rate2=rate2,log=FALSE)
    out <- list(coef=c(rate1,rate2,pi0),pi0=pi0,f0=f0,optout=optout)
}


emBspl <- function(x,m,p,breaks=seq(0,1,length.out=101),k=4,pi.init=rep(1,length(breaks)+k-4)/(length(breaks)+k-4),ncores=1,err.max=1e-5,iter.max=200){  ### estimate the overall model directly via B-spline. remove the basis that only covers bin at the boundary
    p[p==0] <- 1/max(m)^3
    p[p==1] <- 1-1/max(m)^3
  n <- length(p)
  bspl <- bsplfun(range(breaks),breaks=breaks,k=k,ncores=ncores) ## create the break points and coefficients of the B spline functions
  bspl <- bspl[-c(1,length(bspl))]
  nb <- length(bspl)
  dmtx <- evBspl(p,bspl,ncores=ncores) ### get the design matrix
  err <- 1e5*err.max
  iter <- 0
  ll.init <- -Inf
  ll.all <- ll.init
  err.all <- err
  while((is.na(err)|err>err.max)&iter<iter.max){
    dtot <-  dmtx%*%matrix(pi.init,length(pi.init),1) ##### 
    dtot[dtot==0] <- min(dtot[dtot>0])
    post <- dmtx*(matrix(pi.init,n,nb,byrow=TRUE))/matrix(dtot,n,nb,byrow=FALSE)
    pi <- colSums(post)
    pi <- pi/sum(pi)
    ell <- apply(dmtx,1, function(ii) sum(ii*pi))
    ll <- sum(log(ell))
    err <- ll-ll.init
    err <- max(err, err/max(abs(ll.init),abs(ll)))
    ll.all <- c(ll.all,ll)
    err.all <- c(err.all,err)
    pi.init <- pi
    ll.init <- ll
    iter <- iter+1
  } 
  binmtx <- getDesignMtx(bspl,x,m ,ncores=ncores) ### get the design matrix
  out <- list(pi=as.numeric(pi.init),post=post,bspl=bspl,dmtx=dmtx,binmtx=binmtx,f=as.numeric(binmtx%*%pi.init),ll.all=ll.all,err.all=err.all, convergence=list(err=err,err.max=err.max,converged=(err<err.max),niter=iter,ll=ll.init),controls=list(k=k,nbasis=nb,breaks=breaks))
  out
}

 
ebBeta <- function(x,m,p,breaks=seq(0,1,length.out=101),k=4,pi.init=rep(1,length(breaks)+k-4)/(length(breaks)+k-4),pct0=0.45,init=NULL,iter.max=200,err.max=1e-5,ncores=1,lb.opt=c(0,0,0),ub.opt=c(log(1e4),log(1e4),1)){ ## wrapper of EBO and EBE
mod <- emBspl(x,m,p,breaks=breaks,k=k,pi.init=pi.init,ncores=ncores,err.max=err.max,iter.max=iter.max)
f <- mod[['f']]  
null <- betaTrim_mle(x,m,p,pct0,init=init,iter.max=iter.max,err.max=err.max,lb.opt,ub.opt)
f0 <- null[['f0']]
pi0 <- null[['pi0']]
locfdr <- pi0*f0/f
locfdrnorm <- locfdr/max(locfdr)
coefnull <- null[['coef']]
out <- mod
out[['coef.null']] <- list(shape1=coefnull[1],shape2=coefnull[2],pi0=pi0)
out[['f0']] <- f0
out[['f']] <- f
out[['pi0']] <- pi0
out[['locfdr']] <- locfdr
out[['locfdrnorm']] <- locfdrnorm
out[['null']] <- null[['optout']]
out
}

