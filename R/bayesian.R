# Update hyperparameters
hyper_update <- function(hyper.update, wh, hyper, Niter=100, Tol=1e-4){
    
   if(sum(hyper.update)==0) return(hyper)
  
   aw0 <- hyper$aw
   ah0 <- hyper$ah
   lwm <- mean(log(wh$lw))
   lhm <- mean(log(wh$lh))
   ewm <- mean(wh$ew)
   ehm <- mean(wh$eh)
   bw0 <- hyper$bw
   bh0 <- hyper$bh
  
   if(hyper.update[1]+hyper.update[3]>0){
    i <- 1
    while(i<Niter){
      if(hyper.update[1])
        dw <- (log(aw0)-digamma(aw0)-ewm/bw0+1+lwm-log(bw0))/
          (1/aw0-psigamma(aw0,1))
      else dw <- 0
      if(hyper.update[3])
        dh <- (log(ah0)-digamma(ah0)-ehm/bh0+1+lhm-log(bh0))/
          (1/ah0-psigamma(ah0,1))
      else dh <- 0
      aw1 <- aw0 - dw 
      ah1 <- ah0 - dh
      while(aw1<=0){
        dw <- dw/2
        aw1 <- aw0 - dw 
      }
      while(ah1<=0){
        dh <- dh/2
        ah1 <- ah0 - dh
      }
      
      df <- mean((1-aw1/aw0)^2)+mean((1-ah1/ah0)^2)
      if(df<Tol) break
      aw0 <- aw1
      ah0 <- ah1
      i <- i+1
    }
    if(i==Niter) stop('Hyper-parameter update failed to converge')
   } else{
     aw1 <- aw0
     ah1 <- ah0
   }
   if(hyper.update[2]) bw1 <- ewm
   else bw1 <- bw0
   if(hyper.update[4]) bh1 <- ehm
   else bh1 <- ehm
   list(aw=aw1, bw=bw1, ah=ah1, bh=bh1)
}

# Single update step in Bayesian NMF inference
vbnmf_updateR <- function(x, wh, r, estimator, hyper, fudge=NULL){

   x <- as.matrix(x)
   n <- dim(x)[1]
   m <- dim(x)[2]

   lw <- as.matrix(wh$lw)
   lh <- as.matrix(wh$lh)
   ew <- as.matrix(wh$ew)
   eh <- as.matrix(wh$eh)
   aw <- hyper$aw
   bw <- hyper$bw
   ah <- hyper$ah
   bh <- hyper$bh
  
   wth <- lw %*% lh
   sw <- lw*((x/wth) %*% t(lh))
   sh <- lh*(t(lw) %*% (x/wth))
  
   alw <- aw + sw
   bew <- 1/(aw/bw + t(replicate(n, rowSums(eh))))
   ew <- alw*bew             # this update needs to precede lines below
  
   alh <- ah + sh
   beh <- 1/(ah/bh + replicate(m, colSums(ew)))
   eh <- alh*beh

   lw <- exp(digamma(alw))*bew
   lh <- exp(digamma(alh))*beh
   if(is.null(fudge)) fudge <- .Machine$double.eps
   lw[lw < fudge] <- fudge
   lh[lh < fudge] <- fudge
    
   wth <- lw %*% lh
   U1 <- -ew %*% eh - lgamma(x+1) - x*((((lw*log(lw))%*%lh)
          +lw%*%(lh*log(lh)))/wth - log(wth))
   U2 <- -(aw/bw)*ew - lgamma(aw) + aw*log(aw/bw) + 
    alw*(1+log(bew))+lgamma(alw)
   U3 <- -(ah/bh)*eh - lgamma(ah) + ah*log(ah/bh) + 
    alh*(1+log(beh))+lgamma(alh)
   U <- sum(U1) + sum(U2) + sum(U3)
   U <- U/(n*m)  # log evidence per feature per cell
  
   if(estimator=='map'){
     w <- as.matrix(wh$w)  # MAP
     h <- as.matrix(wh$h)  # MAP
     w <- (aw - 1 + w*((x/(w%*%h))%*%t(h)))/(aw/bw + 
                                        t(replicate(n, rowSums(h))))
     h <- (ah - 1 + h*(t(w)%*%(x/(w%*%h))))/(ah/bh + 
                                        replicate(m, colSums(w)))
   }else{
     w <- ew
     h <- eh
   }
    
   list(w=w, h=h, lw=lw, lh=lh, ew=ew, eh=eh, lkh=U)
}

# Initialize bNMF inference
vb_init <- function(nrow,ncol,mat,rank, max=1.0, hyper){
  
   w <- matrix(stats::rgamma(n=nrow*rank, shape=hyper$aw, 
           scale=hyper$bw/hyper$aw), nrow=nrow,ncol=rank)
   rownames(w) <- rownames(mat)
   colnames(w) <- seq_len(rank)
   h <- matrix(stats::rgamma(n=rank*ncol, shape=hyper$ah, 
           scale=hyper$bh/hyper$ah), nrow=rank,ncol=ncol)
   rownames(h) <- seq_len(rank)
   colnames(h) <- colnames(mat)
  
   list(w=w, h=h, lw=w, lh=h, ew=w, eh=h)
}

#' Bayesian NMF inference of count matrix
#' 
#' Perform variational Bayes NMF and store factor matrices in object
#' 
#' The main input is the \code{scNMFSet} object with count matrix.
#' This function performs non-negative factorization using Bayesian algorithm
#' and gamma priors. Slots \code{basis}, \code{coeff}, and \code{ranks} 
#' are filled.
#' 
#' @param object \code{scNMFSet} object containing count matrix.
#' @param ranks Rank for factorization; can be a vector of multiple values.
#' @param nrun No. of runs with different initial guesses.
#' @param verbose The verbosity level: 
#'        3, each iteration output printed;
#'        2, each run output printed; 
#'        1, each randomized sample output printed;
#'        0, silent.
#' @param progress.bar Display progress bar with \code{verbose = 1} for 
#'       multiple runs.
#' @param estimator If \code{'mean'}, mean values of factor matrices 
#'        \code{W} and 
#'        \code{H} in E-step of EM algorithm are stored in \code{object}. 
#'        If \code{estimator = 'max'}, maximum a posteriori (MAP) 
#'        solutions are used instead.
#' @param Itmax Maximum no. of iteration.
#' @param hyper.update Vector of four logicals, each indcating whether
#'        hyperparameters \code{c(aw, bw, ah, bh)} should be optimized.
#' @param gamma.a Gamma distribution shape parameter.
#' @param gamma.b Gamma distribution mean. These two parameters are used for 
#'        fixed hyperparameters with \code{hyper.update} elements \code{FALSE}.
#' @param Tol Tolerance for terminating iteration.
#' @param hyper.update.n0 Initial number of steps in which hyperparameters 
#'        are fixed.
#' @param hyper.update.dn Step intervals for hyperparameter updates.
#' @param connectivity If \code{TRUE}, connectivity and dispersion will
#'        be calculated after each run. Can be turned off to save memory.
#' @param fudge Small positive number used as lower bound for factor matrix 
#'        elements to avoid singularity. If \code{fudge = NULL} (default), 
#'        it will be replaced by \code{.Machine$double.eps}. 
#'        Can be set to 0 to skip 
#'        regularization.
#' @return Object of class \code{scNMFSet} with factorization slots filled.
#' 
#' @details When run with multiple values of \code{ranks}, factorization is 
#'        repeated for each rank and the slot \code{measure} contains 
#'        log evidence and optimal hyperparameters for each rank. 
#'        With \code{nrun > 1}, the solution
#'        with the maximum log evidence is stored for a given rank.
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50,ncol=100,rank=5)
#' s <- scNMFSet(x$x)
#' s <- vb_factorize(s,ranks=seq(2,8),nrun=5)
#' plot(s)
#' @export
vb_factorize <- function(object, ranks=2, nrun=1, verbose=2, progress.bar=TRUE, 
        estimator='mean', Itmax=10000, hyper.update=rep(TRUE,4), 
        gamma.a=1, gamma.b=1, Tol=1e-5, hyper.update.n0=10, 
        hyper.update.dn=1, connectivity=TRUE, fudge=NULL){
  
   mat <- counts(object)  # S4 class scNMFSet
  
   nullr <- sum(Matrix::rowSums(mat)==0)
   nullc <- sum(Matrix::colSums(mat)==0)
   if(nullr>0) stop('Input matrix contains empty rows')
   if(nullc>0) stop('Input matrix contains empty columns')
  
   nrow <- dim(mat)[1]
   ncol <- dim(mat)[2]
  
   nrank <- length(ranks)  # rank values requested
  
   wdat <- vector('list',nrank)
   hdat <- vector('list',nrank)
   cons <- vector('list',nrank)
   awdat <- bwdat <- ahdat <- bhdat <- rdat <- rep(0, nrank)
  
   for(irank in seq_len(nrank)){
    
     rank <- ranks[irank]
     if(verbose > 0) cat('Rank ',rank,'\n',sep='')
    
     aw <- gamma.a[1]
     ah <- gamma.a[length(gamma.a)]
     bw <- gamma.b[1]
     bh <- gamma.b[length(gamma.b)]
    
     hyper <- hyper0 <- list(aw=aw, bw=bw, ah=ah, bh=bh)  
    # list of hyperparameter matrices
    
     if(connectivity){
       npair <- ncol*(ncol-1)/2
       conav <- rep(0,npair)
     }
     rmax <- -Inf
     if(verbose == 1 & progress.bar) 
       pb <- utils::txtProgressBar(style = 3)
     for(irun in seq_len(nrun)){
        
      if(verbose >=2 ) cat('Run #',irun,':\n')
      else if(verbose == 1 & progress.bar) 
        utils::setTxtProgressBar(pb, irun/nrun)
      hyper <- hyper0  
      wh <- vb_init(nrow, ncol, mat, rank, hyper=hyper)
      lk0 <- 0
      for(it in seq_len(Itmax)){
        wh <- vbnmf_updateR(mat, wh, rank, estimator=estimator, 
                            hyper, fudge=fudge)
        if(it > hyper.update.n0 & it%%hyper.update.dn==0) 
          hyper <- hyper_update(hyper.update, wh, hyper, Niter=100, 
                                Tol=1e-3)
        if(is.na(wh$lkh)) break
        if(it>1) if(wh$lkh>=lk0) if(abs(1-wh$lkh/lk0) < Tol) break
        lk0 <- wh$lkh
        if(verbose >= 3) cat(it,': log(evidence) = ',lk0,', aw = ',hyper$aw,
          ', bw = ',hyper$bw,', ah = ',hyper$ah,', bh = ',hyper$bh,
          '\n',sep='')
      }
      if(connectivity){
        cnn <- connectivity(wh$eh)
        conav <- conav + cnn
        disp <- dispersion(conav/irun,ncol)
      }
      if(verbose >= 2){
        if(connectivity) cat('Nsteps =',it,', log(evidence) =',lk0,
                  ', hyper = (',hyper$aw,',',hyper$bw,',',hyper$ah,',',
                  hyper$bh,')', ', dispersion = ',disp,'\n\n',sep='')
        else cat('Nsteps =',it,', log(evidence) =',lk0,
                 ', hyper = (',hyper$aw,',',hyper$bw,',',hyper$ah,',',
                 hyper$bh,')\n\n',sep='')
      }
      if((irun == 1 | lk0 > rmax) & !is.na(lk0)){
        rmax <- lk0
        if(estimator=='mean'){
          wmax <- wh$ew
          hmax <- wh$eh
        }else{
          wmax <- wh$w
          hmax <- wh$h
        }
        abmax <- hyper
      }
    }
    if(verbose == 1 & progress.bar) close(pb)
    
    wdat[[irank]] <- wmax
    hdat[[irank]] <- hmax
    rdat[irank] <- rmax
    awdat[irank] <- abmax$aw
    bwdat[irank] <- abmax$bw
    ahdat[irank] <- abmax$ah
    bhdat[irank] <- abmax$bh
    
    if(verbose >= 1) cat('Max(evidence) =',rmax,'\n\n')
   }   
   object@ranks <- ranks
   object@basis <- wdat
   object@coeff <- hdat
  
   measure(object) <- data.frame(rank=ranks, evidence=rdat, aw=awdat,
                      bw=bwdat, ah=ahdat, bh=bhdat)
   object  
}
