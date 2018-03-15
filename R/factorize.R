# single update step of ML-NMF
nmf_updateR <- function(x,w,h,n,m,r, prior=FALSE, gamma.a, gamma.b){

  x <- as.matrix(x)
  w <- as.matrix(w)
  h <- as.matrix(h)

  up <- h*(t(w) %*% (x/(w %*% h)))
  down <- matrix(rep(colSums(w),m),nrow=r,ncol=m,byrow=FALSE)
  if(prior){
    up <- up + gamma.a-1
    down <- down + gamma.a/gamma.b
  }
  h <- up/down
  h[h < .Machine$double.eps] <- .Machine$double.eps
  
  up <- w*((x/(w %*% h)) %*% t(h))
  down <- matrix(rep(rowSums(h),n),nrow=n,ncol=r,byrow=TRUE)
  if(prior){
    up <- up + gamma.a-1
    down <- down + gamma.a/gamma.b
  }
  w <- up/down
  w[w < .Machine$double.eps] <- .Machine$double.eps
  
  list(ew=w,eh=h)
}

# initialize w and h
init <- function(nrow,ncol,mat,rank,max=1.0){
  w <- matrix(stats::runif(n=nrow*rank),nrow=nrow,ncol=rank)
  rownames(w) <- rownames(mat)
  colnames(w) <- 1:rank
  h <- matrix(stats::runif(n=rank*ncol),nrow=rank,ncol=ncol)
  rownames(h) <- 1:rank
  colnames(h) <- colnames(mat)
  list(ew=w,eh=h)
}

likelihood <- function(mat,w,h){

  wh <- as.vector(w %*% h)
  amat <- as.vector(mat)
  iz <- amat > 0 & wh > 0
  x <- sum(amat[iz]*log(wh[iz]/amat[iz])-wh[iz]+amat[iz])
  x/nrow(mat)/ncol(mat)

}

connectivity <- function(h){

# cid <- apply(h,2,function(x){which(x==max(x))})
  ncol <- dim(h)[2]
  cid <- c()
# for(j in 1:ncol) cid <- c(cid, which(h[,j]==max(h[,j])))
  for(j in 1:ncol) cid <- c(cid, which(h[,j]==max(h[,j]))[1])
  cnn <- outer(cid,cid,'==')
  t(cnn)[lower.tri(t(cnn))]

}

dispersion <- function(cnn,nc){

  con <- sum((cnn-0.5)^2)
  1/nc+8*con/nc^2

}

cophenet <- function(conav, nc){
  
  tmp <- matrix(0,nrow=nc,ncol=nc)
  tmp[lower.tri(tmp)] <- 1-conav
  d <- stats::as.dist(tmp)
  h <- stats::hclust(d)
  coph <- stats::cophenetic(h)
  stats::cor(d,coph)

}

#' Maximum likelihood factorization with randomization
#' 
#' Updates \code{scNMFSet object} by performing factorizations with
#' shuffled data matrices.
#' 
#' @param object \code{scNMFSet} object containing 
#'   count matrix and factorization results
#' @param ... Other parameters to \code{factorize()}.
#' @return Input \code{object} with slot \code{dmeasure} containing 
#'        differences in quality measures between those for original data
#'        and H0 results. 
#'        
#' @details This function updates the main object of class \code{scNMFSet} 
#'          by performing factorizations of randomized count data. It is assumed
#'          that main factorization of the original data has already been 
#'          performed. The same set of rank values in \code{ranks} slot of 
#'          \code{object} is used for factorization. The function 
#'          \code{factorize} 
#'          is then called with appropriate parameters. The resulting quality
#'          measures are stored in the slot
#'          \code{dmeasure}, a data frame with columns \code{likelihood} 
#'          and \code{dispersion}.
#' @examples
#' set.seed(1)
#' x <- simulate_data(nfeatures=10,nsamples=c(20,20,60,40,30))
#' s <- scNMFSet(x)
#' s <- factorize(s,ranks=2:8)
#' s <- factorize_H0(s)
#' plot(s)
#' @export
factorize_H0 <- function(object, ...){
                        
  m1 <- measure(object)
  ranks <- ranks(object)
  h0 <- object
  h0 <- factorize(h0, ranks = ranks, randomize = TRUE, ...)
  m0 <- measure(h0)
  dm <- data.frame(ranks=ranks,likelihood=m0$likelihood,
     r_se=m0$r_se, dispersion=m0$dispersion, d_se=m0$d_se,
     cophenetic=m0$cophenetic, c_se=m0$c_se)
  dmeasure(object) <- dm
  object
}

#' Maximum likelihood factorization
#' 
#' Performs single or multiple rank NMF factorization of count matrix using 
#' maximum likelihood
#' 
#' The main input is the \code{scNMFSet} object with count matrix.
#' This function performs non-negative factorization and fills in the empty
#' slots \code{basis}, \code{coeff}, and \code{ranks}.
#' 
#' @param object \code{scNMFSet} object containing count matrix.
#' @param ranks Rank for factorization; can be a vector of multiple values.
#' @param nrun No. of runs with different initial guess.
#' @param randomize Boolean; if \code{TRUE}, input matrix is randomized.
#' @param nsmpl No. of randomized samples to average over.
#' @param verbose The verbosity level: 
#'        3, each iteration output printed;
#'        2, each run output printed; 
#'        1, each randomized sample output printed;
#'        0, silent.
#' @param progress.bar Display progress bar when \code{nrun > 1 } and 
#'        \code{verbose = 1}.
#' @param Itmax Maximum no. of iteration.
#' @param ncnn.step Minimum no. of steps with no change in connectivity matrix 
#'        to achieve convergence.
#' @param criterion If \code{'likelihood'}, iteration stops when fractional 
#'        changes in likelihood is below tolerance \code{Tol}. If 
#'        \code{criterion = 'connectivity'}, iteration stops when connectivity 
#'        matrix does not change for at least \code{ncnn.step} steps.
#' @param Tol Tolerance for checking convergence with 
#'       \code{criterion = 'likelihood'}.
#'        
#' @return Object of class \code{scNMFSet} with factorization slots filled.
#' 
#' @details When run with multiple values of \code{ranks}, 
#'         factorization is repeated for each rank and the slot \code{measure} 
#'         contains quality measures of the ranks. The quality measure 
#'         \code{likelihood} is negative the KL distance of the fit to the 
#'         target. With \code{nrun > 1}, the likelihood is the maximum 
#'         among all runs. The quality measure \code{dispersion} is the scalar
#'         measure of how far the connectivity matrix is from 0, 1. With 
#'         increasing \code{nrun}, \code{dispersion} decreases from 1. 
#'         \code{nrun} should be chosen such that \code{dispersion} does not 
#'         change appreciably.
#'         With randomization, \code{count} matrix of \code{object} 
#'         is shuffled. 
#'         \code{nsmpl} can be used to average over multiple permutations. This 
#'         averaging applies to each quality measure under a given rank. The 
#'         randomized factorization is called by \code{\link{factorize_H0}} 
#'         to update \code{object}.
#' @examples
#' set.seed(1)
#' x <- simulate_data(nfeatures=10,nsamples=c(20,20,60,40,30))
#' s <- scNMFSet(x)
#' s <- factorize(s,ranks=2:8)
#' plot(s)
#' @export
factorize <- function(object, ranks=2, nrun=20, randomize=FALSE, nsmpl=1, 
      verbose=2, progress.bar=TRUE, Itmax=10000, ncnn.step=40, 
      criterion='likelihood', Tol=1e-5){
  
  mat <- count(object)  # S4 class scNMFSet
  
  nullr <- sum(rowSums(mat)==0)
  nullc <- sum(colSums(mat)==0)
  if(nullr>0) stop('Input matrix contains empty rows')
  if(nullc>0) stop('Input matrix contains empty columns')
  
  nrow <- dim(mat)[1]
  ncol <- dim(mat)[2]
  
  nrank <- length(ranks)  # rank values requested

  wdat <- vector('list',nrank)
  hdat <- vector('list',nrank)
  npair <- ncol*(ncol-1)/2
  rdat <- ddat <- cdat <- vector('list',nrank)
  rave <- dave <- coav <- rste <- cste <- dste <- rep(0,nrank)
  if(randomize) mat0 <- mat
  
  for(irank in 1:nrank){
    
    rank <- ranks[irank]
    if(verbose > 0) cat('Rank ',rank,'\n',sep='')
    
    for(ismpl in 1:nsmpl){
      
      conav <- rep(0,npair)
      if(randomize) mat <- apply(mat0,2,
                     function(x){sample(x,size=length(x),replace=FALSE)})
      if(verbose == 1 & progress.bar) 
        pb <- utils::txtProgressBar(style = 3)
      rmax <- -Inf
      for(irun in 1:nrun){
      
        if(verbose >=2 ){
          if(randomize)
            cat('Rnd.sample #',ismpl,', run #',irun,':\n')
          else
            cat('Run #',irun,':\n')
        }else if(verbose == 1 & progress.bar) 
          utils::setTxtProgressBar(pb, irun/nrun)
        
        wh <- init(nrow,ncol,mat,rank)
        
        zstep <- 0; lkold <- -Inf
        for(it in 1:Itmax){
          wh <- nmf_updateR(mat,wh$ew,wh$eh,nrow,ncol,rank)
          lk0 <- likelihood(mat,wh$ew,wh$eh)
          if(criterion=='connectivity'){
            cnn <- connectivity(wh$eh)
            if(it == 1) nchange <- npair
            else nchange <- sum(cnn!=cnn0)
          
            if(verbose >= 3)
              cat(it,': likelihood = ',lk0,', connectivity change = ', 
                  nchange,'\n')
            if(nchange == 0) zstep <- zstep+1
            else zstep <- 0
            if(zstep == ncnn.step) break
            cnn0 <- cnn
          }
          else if(criterion=='likelihood'){
            if(abs(lkold-lk0)<Tol*abs(lkold)) break
            if(verbose >=3) cat(it,': likelihood = ',lk0,'\n')
            lkold <- lk0
          }
          else stop('Unknown stopping criterion.')
        }
        cnn <- connectivity(wh$eh)
        conav <- conav + cnn
        disp <- dispersion(conav/irun,ncol)
        if(verbose >= 2) cat('Nsteps =',it,
                ', likelihood =',lk0,', dispersion =',disp,'\n\n')
        if((irun == 1 | lk0 > rmax) & !is.na(lk0)){
          rmax <- lk0
          wmax <- wh$ew
          hmax <- wh$eh
        }
      }
      if(verbose == 1 & progress.bar) close(pb)
      coph <- cophenet(conav/nrun,ncol)
      if(verbose >= 1) cat('Sample#',ismpl,': Max(likelihood) =',rmax,
                ', dispersion =',disp,', cophenetic =',coph,'\n')
      if(ismpl==1){
        wdat[[irank]] <- wmax
        hdat[[irank]] <- hmax
        rdat[[irank]] <- rmax
        ddat[[irank]] <- disp
        cdat[[irank]] <- coph
      } else{
        wdat[[irank]] <- wdat[[irank]]+wmax
        hdat[[irank]] <- hdat[[irank]]+hmax
        rdat[[irank]] <- c(rdat[[irank]], rmax)
        ddat[[irank]] <- c(ddat[[irank]], disp)
        cdat[[irank]] <- c(cdat[[irank]], coph)
      }
    }
    wdat[[irank]] <- wdat[[irank]]/nsmpl
    hdat[[irank]] <- hdat[[irank]]/nsmpl
    if(nsmpl>1){
      denom <- sqrt(nsmpl-1)
      rste[irank] <- stats::sd(rdat[[irank]])/denom
      dste[irank] <- stats::sd(ddat[[irank]])/denom
      cste[irank] <- stats::sd(cdat[[irank]])/denom
    }
    else{
      rste[irank] <- dste[irank] <- cste[irank] <- NA
    }
    rave[irank] <- mean(rdat[[irank]])
    dave[irank] <- mean(ddat[[irank]])
    coav[irank] <- mean(cdat[[irank]])
    if(verbose >= 1 & randomize) cat('Mean(likelihood) = ',rave[irank],
                    ', Mean(dispersion) =',dave[irank],
                    ', Mean(cophenetic) =',coav[irank],'\n\n')
  } 
  ranks(object) <- ranks
  basis(object) <- wdat
  coeff(object) <- hdat
  if(randomize) measure(object) <- data.frame(rank=ranks, likelihood=rave, 
              r_se=rste, dispersion=dave, d_se=dste, cophenetic=coav, c_se=cste)
  else measure(object) <- data.frame(rank=ranks, likelihood=rave, 
                                     dispersion=dave, cophenetic=coav)
  object  
}