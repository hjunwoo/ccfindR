#' Write meta genes to a file
#' 
#' Write a csv file of meta gene lists from input list
#' @param meta List of meta genes output from \code{meta_genes}
#' @param file Output file name
#' @return \code{NULL}
#' @export
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50, ncol=100, rank=5)
#' s <- scNMFSet(x$x)
#' s <- vb_factorize(s, ranks=seq(2,8), nrun=5)
#' plot(s)
#' m <- meta_genes(s, rank=5)
#' write_meta(m, file='meta.csv')
write_meta <- function(meta,file){
  
  rank <- length(meta)
  nmax <- 0
  for(k in seq_len(rank))
    if(length(meta[[k]])> nmax) nmax <- length(meta[[k]])
  dat <- matrix('',nrow=nmax,ncol=rank)
  for(k in seq_len(rank))
    dat[seq_len(length(meta[[k]])),k] <- meta[[k]]
  rownames(dat) <- seq_len(nmax)
  colnames(dat) <- seq_len(rank)
  utils::write.csv(dat,file=file)
  return(invisible(meta))
}

#' Determine optimal rank
#' 
#' Takes as main argument \code{scNMFSet} object containing factorized output
#' and estimate the optimal rank.
#' 
#' @param object \code{scNMFSet} object containing factorization output, or 
#'        data frame containing the rank-evidence profile.
#' @param df Degrees of freedom for split fit. Upper bound is the total number of
#'        data points (number of rank values scanned).
#' @param BF.threshold Bayes factor threshold for statistical threshold. 
#' @param type \code{c(1,2)}. Type 1 is where there is a clear maximum. Type 2 
#'        is where marginal likelihood reaches a maximal level and stays constant.
#'        If omitted, the type will be inferred from data.
#' @param m Number of features (e.g., genes) in the count matrix. Only necessary when
#'        \code{object} is of type \code{data.frame}.
#' @return List containing \code{type} and \code{ropt} (optimal rank).
#' 
#' @details The input object is used along with Bayes factor threshold to determine the
#'          heterogeneity type (1 or 2) and the optimal rank. 
#'          If evidence(rank 1)/evidence(rank2) > \code{BF.treshold}, rank 1 is favorable than rank 2.
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50, ncol=100, rank=5)
#' s <- scNMFSet(x$x)
#' s <- vb_factorize(s, ranks=seq(2,8), nrun=5)
#' plot(s)
#' optimal_rank(s)
#' @export
optimal_rank <- function(object, df=10, BF.threshold=3, type=NULL, m=NULL){
  
  if(is(object,'scNMFSet')){
    me <- measure(object)[,c(1,2)]
    m <- nrow(object)
  }
  else if(is(object,'data.frame')){
    me <- object[,c(1,2)]
    if(is.null(m)) stop('No. of rows unknown')
  }
  else stop('Inappropriate class of object')
  
  df <- min(df,nrow(me))
  fs <- stats::smooth.spline(x=me[,1],y=me[,2],df=df)
  rst <- fs$x[which.max(fs$y)]     # arg max_r (L)
  bf <- log(BF.threshold)/m
  
  if(is.null(type)){    # determine the type
    rmax <- max(me[,1])
    flag <- abs(fs$y-max(fs$y)) <= bf
    range <- fs$x[flag]  # rank values within the threshold around optimum
    if(rmax %in% range) type <- 2
    else type <- 1
  }
  
  if(type==1) ropt <- rst
  else{
    sl <- slope(fs$y,fs$x)
    if(sum(sl < bf)>0) 
      idx <- which(sl < bf)[1]
    else
      idx <- nrow(me)
    ropt <- fs$x[idx]
  }
  
  return(list(type=type, ropt=ropt))
}

slope <- function(y,x){
  
  n <- length(x)
  s <- rep(0,n)
  s[1] <- (y[2]-y[1])/(x[2]-x[1])
#  for(i in seq(2,n-1))
#    s[i] <- mean(c((y[i+1]-y[i])/(x[i+1]-x[i]),
#                   (y[i]-y[i-1])/(x[i]-x[i-1])))
#  s[n] <- (y[n]-y[n-1])/(x[n]-x[n-1])
  for(i in seq(2,n-1))
    s[i] <- (y[i+1]-y[i])/(x[i+1]-x[i])
  s[n] <- s[n-1]
  
  return(s)
}

#' Meta gene table with CV
#'
#' Generates meta gene table with coefficient of variation 
#' @param object Main object containing factorization outcome
#' @param rank Rank for which meta gene is to be found
#' @param basis.matrix Basis matrix to work with. Only necessary when 
#'        \code{object} is \code{NULL}.
#' @param dbasis Variance of basis matrix. Only necessary when 
#'        \code{object} is \code{NULL}.
#' @param max.per.cluster Maximum meta genes per cluster.
#' @param gene_names Name of genes. If \code{NULL}, will be taken from row names.
#' @param subtract.mean Standardize magnitudes of basis elements by subtracting mean
#' @param log Use geometric mean.
#' @param cv.max Upper bound for CV in selecting meta genes.
#' @return Data frame with meta genes and their CV in each column.
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50, ncol=100, rank=5)
#' s <- scNMFSet(x$x)
#' s <- vb_factorize(s, ranks=seq(2,8), nrun=5)
#' plot(s)
#' meta_gene.cv(s, rank=5)
#' @export
meta_gene.cv <- function(object=NULL, rank, basis.matrix=NULL, 
                         dbasis=NULL, max.per.cluster=100,
                         gene_names=NULL, subtract.mean=TRUE,
                         log=TRUE, cv.max=1){
  
  if(is.null(basis.matrix)){
    idx <- ranks(object)==rank
    w <- basis(object)[idx][[1]]
    cw <- dbasis(object)[idx][[1]]/w
    if(subtract.mean){
      if(log) w <- log10(w)
      w <- w - rowMeans(w)
      if(log) w <- 10^w
    }
  } else{
    w <- basis.matrix
    rank <- ncol(w)
    cw <- dbasis/w
  }
  if(!is.null(gene_names)) rownames(w) <- rownames(cw) <- gene_names
  nmax <- min(max.per.cluster, nrow(w))

  maxrow <- 0  
  for(k in seq_len(rank)){
    x <- data.frame(rep('',nmax), rep(NA,nmax), rep(NA,nmax),
                    stringsAsFactors=FALSE)
    names(x) <- c(paste0('Gene_',k), paste0('W_',k),paste0('CV_',k))
    idx <- order(w[,k],decreasing=TRUE)[seq_len(nmax)]
    v <- w[idx,]
    cv <- cw[idx,]
    sig <- which(cv[,k] <= cv.max)
    nsig <- seq_len(length(sig))
    x[nsig,1] <- rownames(v)[sig]
    x[nsig,2] <- v[sig,k]
    x[nsig,3] <- cv[sig,k]
    if(length(nsig)>maxrow) maxrow <- length(nsig)
    if(k==1) dat <- x
    else dat <- cbind(dat, x)
  }
  dat <- dat[seq_len(maxrow),]

  return(dat)
}
