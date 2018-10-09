#' Write meta genes to a file
#' @export
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
  write.csv(dat,file=file)
}

#' Determine optimal rank
#' @export
optimal_rank <- function(object, df=10, type=1, max.slope=1e-4){
  
  if(class(object)=='scNMFSet')
    me <- measure(object)[,1:2]
  else if(class(object)=='data.frame')
    me <- object[,1:2]
  else stop('Inappropriate class of object')
  
  df <- min(df,nrow(me))
  fs <- smooth.spline(x=me[,1],y=me[,2],df=df)
  if(type==1) rst <- fs$x[which.max(fs$y)]
  else{
    sl <- slope(fs$y,fs$x)
    if(sum(sl < max.slope)>0) 
      idx <- which(sl<max.slope)[1]
    else
      idx <- nrow(me)
    rst <- fs$x[idx]
  }
  
  return(rst)
}

slope <- function(y,x){
  n <- length(x)
  s <- rep(0,n)
  s[1] <- (y[2]-y[1])/(x[2]-x[1])
  for(i in seq(2,n-1))
    s[i] <- mean(c((y[i+1]-y[i])/(x[i+1]-x[i]),
                   (y[i]-y[i-1])/(x[i]-x[i-1])))
  s[n] <- (y[n]-y[n-1])/(x[n]-x[n-1])
  
  return(s)
}

#' Meta gene table with CV
#' @export
meta_gene.cv <- function(object, rank, basis.matrix=NULL, dbasis=NULL,
                       max.per.cluster=100,gene_names=NULL,
                       subtract.mean=TRUE,log=TRUE,
                       scheme='max', cv.min=0.1){
  
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
    sig <- which(cv[,k] <= cv.min)
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