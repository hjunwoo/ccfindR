#' Plot signatures
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=96,ncol=100,rank=3)
#' mut <- mut_list(3)
#' rownames(x$x) <- mut[1:nrow(x$x)]
#' s <- scNMFSet(count=x$x)
#' s <- vb_factorize(s,ranks=2:5,normalize.signature=TRUE)
#' plot_signatures(s,rank=3)
#' 
#' @export
plot_signatures <- function(obj,dobj=NULL,rank, 
                            col=NULL, kmer.size=3, pmax=NULL,
                            cosmic=NULL, similarity.cutoff=0.5,
                            mfrow=c(rank,1),show.sd=TRUE){

  old.par <- graphics::par(no.readonly=TRUE)
  # read cosmic signature list
  if(is.null(cosmic))
    cosmic <- system.file('extdata','PCAWG_SBS_v3.txt',
                          package='ccfindR')
  cosmic <- read.table(cosmic,header=TRUE,sep=' ')
  # rownames = mutation context (3-mer), colnames=signature #

  if(class(obj)=='scNMFSet'){
    id <- which(ranks(obj)==rank)
    P <- basis(obj)[[id]]    # signatures
    dP <- dbasis(obj)[[id]]
    E <- coeff(obj)[[id]]    # mutation loads
  }
  else if(class(obj)=='matrix'){
    P <- obj
    rank <- ncol(P)
    dP <- dobj
  }
  
  graphics::par(mfrow=mfrow,lwd=0.1)
  
  cosim <- cos_sim(P,cosmic,kmer.size) 
                  # cosine similarity to cosmic signatures
  
  kd <- (kmer.size-1)/2
# if(is.null(pmax)) pmax <- max(P)
  
  if(is.null(col)) col <- c('blue','gray60','red','orange','green','cornflowerblue')
  x <- rownames(P)
  central <- vapply(x,function(x){substr(x,start=kd+2,stop=kd+4)},
                    character(1))
  fl <- levels(factor(central))
  xlabel <- vapply(fl, function(x){mean(which(central==x))},numeric(1))  
  # position of mutation label
  col2 <- c()
  for(i in 1:6) col2 <- c(col2, rep(col[i],table(central)[fl][i]))
  
  for(k in seq_len(rank)){
    if(is.null(pmax)) pmaxk <- max(P[,k])*1.2
    else pmaxk <- pmax
    bp <- barplot(P[,k], names.arg='', col=col2, border=NA,
            ylab='Probability',ylim=c(0,pmaxk),las=1)
    if(show.sd){
      segments(x0=bp,x1=bp,y0=P[,k]-dP[,k],y1=P[,k]+dP[,k],col='black',lty=1,lwd=0.5)
      dx <- (bp[2]-bp[1])*0.3
      segments(x0=bp-dx,x1=bp+dx,y0=P[,k]-dP[,k],y1=P[,k]-dP[,k],col='black',lty=1,lwd=0.5)
      segments(x0=bp-dx,x1=bp+dx,y0=P[,k]+dP[,k],y1=P[,k]+dP[,k],col='black',lty=1,lwd=0.5)
    }  
    sid <- colnames(cosmic)[which.max(cosim[k,])]
    cos <- max(cosim[k,])
    if(cos>similarity.cutoff)
      title(adj=0, main=paste0(
        'Process ', k,': ',sid,'-like (',signif(cos,3),')'),cex.main=1)
    else
      title(adj=0, main=paste0('Process ',k),cex.main=1)
    text(x=xlabel*1.2, label=fl, y=-pmaxk*0.15, xpd=NA, col=col)
  }
  graphics::par(old.par)
  invisible()
}

#' Generate vector of mutation context names
#' @export
mut_list <- function(kmer.size=3){
  
  nt <- c('A','C','G','T')
  dk <- (kmer.size-1)/2
  stype <- c('[C>A]','[C>G]','[C>T]','[T>A]','[T>C]','[T>G]')
  z <- vector('list',2*dk+1)
  for(k in seq_len(2*dk)) z[[k]] <- nt
  z[[2*dk+1]] <- stype
  x <- expand.grid(z)
  x <- x[,c(seq(dk+1,2*dk),2*dk+1,seq(1,dk))]
  subtypes <- apply(x,1,paste,collapse='')

  return(subtypes)
}

#' Compute cosine similarity scores
#' @export
cos_sim <- function(X,Y,kmer.size=NULL){
  
# if(nrow(X)!=nrow(Y))  # reduce 5-mer into 3-mer
  if(!is.null(kmer.size)) if(kmer.size==5)
    X <- penta2trimer(X)
  
  if(nrow(X)!=nrow(Y)){
    if(is.null(rownames(X)) | is.null(rownames(Y)))
      stop('X and Y not comparable in cos_sim')
    Y <- Y[rownames(Y)%in% rownames(X),]
  }
  if(!is.null(rownames(X)) & !is.null(rownames(Y)))
    Y <- Y[match(rownames(X),rownames(Y)),]
  m <- nrow(X)
  nx <- ncol(X)
  ny <- ncol(Y)
  cr <- matrix(0, nrow=nx, ncol=ny)
  colnames(cr) <- colnames(Y)
  rownames(cr) <- colnames(X)
  for(i in seq(1,nx)) for(j in seq(1,ny))
    cr[i,j] <- sum(X[,i]*Y[,j])/sqrt(sum(X[,i]^2)*sum(Y[,j]^2))
  
  return(cr)
}

#' Plot signature overlap with database
#' @export
signature_map <- function(obj, rank, cosmic=NULL, col=NULL, mar=c(4,4,4,7),
                          cex.axis=0.8){
  
  old.par <- graphics::par(no.readonly=TRUE)
  graphics::par(mar=mar,cex.axis=cex.axis)
  # read cosmic signature list
  if(is.null(col))
    col <- RColorBrewer::brewer.pal(n=9,'YlOrRd')
  if(is.null(cosmic))
    cosmic <- system.file('extdata','signatures_probabilities.txt',
                          package='ccfindR')
  cosmic <- read.table(cosmic,header=TRUE,sep='\t')
  x <- cosmic$Somatic.Mutation.Type
  cosmic <- cosmic[,seq(4,33)]
  rownames(cosmic) <- x
  
  id <- which(ranks(obj)==rank)
  P <- basis(obj)[[id]]    # signatures
  cosim <- cos_sim(P,cosmic) # cosine similarity to cosmic signatures
  
  nrow <- nrow(cosim)
  ncol <- ncol(cosim)
  graphics::image(seq_len(nrow),seq_len(ncol),cosim,col=col, xlab='Process',
        ylab='Signatures',las=1)
  
  ym <- ncol(cosmic)/7
  x0 <- rank+1
  x1 <- rank+1+rank/10
  y0 <- seq(1,ym,length.out=length(col))
  dy <- y0[2]-y0[1]
  rect(xleft=x0,xright=x1,ybottom=y0,ytop=y0+dy,col=col,border=NA,xpd=NA)
  text(x=rank+1+(rank+2)/10,y=1,label=round(min(cosim),1),xpd=NA,cex=0.7,adj=0)
  text(x=rank+1+(rank+2)/10,y=ym,label=round(max(cosim),1),xpd=NA,cex=0.7,adj=0)
  
  graphics::par(old.par)  # restore graphical parameters
  invisible()
}

# reduce 5-mer signature into 3-mer
penta2trimer <- function(X){
  
  central <- vapply(rownames(X),function(x){substr(x,start=2,stop=8)},
                    character(1))
  X <- aggregate(X,by=list(match(central,central[!duplicated(central)])),
                 FUN=sum)
  X <- X[,-1]
  rownames(X) <- central[!duplicated(central)]
  
  return(X)
}

#' Plot mutation load
#' @export
plot_mutload <- function(obj,rank, sample.id=NULL, col=NULL,labRow=''){
  
  if(is.null(col))
    col <- RColorBrewer::brewer.pal(n=9,'YlOrRd')
  
  id <- which(ranks(obj)==rank)
  H <- coeff(obj)[[id]]    # signatures
  if(!is.null(sample.id)) 
    H <- H[,sample.id]
  stats::heatmap(log(t(H)),revC=TRUE,scale='none',col=col,Rowv=NULL,Colv=NA,
                 labRow=labRow)
  invisible()
}

#' @export
top_genes <- function(object,basis.matrix=NULL,rank,ntop=5,gene_names=NULL){
  
  if(is.null(basis.matrix))
    w <- basis(object)[[which(ranks(object)==rank)]]
  else
    w <- basis.matrix
  if(is.null(gene_names)) gene_names <- rownames(w)
  
  for(k in seq_len(rank)){
    v <- w[,k]
    g <- order(v,decreasing=TRUE)[seq_len(ntop)]
    if(k==1) dat <- data.frame(gene_names[g])
    else dat <- cbind(dat,data.frame(gene_names[g]))
  }
  names(dat) <- seq_len(rank)
  
  return(dat)
  
}