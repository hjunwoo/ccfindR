#' Cell type assignment via GSEA
#'
#' Computes GSEA enrichment score of marker sets in meta gene list
#' 
#' If \code{obj} is of clas \code{scNMFSet}, it computes meta gene list using
#' \code{\link{meta_gene.cv}}. Otherwise, \code{obj} is expected to be 
#' a data frame of the same structure as the output of \code{meta_gene.cv};
#' the number of rows same as the total number of metagenes per cluster, 
#' three columns per each cluster (gene name, meta-gene score, and coefficient of variation).
#' The argument \code{gset} is a list of gene sets to be checked for enrichment in 
#' each cluster meta gene list. The enrichment score is computed using 
#' the GSEA algorithm \insertCite{subramanian_etal}{ccfindR}.
#'
#' @param obj Object of class \code{scNMFSet}.
#' @param rank Rank to examine
#' @param gset List of gene sets to be used as markers
#' @param gene_names Names of genes to be used for meta-gene identification
#' @param p Enrichment score exponent.
#' @param remove.na Remove gene sets with no overlap
#' @param p.value Estimatte p values using permutation
#' @param nperm No. of permutation replicates
#' @param progress.bar Display progress bar for p value computation
#' @param grp.prefix Gene name prefix to search for with wildcard matches in query
#' @return Matrix of enrichment score statistics with cell types in rows and
#'         clusters in columns
#' @references
#' \insertRef{subramanian_etal}{ccfindR}
#' @examples
#' dir <- system.file('extdata',package='ccfindR')
#' pbmc <- read_10x(dir)
#' pbmc <- vb_factorize(pbmc, ranks=5)
#' meta <- meta_gene.cv(pbmc,rank=5, gene_names=rowData(pbmc)[,2])
#' markers <- list('B cell'=c('CD74','IG','HLA'),
#'                 'CD8+ T'=c('CD8A','CD8B','GZMK','CCR7','LTB'),
#'                 'CD4+ T'=c('CD3D','CD3E','IL7R','LEF1'),
#'                 'NK'=c('GNLY','NKG7','GZMA','GZMH'),
#'                 'Macrophage'=c('S100A8','S100A9','CD14','LYZ','CFD'))
#' gsea <- assignCelltype(meta, rank=5, gset=markers, grp.prefix=c('IG','HLA'))
#' gsea
#' @export
assignCelltype <- function(obj, rank, gset, gene_names=NULL, p=0, remove.na=FALSE,
                 p.value=FALSE, nperm=1000, progress.bar=TRUE,
                 grp.prefix=c('IG')){
  
  if(is(obj, 'scNMFSet'))
    meta <- meta_gene.cv(object=obj, rank=rank, max.per.cluster=Inf,
                       gene_names=gene_names, subtract.mean=TRUE,
                       log=TRUE, cv.max=Inf)
  else if(is(obj, 'data.frame')){
    meta <- obj
    if(NCOL(meta)!=rank*3) stop('Incorrect dimension of meta')
  } else stop('Incorrect input type of obj')
    
  glist <- meta[,seq(1,3*rank,by=3)]
  gwgt <- meta[,seq(2,3*rank+1,by=3)]
  
  ES <- gsea(glist=glist, gwgt=gwgt, gset=gset, p=p, remove.na=remove.na,
             grp.prefix=grp.prefix)
  
  if(!p.value) return(ES)
  N <- NROW(glist)
  nS <- length(gset)
  Ep <- matrix(0, nrow=NROW(ES), ncol=rank)
  if(progress.bar) pb <- txtProgressBar(style=3)
  for(k in seq_len(nperm)){
    perm <- sample(N)
    x <- gsea(glist=glist[perm,], gwgt=gwgt[perm,], gset=gset, p=p, 
              remove.na=remove.na, grp.prefix=grp.prefix)
    Ep <- Ep + (ES < x)
    if(progress.bar & k %% 10==0)
      setTxtProgressBar(pb, value=k/nperm)
  }
  if(progress.bar) close(pb)
  Ep <- Ep / nperm
  
  return(list(ES=ES, pvalue=Ep))
}

gsea <- function(glist, gwgt, gset, p, remove.na, grp.prefix=c('IG','KRT')){
  
  rank <- NCOL(glist)
  nS <- length(gset)
  N <- NROW(glist)
  
  ES <- matrix(0, nrow=nS, ncol=rank)
  colnames(ES) <- seq_len(rank)
  rownames(ES) <- names(markers)
  for(k in seq_len(rank)){
    for(is in seq_len(nS)){
      gs <- gset[[is]]
      gl <- glist[,k]
      flag <- !is.na(gwgt[,k])
      gl <- gl[flag]
      gw <- gwgt[flag,k]
      x <- overlap(query=gl, glist=gs, grp.prefix=grp.prefix)
      if(sum(x)==0){
        ES[is,k] <- NA
        next()
      }
      ph <- gw^p
      ph <- cumsum(x * ph)
      phit <- ph / ph[length(ph)]
      
      y <- ! (gl %in% gs)
      pm <- cumsum(y)
      pmiss <- pm / pm[length(pm)]

      ES[is,k] <- max(phit-pmiss)
    }
  }
  
  if(remove.na) ES <- ES[!is.na(ES[,1]),]
  
  return(ES)
}

overlap <- function(query, glist, grp.prefix=c('IG')){
  
  glist0 <- glist[!glist %in% grp.prefix] # non-IG genes
  x <- query %in% glist0
  grp <- grp.prefix[grp.prefix %in% glist]
  for(i in seq_along(grp)){
    gr <- grp[i]
    x1 <- substr(query, start=1, stop=nchar(gr)) %in% gr
    x <- x | x1
  }
  return(x)
}
