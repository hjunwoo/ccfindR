#' Read 10x data and generate scNMF object
#' 
#' Read count, gene, and barcode annotation data in 10x format 
#' and create an object of class \code{scNMFSet}.
#' 
#' Files for \code{count}, \code{genes}, and \code{barcodes} are
#' assumed to be present in \code{dir}. 
#' Count data are in sparse "Matrix Market" format 
#' (\url{https://math.nist.gov/MatrixMarket/formats.html}).
#' 
#' @param dir Name of directory containing data files.
#' @param count Name of count matrix file.
#' @param genes Name of gene annotation file.
#' @param barcodes Name of cell annotation file.
#' @param remove.zeros If \code{TRUE}, empty rows/columns are 
#'        removed.
#' @examples
#' library(S4Vectors)
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' rowData(s) <- DataFrame(seq_len(4))
#' colData(s) <- DataFrame(seq_len(3))
#' write_10x(s,dir='.')
#' s <- read_10x(dir='.')
#' s
#' @return Object of class \code{scNMFSet}
#' @export
#' @import SingleCellExperiment
read_10x <- function(dir, count='matrix.mtx', genes='genes.tsv',
                     barcodes='barcodes.tsv', remove.zeros=TRUE){
  
  if(!dir.exists(dir)) stop(cat('Input directory',dir,'does not exist\n'))
  count <- paste0(dir, '/', count)
  if(!file.exists(count)) stop(cat('Count file',count,'does not exist\n'))
  Mat <- as(Matrix::readMM(count),'dgCMatrix')
  genes <- paste0(dir, '/', genes)
  if(!file.exists(genes)) 
    stop(cat('Count file',genes,'does not exist\n'))
  glist <- utils::read.table(genes,stringsAsFactors=FALSE)
  barcodes <- paste0(dir,'/', barcodes)
  if(!file.exists(barcodes)) 
    stop(cat('Count file',barcodes,'does not exist\n'))
  clist <- utils::read.table(barcodes,stringsAsFactors=FALSE)
  rownames(Mat) <- rownames(glist) <- glist[,1]
  colnames(Mat) <- rownames(clist) <- clist[,1]
  
  x <- scNMFSet(assays=list(counts=Mat), remove.zeros=FALSE)
  rowData(x) <- as(glist, 'DataFrame')
  rownames(x) <- rownames(glist)
  clist <- as(clist, 'DataFrame')
  colData(x) <- clist
  colnames(x) <- rownames(clist)
  if(remove.zeros) x <- remove_zeros(x)
  return(x)
}

#' Filter cells with quality control criteria
#' 
#' Remove low quality cell entries from object
#' 
#' Takes as input \code{scNMFSet} object and plots histogram of UMI 
#' counts for each 
#' cell. Optionally, cells are filtered using minimum and maximum UMI counts. 
#' The resulting object is returned after removing empty rows and columns, 
#' if any.
#' 
#' @param object \code{scNMFSet} object
#' @param umi.min Minimum UMI count for cell filtering
#' @param umi.max Maximum UMI count for cell filtering
#' @param plot If \code{TRUE}, the UMI count distribution of all cells 
#'   will be displayed. Cells selected are colored red.
#' @param remove.zeros Remove rows/columns containing zeros only
#' @return \code{scNMFSet} object with cells filtered.
#' @examples
#' set.seed(1)
#' s <- scNMFSet(matrix(stats::rpois(n=1200,lambda=3),40,30))
#' s <- filter_cells(s,umi.min=10^2.0,umi.max=10^2.1)
#' @export
filter_cells <- function(object, umi.min=0, umi.max=Inf, plot=TRUE, 
                         remove.zeros=TRUE){
  
  umi.count <- Matrix::colSums(counts(object))
  selected_cells <- umi.min <= umi.count & umi.count <= umi.max
  
  h <- graphics::hist(log10(umi.count), plot=FALSE)
  n <- length(h$counts)
  color <- rep('white',n)
  if(umi.min | umi.max<Inf)
    color[log10(umi.min) <= h$mids & h$mids <= log10(umi.max)] <- 'red'
  if(plot) plot(h, col=color, main='Cells')
  if(umi.min | umi.max<Inf)
    cat(sum(selected_cells),'cells out of',dim(object)[2],'selected\n')
  
  object <- object[, selected_cells]
  return(remove_zeros(object))
}

#' Filter genes with quality control criteria
#' 
#' Select genes with high relative variance in count data for 
#' further analysis
#' 
#' Takes as input \code{scNMFSet} object and scatterplot no. of cells 
#' expressed versus VMR (variance-to-mean ratio)
#' for each gene. Optionally, genes are filtered using minimum VMR 
#' together with a range of no. of cells expressed.
#'
#' @param object \code{scNMFSet} object.
#' @param markers A vector containing marker genes to be selected. 
#'    All rows in 
#' \code{rowData} that contain columns matching this set will 
#'    be selected. 
#' @param vmr.min Minimum variance-to-mean ratio for gene filtering.
#' @param min.cells.expressed Minimum no. of cells expressed for gene 
#'    filtering.
#' @param max.cells.expressed Maximum no. of cells expressed for gene 
#'    filtering.
#' @param rescue.genes Selected additional genes whose (non-zero)
#'        count distributions have at least one mode. 
#' @param plot Plot the distribution of no. of cells expressed vs. VMR.        
#' @param log Axis in log-scale, \code{c('x','y','xy')}.
#' @param progress.bar Display progress of mode-gene scan or VMR
#'    calculation with \code{save.memory = TRUE}.
#' @param save.memory For a very large number of cells, calculate VMR 
#'    row by row while avoiding calls to \code{as.matrix()}.
#'    Progress bar will be displayed unless \code{progress.bar=FALSE}.
#' @param cex Symbol size for each gene in the plot. 
#' @return Object of class \code{scNMFSet}.
#' @examples
#' set.seed(1)
#' s <- scNMFSet(matrix(stats::rpois(n=1200,lambda=3),40,30))
#' s <- filter_genes(s,vmr.min=1.0,min.cells.expressed=28,
#'         rescue.genes=FALSE)
#' @export
filter_genes <- function(object, markers= NULL, vmr.min=0, 
                         min.cells.expressed=0, max.cells.expressed=Inf, 
                         rescue.genes=FALSE, progress.bar=TRUE,
                         save.memory=FALSE, plot=TRUE, log='xy',cex=0.5){

  ncexpr <- Matrix::rowSums(counts(object) > 0) 
  # no. of cells expressing each gene
  count <- counts(object)
  count <- count[ncexpr>0,]
  ncexpr <- ncexpr[ncexpr>0]
  ngenes <- nrow(count)
  
  selected_genes <- variable_genes <- rep(FALSE, ngenes)
  
  if(is.null(markers)) marker_genes=NULL
  else{
    for(k in seq_len(ncol(rowData(object))))
      selected_genes <- selected_genes | 
        (mk <- rowData(object)[,k] %in% markers)
    marker_genes <- selected_genes
  } 
    
  vmr <- calc_vmr(count, save.memory=save.memory, 
                  progress.bar=progress.bar)
  
  variable_genes <- vmr.min < vmr & 
            min.cells.expressed <= ncexpr & 
            ncexpr <= max.cells.expressed
  if(rescue.genes & sum(variable_genes) < ngenes){
    mode_genes <- c()
    cat('Looking for genes with modes ...\n')
    if(progress.bar) pb <- utils::txtProgressBar(style=3)
    for(i in seq_len(ngenes)){
      if(variable_genes[i]) mode <- FALSE
      else mode <- has_mode(counts(object)[i,])
      mode_genes <- c(mode_genes,mode)
      if(progress.bar) utils::setTxtProgressBar(pb,i/ngenes)
    }
    if(progress.bar) close(pb)
    selected_genes <- selected_genes | variable_genes | mode_genes
  }
  else{ 
    selected_genes <- selected_genes | variable_genes
    mode_genes <- NULL
  }
  
  if(!is.null(markers)) if(sum(marker_genes) > 0)
    cat(sum(marker_genes),'marker genes found\n')
  if(vmr.min>0 | min.cells.expressed>0 | max.cells.expressed <Inf){
    cat(sum(variable_genes),'variable genes out of',dim(object)[1],'\n')
    if(rescue.genes) cat(sum(selected_genes & ! variable_genes),
                         'additional genes rescued\n')
    cat(sum(selected_genes),'genes selected\n')
  }
  
  if(plot) plot_genes(object, vmr=vmr, ncexpr=ncexpr, 
    selected_genes=selected_genes, variable_genes=variable_genes,
    mode_genes=mode_genes, marker_genes=marker_genes, log=log,cex=cex) 
  
  return(object[selected_genes,])
}

# Calculate variance-to-mean ratio
calc_vmr <- function(count, save.memory=FALSE, progress.bar=TRUE){
  
  gmean <- Matrix::rowMeans(count)
  ngenes <- nrow(count)
  if(!save.memory) 
    var <- rowVars(count,means=gmean) # memory intensive
  else{
    if(progress.bar){
      cat('Calculating vmr...\n')
      pb <- utils::txtProgressBar(style=3)
    }
    var <- c()
    denom <- ncol(count)  
    for(i in seq_len(nrow(count))){
      g <- count[i,]-gmean[i]
      var <- c(var, sum(g*g)/denom)
      if(progress.bar) utils::setTxtProgressBar(pb, i/ngenes)
    }
    if(progress.bar) close(pb)
  }
  return(var/gmean)
}

#' Plot gene variance distributions
#' 
#' Gene variance to mean ratio and the number of expressing cells are
#' plotted.
#' 
#' This function can be called separately or is also called within
#' \code{\link{filter_genes}} by default. 
#' In the latter case, parameters other
#' than \code{object} will have been already filled. If called separately
#' with \code{NULL} gene lists, VMR is recalculated but gene selection
#' is not done.
#' 
#' @param object Object containing count data
#' @param vmr Variance to mean ratio (VMR)
#' @param ncexpr Number of cells expressing each gene
#' @param selected_genes Logical vector specifing genes selected
#' @param variable_genes Logical vector specifing genes with high VMR
#' @param mode_genes Logical vector specifying genes with nonzero modes
#' @param marker_genes Logical vector specifying marker genes
#' @param save.memory If \code{TRUE}, calculate VMR using slower
#'        method to save memory. Not used when gene lists are supplied.
#' @param progress.bar Display progress bar for VMR calculation.
#'        Not used when gene lists are supplied.
#' @param cex Symbol size for genes (supplied to \code{plot()}).
#' @param log Axis in log-scale, \code{c('x','y','xy')}.
#' @return \code{NULL}
#' @examples
#' set.seed(1)
#' s <- scNMFSet(matrix(stats::rpois(n=1200,lambda=3),40,30))
#' plot_genes(s)
#' @export
plot_genes <- function(object, vmr=NULL, ncexpr=NULL, selected_genes=NULL, 
                  variable_genes=NULL, mode_genes=NULL, 
                  marker_genes=NULL, save.memory=FALSE, progress.bar=TRUE, 
                  log='xy', cex=0.5){
  
  if(is.null(ncexpr) | is.null(vmr)){  # no. of cells expressing each gene
    ncexpr <- Matrix::rowSums(counts(object) > 0) 
    count <- counts(object)
    count <- count[ncexpr>0,]
    ncexpr <- ncexpr[ncexpr>0]
  }
  if(is.null(vmr))
    vmr <- calc_vmr(count, save.memory=save.memory, progress.bar=progress.bar)
  
  ncexpr <- ncexpr[vmr>0]
  vmr <- vmr[vmr>0]
  ngenes <- length(vmr)
  if(is.null(selected_genes)) selected_genes <- rep(FALSE, ngenes)
  
  xlim=c(min(ncexpr),max(ncexpr))
  if(sum(selected_genes) < ngenes | !is.null(variable_genes)){
    graphics::plot(x=ncexpr[!selected_genes], y=vmr[!selected_genes], 
         xlim=xlim,ylim=c(min(vmr),max(vmr)), log=log, pch=21, 
        col='gray', bg='white', xlab='No. of cells expressed', ylab='VMR', 
        main='Genes', cex=cex, bty='n')
    if(!is.null(variable_genes))
      graphics::points(x=ncexpr[variable_genes], y=vmr[variable_genes], 
                       pch=21, bg='red', cex=cex, lwd=0.5)
    if(!is.null(mode_genes))
        graphics::points(x=ncexpr[mode_genes], y=vmr[mode_genes], 
              pch=21, bg='blue', cex=cex, lwd=0.5)
    if(!is.null(marker_genes))
      if(sum(marker_genes) > 0)
        points(x=ncexpr[marker_genes], y=vmr[marker_genes], pch=21, 
                 bg='orange', cex=cex, lwd=0.5)
  }else{
    graphics::plot(x=ncexpr[variable_genes], y=vmr[variable_genes], 
        xlim=xlim,ylim=c(min(vmr),max(vmr)), pch=21, bg='white', 
        cex=cex, lwd=0.5, log=log, col='gray', xlab='No. of cells expressed', 
        ylab='VMR', main='Genes', bty='n')
    if(!is.null(marker_genes)) if(sum(marker_genes) > 0)
        graphics::points(x=ncexpr[marker_genes], y=vmr[marker_genes], pch=21, 
               bg='orange', cex=cex, lwd=0.5)
  }
  return(invisible(object))
}

#' Normalize count data
#' 
#' Rescale count matrix entries such that all cells have the same library size.
#' 
#' For analysis purposes, it is sometimes useful to rescale integer count
#' data into floats such that all cells have the same median counts. 
#' This function will calculate the median of all UMI counts of cells (total
#' number of RNAs derived from each cell). All count data are then rescaled
#' such that cells have uniform UMI count equal to the median. 
#' 
#' @param object \code{scNMFSet} object.
#' @return \code{scNMFSet} object with normalized count data.
#' @examples
#' library(Matrix)
#' set.seed(1)
#' s <- scNMFSet(count=matrix(rpois(n=1200,lambda=3),40,30))
#' colMeans(counts(s))
#' s <- normalize_count(s)
#' colMeans(counts(s))
#' @export
normalize_count <- function(object){
  
  count <- counts(object)
  umi.count <- Matrix::colSums(count)
  count <- t(t(count)/umi.count)
  med <- stats::median(umi.count)
  count <- count*med
  counts(object) <- count
  object
}

has_mode <- function(g){
  
  tb <- table(g)
  flag <- FALSE
  tb <- tb[names(tb)]    # remove zero count
  n <- length(tb)
  if(n<2) return(FALSE)
  for(k in seq_len(n-1)) if(tb[k] < tb[k+1]) return(TRUE)
  return(FALSE)

}

rowVars <- function(x, means){
  if(missing(means)) means <- Matrix::rowMeans(x)
  Matrix::rowSums((x-means)^2)/(dim(x)[2]-1)
}

#' Plot heatmap of metagene matrix
#' 
#' Generate heatmap of metagenes derived from factorization of count data.
#' 
#' Wrapper for \code{\link[stats]{heatmap}} to display metagenes and
#' associated basis matrix element magnitudes. Factorization results inside
#' an object specified by its rank value will be retrieved, and metagene
#' sets identified from clusters.
#' 
#' @param object Object of class \code{scNMFSet}.
#' @param rank Rank value for which the gene map is to be displayed. 
#'   The object must contain the corresponding slot (one element of 
#'        \code{basis(object)[[k]]} for which \code{ranks(object)[[k]]==rank}.
#' @param markers Vector of gene names containing markers to be included 
#'         in addition to the metagenes. All entries of \code{rowData(object)}
#'         matching them will be added to the metagene list.
#' @param subtract.mean Process each rows of basis matrix \code{W} by 
#'        standardization using the mean of elements within the row.
#' @param log If \code{TRUE}, \code{subtract.mean} uses geometric mean
#'        and division. Otherwise, use arithmetic mean and subtraction.
#' @param max.per.cluster Maximum number of metagenes per cluster.
#' @param Colv \code{NA} suppresses reordering and dendrogram of clusters
#'        along the column. See \code{\link[stats]{heatmap}}.
#' @param col Colors for the cluster panels on the left and top.
#' @param gene.names Names to be used in the plot for genes.
#' @param main Title of plot.
#' @param ... Other arguments to be passed to \code{\link{heatmap}}, 
#'            \code{\link{image}}, and \code{\link{plot}}.
#'            
#' @details If \code{object} contains multiple ranks, only the requested 
#'   rank's basis matrix W will be displayed. The genes displayed in rows
#'      are selected by "max" scheme
#'     [Carmona-Saez, BMC Bioinformatics (2006), 
#'    \url{https://doi.org/10.1186/1471-2105-7-54}]:
#'      for each cluster (\code{k in 1:ncol}), 
#'      rows of W are sorted by decreasing order
#'      of \code{W[,k]}. Marker genes for \code{k} are those among the top
#'      \code{nmarker} for which \code{W[,k]} is maximum within each row.
#' 
#' @return \code{NULL}
#' @examples  
#' set.seed(1)
#' x <- simulate_data(nfeatures=10,nsamples=c(20,20,60))
#' rownames(x) <- seq_len(10)
#' colnames(x) <- seq_len(100)
#' s <- scNMFSet(count=x,rowData=seq_len(10), colData=seq_len(100))
#' s <- vb_factorize(s,ranks=seq(2,5))
#' plot(s)
#' gene_map(s, rank=3)
#' @export
gene_map <- function(object, rank, markers=NULL, subtract.mean=TRUE,
                     log=TRUE, max.per.cluster = 10, Colv=NA,
                     gene.names=NULL, main='Genes', col=NULL, ...){
  
  if(missing(rank)) rank <- ranks(object)[1] # by default the first rank
  w <- basis(object)[ranks(object)==rank][[1]]
  if(subtract.mean){
    if(log) w <- log10(w)
    w <- w - rowMeans(w)
    if(log) w <- 10^w
  }
  colnames(w) <- seq_len(rank)
  if(!is.null(gene.names)) rownames(w) <- gene.names
  if(dim(w)[1] <= max.per.cluster) select <- rownames(w)
  else select <- gene_select(w, markers, max.per.cluster= max.per.cluster)
  w <- w[select,]
  if(is.null(col)) ccol <- grDevices::rainbow(n = dim(w)[2])
  else ccol <- col
  gid <- apply(w,1,which.max)
  stats::heatmap(w, Colv = Colv, RowSideColors = ccol[gid], 
                 ColSideColors = ccol, revC = TRUE, main =main, 
                 col = RColorBrewer::brewer.pal(n=9,'YlOrRd'),...)
}

#' Plot heatmap of basis matrix
#' 
#' Generate heatmap of features derived from factorization of count data.
#' 
#' This function uses \code{image()} and is more flexible than 
#' \code{gene_map}.
#' @param object Object of class \code{scNMFSet}.
#' @param basis.matrix Basis matrix can be supplied instead of \code{object}.
#' @param rank Rank value for which the gene map is to be displayed. 
#'   The object must contain the corresponding slot (one element of 
#'        \code{basis(object)[[k]]} for which \code{ranks(object)[[k]]==rank}.
#' @param markers Vector of gene names containing markers to be included 
#'         in addition to the metagenes. All entries of \code{rowData(object)}
#'         matching them will be added to the metagene list.
#' @param subtract.mean Process each rows of basis matrix \code{W} by 
#'        standardization using the mean of elements within the row.
#' @param log If \code{TRUE}, \code{subtract.mean} uses geometric mean
#'        and division. Otherwise, use arithmetic mean and subtraction.
#' @param max.per.cluster Maximum number of metagenes per cluster.
#' @param feature.names Names to be used in the plot for features.
#' @param perm Permutation of cluster IDs.
#' @param main Main title.
#' @param cscale Colors for heatmap.
#' @param cex.cluster Cluster ID label size. 
#' @param cex.feature Feature ID label size.
#' @param mar Margins for \code{graphics::par}.
#' @param ... Other arguments to be passed to \code{\link{image}}, and 
#'        \code{\link{plot}}.
#'            
#' @details If \code{object} contains multiple ranks, only the requested 
#'   rank's basis matrix W will be displayed. As in \code{gene_map}, the features
#'   displayed in rows are selected by "max" scheme
#' @return \code{NULL}
#' @examples  
#' set.seed(1)
#' x <- simulate_data(nfeatures=10,nsamples=c(20,20,60))
#' rownames(x) <- seq_len(10)
#' @examples  
#' set.seed(1)
#' x <- simulate_data(nfeatures=10,nsamples=c(20,20,60))
#' rownames(x) <- seq_len(10)
#' colnames(x) <- seq_len(100)
#' s <- scNMFSet(count=x,rowData=seq_len(10), colData=seq_len(100))
#' s <- vb_factorize(s,ranks=seq(2,5))
#' plot(s)
#' feature_map(s, rank=3)
#' @export
feature_map <- function(object, basis.matrix=NULL, rank, markers=NULL, 
                        subtract.mean=TRUE, log=TRUE, 
                        max.per.cluster = 10, feature.names=NULL, perm=NULL,
                        main='Feature map', cscale=NULL, 
                        cex.cluster=1, cex.feature=0.5, mar=NULL, ...){
  
  if(is.null(cscale)) cscale <- RColorBrewer::brewer.pal(n=9,'YlOrRd')
  if(missing(rank)){
     if(!is.null(basis.matrix)) rank <- ncol(basis.matrix)
     else rank <- ranks(object)[1] # by default the first rank
  }
  if(is.null(perm)) perm <- seq_len(rank)

  if(is.null(basis.matrix)){
    w <- basis(object)[ranks(object)==rank][[1]][,perm]
    meta <- meta_genes(object, rank=rank, subtract.mean=subtract.mean,
                       gene_names=feature.names,
                       log=log, max.per.cluster=max.per.cluster)
    colnames(w) <- seq_len(rank)
  }
  else{
    w <- basis.matrix
    meta <- meta_genes(basis.matrix=w, rank=rank, subtract.mean=subtract.mean,
                       gene_names=feature.names,
                       log=log, max.per.cluster=max.per.cluster)
  }

  if(!is.null(feature.names)) rownames(w) <- feature.names
#  if(subtract.mean){
    if(log) w <- w/exp(rowMeans(log(w)))
    else w <- w - rowMeans(w)
#  }

  meta <- meta[perm]
  if(!is.null(markers)){
    markers <- markers[markers %in% rownames(w)]
    markers <- markers[!(markers %in% unlist(meta))] # markers in meta-gene list
  }
  gid <- apply(w[markers,],1, which.max)
  
  w1 <- w[c(unlist(meta),markers),seq_len(rank)]
  idx <- c()
  step <- c()
  for(k in seq_len(rank)){
    tmp <- meta[[k]]
    tmp <- c(tmp, markers[gid==k])
    step <- c(step, length(tmp))
    idx <- c(idx,tmp)
  }
  w1 <- w1[idx,]
  nc <- ncol(w1)
  nr <- nrow(w1)
  
  x <- sweep(w1, 1, rowMeans(w1), check.margin=TRUE)
  sx <- apply(x,1,stats::sd)
  x <- sweep(x, 1, sx, '/', check.margin=FALSE)
  
  if(is.null(mar)) mar <- c(5.1,4.1,4.1,4)
  par(mar=mar)
  graphics::image(seq_len(nc), seq_len(nr), t(x)[,seq(nr,1)], xlim=0.5+c(0,nc),
        ylim=0.5+c(0,nr), axes=FALSE, xlab='',ylab='',col=cscale)
  graphics::axis(1,seq_len(nc), labels=seq_len(nc), las=1, line=-1.0,tick=0, 
       cex.axis=cex.cluster)
  col <- rep('black',nr)
  graphics::text(x=rank+0.7,y=seq(nr,1), labels=rownames(w1), col=col,xpd=NA, 
       cex=cex.feature, adj=0)
  y <- nrow(w1) + 0.5
  graphics::segments(x0=0.5, x1=rank+2, lty=2, y0=y, y1=y, xpd=NA, lwd=0.5)
  for(k in seq_len(rank)){
    graphics::text(x=-0.1, y=y-1, label=k, cex=cex.cluster, xpd=NA)
    y <- y-step[k]
    graphics::segments(x0=0.5, x1=rank+2, lty=2, y0=y, y1=y, xpd=NA, lwd=0.5)
  }
  graphics::title(adj=0.5,main=main)
}
#' Plot heatmap of clustering coefficient matrix
#' 
#' Retrieve a coefficient matrix \code{H}
#' derived from factorization by rank value
#' and generate heatmap of its elements.
#' 
#' @param object Object of class \code{scNMFSet}.
#' @param rank Rank value for which the cell map is to be displayed. 
#'   The object must contain the corresponding slot: one element of 
#'             \code{coeff(object)[[k]]} for which 
#'             \code{ranks(object)[[k]]==rank}.
#' @param main Title of plot.
#' @param ... Other arguments to be passed to \code{\link{heatmap}}, 
#' \code{\link{image}},
#'            and \code{\link{plot}}.
#' @return \code{NULL}
#' @examples
#' set.seed(1)
#' x <- simulate_data(nfeatures=10,nsamples=c(20,20,60))
#' rownames(x) <- seq_len(10)
#' colnames(x) <- seq_len(100)
#' s <- scNMFSet(count=x,rowData=seq_len(10),colData=seq_len(100))
#' s <- vb_factorize(s,ranks=seq(2,5))
#' plot(s)
#' cell_map(s, rank=3)
#' @export
cell_map <- function(object, rank, main = 'Cells', ...){
  
  if(missing(rank)) rank <- ranks(object)[1]
  h <- coeff(object)[ranks(object) == rank][[1]]
  nrow <- nrow(h)
  ccol <- grDevices::rainbow(n = nrow)
  cid <- apply(h,2, which.max)
  labCol <- if(dim(h)[2]>10) "" else colnames(h)
  stats::heatmap(h, Rowv = NA, ColSideColors = ccol[cid], RowSideColors = ccol, 
    revC = TRUE,labCol = labCol, main = main, 
    col = RColorBrewer::brewer.pal(n=9,'YlOrRd'), ...)
}

#' Find metagenes from basis matrix
#' 
#' Retrieve a basis matrix from an object and find metagenes.
#' 
#' @param object Object of class \code{scNMFSet}.
#' @param rank Rank value for which metagenes are to be found.
#' @param basis.matrix Instead of an object containing basis 
#'     matrices, the matrix itself can be provided.
#' @param max.per.cluster Maximum number of metagenes per cluster.
#' @param gene_names Names of genes to replace row names of basis matrix.
#' @param subtract.mean Standardize the matrix elements with means 
#'        within each row.
#' @param log Use geometric mean and division instead of arithmetic 
#'        mean and subtraction with \code{subtract.mean}.
#' @return List of vectors each containing metagene names of clusters.
#' @examples
#' set.seed(1)
#' x <- simulate_data(nfeatures=10,nsamples=c(20,20,60))
#' rownames(x) <- seq_len(10)
#' colnames(x) <- seq_len(100)
#' s <- scNMFSet(count=x,rowData=seq_len(10),colData=seq_len(100))
#' s <- vb_factorize(s,ranks=seq(2,5))
#' meta_genes(s, rank=4)
#' @export
meta_genes <- function(object, rank, basis.matrix=NULL, 
                       max.per.cluster=10,gene_names=NULL,
                       subtract.mean=TRUE,log=TRUE){
  
  if(is.null(basis.matrix)){
    idx <- ranks(object)==rank
    w <- basis(object)[idx][[1]]
    if(subtract.mean){
      if(log) w <- log10(w)
      w <- w - rowMeans(w)
      if(log) w <- 10^w
    }
  } else{
    w <- basis.matrix
    rank <- ncol(w)
  }
  if(!is.null(gene_names)) rownames(w) <- gene_names
  nmax <- min(max.per.cluster, nrow(w))
  select <- vector('list',rank)
  for(k in seq_len(rank)){
    idx <- order(w[,k],decreasing=TRUE)
    v <- w[idx,]
    itmp <- NULL
    for(i in seq_len(nrow(v))){
        x <- v[i,] 
        ix <- order(x,decreasing=TRUE)
        flag <- k==ix[1] 
        if(!flag) next
        xp <- log(mean(exp(x[ix[-1]])))
        itmp <- c(itmp,i)
    }
    tmp <- rownames(v)[itmp]
    if(length(tmp) > nmax) tmp <- tmp[seq_len(nmax)]
    select[[k]] <- tmp
  }
  select
}

gene_select <- function(w, markers = NULL, max.per.cluster = 10){
  
  rank <- dim(w)[2]
  select <- c()
  for(k in seq_len(rank)){
    nmax <- min(max.per.cluster, dim(w)[1])
    if(!is.null(markers)){
      select <- c(select, markers[markers %in% rownames(w)])
      w <- w[!rownames(w) %in% markers, ]
    }
    v <- w[order(w[,k],decreasing=TRUE),]
    flag <- apply(v,1,function(x){x[k]==max(x)})
    tmp <- rownames(v)[which(flag)]
    if(length(tmp) > nmax) tmp <- tmp[seq_len(nmax)]
    select <- c(select, tmp)
  }
  unlist(select)
}

#' Visualize clusters
#' 
#' Use tSNE to generate two-dimensional map of coefficient matrix.
#' 
#' It retrieves a coefficient matrix \code{H} 
#' from an object and use its elements
#' to assign each cell into clusters. 
#' t-Distributed Stochastic Neighbor Embedding (t-SNE; 
#' \url{https://lvdmaaten.github.io/tsne/}) is used to visualize the 
#' clustering
#' in 2D. Also plotted is the distribution of cell counts for all clusters.
#' 
#' @param object \code{scNMF} object.
#' @param rank Rank value to extract from \code{object}.
#' @param verbose Print tSNE messages.
#' @param cex Symbol size in tSNE plot
#' @param cex.names Font size of labels in count barplot.
#' @param ... Other parameters to send to \code{Rtsne}.
#' @return \code{NULL}
#' @examples
#' set.seed(1)
#' x <- simulate_data(nfeatures=10,nsamples=c(20,20,60,40,30))
#' rownames(x) <- seq_len(10)
#' colnames(x) <- seq_len(170)
#' s <- scNMFSet(count=x,rowData=seq_len(10),colData=seq_len(170))
#' s <- vb_factorize(s,ranks=seq(2,5))
#' visualize_clusters(s,rank=5)
#' @export
#' @importFrom graphics par
#' @importFrom graphics points
visualize_clusters <- function(object, rank, verbose=FALSE, cex=1, 
                               cex.names = 0.7, ...){
  
  if(missing(rank)){
    rid <- 1
    rank <- ranks(object)[1]
  }
  h <- coeff(object)[ranks(object) == rank][[1]]
  tsne <- Rtsne::Rtsne(t(h), verbose=verbose, ...)
  cid <- apply(h,2,which.max)       # cluster assignment
  col <- grDevices::rainbow(n=dim(h)[1])
  color <- col[cid]
  
  par(mfrow=c(1,2))
  plot(x=tsne$Y, pch=21,bg=color,main='Clusters', xlab='tSNE1',ylab='tSNE2',
       lwd=0.5, cex=cex)
  graphics::barplot(table(cid), names.arg=rownames(h),
          col=col,log='y', main='Cell counts', cex.names=cex.names,
          las=2)
  return(invisible(object))
}

#' Generate simulated data for factorization
#' 
#' Use one of two schemes to generate simulated data suitable for 
#' testing factorization.
#' 
#' In one scheme (\code{generate.factors = TRUE}), simulated factor 
#'   matrices 
#' \code{W} and \code{H} are used to build count data \code{X = WH}.
#' In the second scheme, factor matrices are not used and \code{X} is 
#'   sampled directly from \code{r} (rank requested) sets of 
#'   multinomial distributions.
#' 
#' @param nfeatures Number of features \code{m} (e.g., genes).
#' @param nsamples Vector of sample sizes in each cluster. 
#'    Rank \code{r} is equal to 
#'        the length of this vector. Sum of elements is the total 
#'        sample size \code{n}.
#' @param generate.factors Generate factor matrices \code{W} and \code{H},
#'        each with dimension \code{n x r} and \code{r x n}. If \code{FALSE},
#'        factor matrices are not used and count data are generated 
#'        directly from
#'        \code{r} multinomials for \code{m} genes.
#' @param nfactor Total RNA count of multinomials for each cluster with 
#'        \code{generate.factors = FALSE}. Small \code{nfactor} will 
#'        yield sparse count matrix.
#' @param alpha0 Variance parameter of Dirichlet distribution from which 
#'    multinomial probabilities are sampled with 
#'    \code{generate.factors = FALSE}.
#' @param shuffle Randomly permute rows and columns of count matrix.
#' @return If \code{generate.factors = TRUE}, list of components 
#'   \code{w} (basis matrix, \code{nfeatures x rank}),
#'         \code{h} (coefficient matrix, \code{rank x ncells}, where 
#'         \code{ncells}
#'         is equal to \code{n}, the sum of \code{nsamples}), and 
#'         \code{x}, a matrix of Poisson deviates with mean \code{W x H}.
#'         If \code{generate.factors = FALSE}, only the count matrix 
#'         \code{x} is in the list.
#' @examples
#' set.seed(1)
#' x <- simulate_data(nfeatures=10,nsamples=c(20,20,60,40,30))
#' s <- scNMFSet(x)
#' s
#' @export
simulate_data <- function(nfeatures, nsamples, generate.factors=FALSE, 
                          nfactor=10, alpha0=0.5, shuffle=TRUE){
  rank <- length(nsamples)
  if(generate.factors){
    if(nfeatures<rank) stop('No. of features too small for rank requested.')
    a <- rep(floor(nfeatures/rank),rank-1)
    a <- c(a,nfeatures-sum(a))
    H <- W <- NULL
    for(k in seq_len(rank)){
      tmp <- matrix(0, nrow=rank, ncol=nsamples[k])
      tmp[k,] <- 1
      H <- cbind(H, tmp)
      if(a[k]==0) next
      tmp2 <- matrix(0, nrow=a[k], ncol=rank)
      tmp2[, k] <- stats::rmultinom(n=1, size=a[k]*5, 
                    prob=rep(1,a[k]))
      W <- rbind(W, tmp2)
    }
    wh=W %*% H
    x <- apply(wh, c(1,2), function(x){stats::rpois(n=1, lambda=x)})
    if(shuffle){
      cel <- sample(seq_len(ncol(x)),size=ncol(x),replace=FALSE)
      x <- x[,cel]
      H <- H[,cel]
      gen <- sample(seq_len(nrow(x)),size=nrow(x),replace=FALSE)
      x <- x[gen,]
      W <- W[gen,]
    }
    return(list(w=W, h=H, x=matrix(x,nrow=nfeatures, ncol=sum(nsamples))))
  } else{
    x <- NULL
    q <- gtools::rdirichlet(n=rank,alpha=rep(alpha0,nfeatures))
    for(k in seq_len(rank)){
      x <- cbind(x, stats::rmultinom(n=nsamples[k], 
                             size=nfeatures*nfactor, prob=q[k,]))
    }
    if(shuffle) 
      x <- x[,sample(seq_len(ncol(x)),size=ncol(x),replace=FALSE)]
    return(x)
  }
}

#' Simulate factor matrices and data using priors
#' 
#' Under Bayesian formulation, use prior distributions of factor matrices
#' and generate simulated data
#' 
#' Basis \code{W} and coefficient matrices \code{H} are sampled from 
#' gamma distributions (priors) with shape (\code{aw,ah}) and mean
#' (\code{bw,bh}) parameters. Count data \code{X} are sampled from Poisson
#' distribution with mean values given by \code{WH}.
#' 
#' @param nrow Number of features (genes).
#' @param ncol Number of cells (samples).
#' @param rank Rank (ncol of W, nrow of H).
#' @param aw Shape parameter of basis prior.
#' @param bw Mean of basis prior. Scale parameter is equal to \code{aw/bw}.
#' @param ah Shape parameter of coefficient prior.
#' @param bh Mean of coefficient prior. Scale parameter is equal to 
#'    \code{ah/bh}.
#' @return List with elements \code{w}, \code{h}, and \code{x}, each 
#'    containing basis, coefficient, and count matrices.
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50,ncol=100,rank=5)
#' s <- scNMFSet(count=x$x)
#' s <- vb_factorize(s,ranks=seq(2,8),nrun=5)
#' plot(s)  
#' @export
simulate_whx <- function(nrow, ncol, rank, aw=0.1, bw=1, ah=0.1, bh=1){
  
  w <- stats::rgamma(n=nrow*rank, shape=aw, scale=bw/aw)
  w <- matrix(w, nrow=nrow, ncol=rank)
  
  h <- stats::rgamma(n=rank*ncol, shape=ah, scale=bh/ah)
  h <- matrix(h, nrow=rank, ncol=ncol)
  
  lambda <- w%*%h
  x <- apply(lambda, seq(1,2), 
             function(x){stats::rpois(n=1,lambda=x)})
  i <- rowSums(x) > 0
  j <- colSums(x) > 0
  x <- x[i,j]
  w <- w[i,]
  h <- h[,j]
  rownames(x) <- seq_len(nrow(x))
  colnames(x) <- seq_len(ncol(x))
  
  list(w=w, h=h, x=x)
}

#' Write 10x data files
#' 
#' Use an object and write count and annotation files in 10x format.
#'  
#' @param object Object of class \code{scNMFSet} containing count data
#' @param dir Directory where files are to be written. 
#' @param count File name for count matrix. 
#' @param genes File name for gene annotation.
#' @param barcodes File name for cell annotation.
#' @param quote Suppress quotation marks in output files.
#' @return \code{NULL}
#' @examples
#' set.seed(1)
#' x <- matrix(rpois(n=12,lambda=3),4,3)
#' rownames(x) <- seq_len(4)
#' colnames(x) <- seq_len(3)
#' s <- scNMFSet(count=x,rowData=seq_len(4),colData=seq_len(3))
#' write_10x(s,dir='.')
#' @export
write_10x <- function(object, dir, count='matrix.mtx',
                      genes='genes.tsv', barcodes='barcodes.tsv',
                      quote=FALSE){
  
  count <- paste0(dir, '/', count)
  genes <- paste0(dir, '/', genes)
  barcodes <- paste0(dir, '/', barcodes)
  x <- counts(object)
  x <- as(x,'sparseMatrix')
  Matrix::writeMM(as(x,'Matrix'),file=count)
  g <- rowData(object)
  rownames(g) <- rownames(rowData(object))
  utils::write.table(g,file=genes,col.names=FALSE,quote=quote,sep=' ',
              row.names=FALSE)
  utils::write.table(colData(object),file=barcodes,col.names=FALSE,
                     quote=FALSE,sep=' ',row.names=FALSE)
  return(invisible(object))
}

#' Assign cells into clusters
#' 
#' Use factorization results in an object to assign cells into clusters.
#' 
#' @param object Object of class \code{scNMFSet}
#' @param rank Rank value whose factor matrices are to be used for 
#'        assignment.
#' @return Vector of length equal to the number of cells containing 
#'      cluster ID numbers of each cell.
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50,ncol=100,rank=5)
#' s <- scNMFSet(count=x$x)
#' s <- vb_factorize(s,ranks=seq(2,8),nrun=5)
#' cid <- cluster_id(s, rank=5)
#' table(cid)
#' @export
cluster_id <- function(object, rank=2){
  
  h <- coeff(object)[ranks(object) == rank][[1]]
  cid <- apply(h, 2, which.max)
  names(cid) <- colnames(object)
  cid
}
