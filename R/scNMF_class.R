#' Class for storing input data and results
#'
#' \code{S4} class that can store single-cell count matrix, 
#' gene and cell 
#' annotation data frames, and factorization factors as well 
#' as quality measures
#' for rank determination.
#' 
#' @slot count Matrix containing count data. 
#' @slot genes Data frame with row names containing gene (feature) 
#'    names and annotations in columns.
#' @slot cells Data frame with cell IDs in row names and other 
#'    annotations in columns (e.g., barcodes, cell types).
#' @slot ranks Vector for rank values for which factorization has been 
#'   performed.
#' @slot basis List (of length equal to that of \code{ranks}) of 
#'   basis matrices \strong{W} from factorization; 
#'   dimension \code{nrow} x \code{rank},
#'       where \code{nrow} is no. of rows in \code{count}.
#' @slot coeff List (of length equal to that of \code{ranks}) of 
#'  coefficient matrices \strong{H} from factorization; 
#'       dimension \code{rank} x \code{ncol},
#'       where \code{ncol} is no. of columns in \code{count}.
#' @slot measure Data frame of factorization quality measures for 
#'   each rank (\code{likelihood} and \code{dispersion}).
#' @slot dmeasure Data frame of quality measures with randomized
#'       matrix counts (H0).
#' 
#' @examples
#' # toy matrix
#' ngenes <- 8; ncells <- 5
#' mat <- matrix(rpois(n=ngenes*ncells,lambda=3),ngenes,ncells)
#' 
#' abc <- letters[1:ngenes]
#' ABC <- LETTERS[1:ncells] 
#' genes <- data.frame(gene_id=abc)
#' cells <- data.frame(cell_id=ABC)
#' rownames(mat) <- rownames(genes) <- abc
#' colnames(mat) <- rownames(cells) <- ABC
#' 
#' # create scNMFSet object
#' s <- scNMFSet(count=mat,genes=genes,cells=cells)
#' # alternative ways (genes and cells slots will be filled from count)
#' s <- scNMFSet(count=mat)
#' s <- scNMFSet(mat)        
#'
#' # show dimensions
#' s
#' 
#' # show slots
#' head(genes(s))
#' 
#' # modify slots
#' cells(s) <- data.frame(cell_id=1:ncells,
#'             cell_type=c(rep('tissue1',2),
#'             rep('tissue2',ncells-2)))
#' # slots can be accessed directly
#' head(cells(s))  
#'
#' @return Object of class \code{scNMFSet}
#' @export scNMFSet
#' @import Matrix
#' @import methods
scNMFSet <- setClass('scNMFSet', 
         slots=c(count='Matrix',
                 cells='data.frame',
                 genes='data.frame',
                 ranks='vector',
                 basis='list',
                 coeff='list',
                 measure='data.frame',
                 dmeasure='data.frame'
         ))
#if cells or genes are missing, fill them in from the count matrix
setMethod('initialize',signature=('scNMFSet'),
         definition=function(.Object, count, remove.zeros=TRUE, ...){
           if(min(count)<0) stop('Count data contains negative values.')
           if(class(count)=='data.frame' | class(count)=='matrix') 
             count <- as(as.matrix(count),'dgCMatrix') # coerce count
           .Object <- callNextMethod(.Object,count=count,...)
           if(is.null(colnames(.Object@count))) 
             colnames(.Object@count) <- 1:ncol(.Object)
           if(is.null(rownames(.Object@count)))
             rownames(.Object@count) <- 1:nrow(.Object)
           if(is(count,'Matrix')){
             if(sum(dim(.Object@cells))==0){ # cells empty
               cells(.Object) <- data.frame(
                 cell_id=colnames(.Object@count))
               rownames(cells(.Object)) <- colnames(.Object@count)
             }
             if(sum(dim(.Object@genes))==0){ # genes empty
               genes(.Object) <- data.frame(gene_id=rownames(.Object@count))
               rownames(genes(.Object)) <- rownames(.Object@count)
             }
           }
           if(remove.zeros) .Object <- remove_zeros(.Object)
           .Object
      })

#' Display object
#' 
#' Display the class and dimension of an object
#'  
#' Object name itself on command line or \code{(show(object))} will display 
#' class and dimensionality
#' 
#' @param object Object of class \code{scNMFSet}
#' @return \code{NULL}
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' show(s)
#' @export
setMethod('show', signature='scNMFSet',
           definition=function(object){
             cat('An object of class ', class(object), '\n', sep='')
             cat(' ', nrow(object@genes), ' genes by ',
                 nrow(object@cells), ' cells.\n', sep='')
             invisible(NULL)
           })

#' Dimension of an object
#'
#' Retrieve dimension of an object.
#'
#' \code{(dim(object))} will return dimesions of \code{object}
 
#' @param x Object of class \code{scNMFSet}.
#' @return Vector of dimensions.
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' dim(s)
#' @export
setMethod('dim', signature='scNMFSet',
          definition=function(x){
             dim(count(x))
          })

#' Count matrix of Object
#'
#' Retrieve or set the count matrix of an object
#'
#' \code{count(object)} will return the count matrix. 
#' \code{count(object) <- value}
#' can be used to assign to the count slot of \code{object}.
#' 
#' @param object Object of class \code{scNMFSet}.
#' @return Either \code{NULL} or count matrix of class 
#' \code{dgCMatrix}.
#' @return Count matrix in sparse format
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' head(count(s))
#' @export
setGeneric('count', function(object) standardGeneric('count'))
#' Count matrix accessor
#' @param object Object containing count matrix
#' @return Count matrix
setMethod('count', signature='scNMFSet', function(object) object@count)

#' Feature set of Object
#' 
#' Retrieve or set the gene annotation of object
#' 
#' \code{genes(object)} will return the gene annotation. 
#' \code{genes(object) <- value}
#' can be used to assign to the \code{genes} slot of \code{object}.
#' @param object Object of class \code{scNMFSet}.
#' @return Either \code{NULL} or data frame.
#' @examples
#' mat <- matrix(rpois(n=12,lambda=3),4,3)
#' rownames(mat) <- letters[1:4]
#' s <- scNMFSet(mat)
#' head(genes(s))
#' @export  
setGeneric('genes', function(object) standardGeneric('genes'))
#' Gene annotation accessor
#' @param object Object containing gene annotation
#' @return Gene annotation data frame
setMethod('genes', signature='scNMFSet', function(object) object@genes)

#' Cell annotation of Object
#' 
#' Retrieve or set the cell annotation of object
#' 
#' \code{cells(object)} will return the cell annotation. 
#' \code{cells(object) <- value}
#' can be used to assign to the \code{cells} slot of \code{object}. 
#' @param object Object of class \code{scNMFSet}.
#' @return Either \code{NULL} or data frame.'
#' @examples
#' mat <- matrix(rpois(n=12,lambda=3),4,3)
#' colnames(mat) <- LETTERS[1:3]
#' s <- scNMFSet(mat)
#' head(cells(s))
#' @export
setGeneric('cells', function(object) standardGeneric('cells'))
#' Cell annotation accessor
#' @param object Object containing cell annotation
#' @return Cell annotation data frame
setMethod('cells', signature='scNMFSet', function(object) object@cells)

#' Rank values in an Object
#' 
#' Retrieve or set the rank values in an object
#' 
#' Ranks for which factorization has been performed are stored 
#' in slot \code{ranks}
#' of object of class \code{scNMFSet}.
#' \code{ranks(object)} will return the rank vector. 
#' \code{ranks(object) <- value}
#' can be used to modify it.
#' @param object Object of class \code{scNMFSet}.
#' @return Either \code{NULL} or vector.
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s,ranks=2:4)
#' ranks(s)
#' @export
setGeneric('ranks', function(object) standardGeneric('ranks'))
#' Rank accessor
#' @param object Object containing rank values
#' @return Vector of rank values
setMethod('ranks', signature='scNMFSet', function(object) object@ranks)

#' Basis matrices in an Object
#' 
#' Retrieve or set the basis matrices \code{W} from factorization 
#' in an object
#' 
#' After factorization, basis matrices corresponding to each rank value are
#' stored as elements of a list, which is in slot \code{basis} of object of
#' class \code{scNMFSet}.
#' \code{basis(object)} will return the list of matrices. 
#' \code{basis(object) <- value}
#' can be used to modify it.
#' @param object Object of class \code{scNMFSet}
#' @return Either \code{NULL} or a list of same length as 
#' \code{ranks(object)}, whose
#' elements are basis matrices derived from factorization 
#' under each rank value.
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s,ranks=2:4)
#' basis(s)[[1]]
#' @export
setGeneric('basis', function(object) standardGeneric('basis'))

#' Basis matrix accessor
#' 
#' @param object Object containing basis matrix
#' @return List of basis matrices
setMethod('basis', signature='scNMFSet', function(object) object@basis)

#' Coefficient matrices in an Object
#' 
#' Retrieve or set the coefficient matrices from factorization in an 
#' object
#' 
#' After factorization, coefficient matrices \code{H} corresponding 
#' to each rank value are
#' stored as elements of a list, which is in slot \code{coeff} of object of
#' class \code{scNMFSet}.
#' \code{coeff(object)} will return the list of matrices. 
#' \code{coeff(object) <- value}
#' can be used to modify it.
#' @param object Object of class \code{scNMFSet}.
#' @return Either \code{NULL} or a list of same length as 
#' \code{ranks(object)}, whose
#' elements are coefficient matrices derived from factorization 
#' under each rank value.
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s,ranks=2:4)
#' coeff(s)[[1]]
#' @export
setGeneric('coeff', function(object) standardGeneric('coeff'))
#' Coefficient matrix accessor
#' 
#' @param object Object containing coefficient matrix
#' @return List of coefficient matrices
setMethod('coeff', signature='scNMFSet', function(object) object@coeff)

#' Factorization measures in an Object
#' 
#' Retrieve or set factorization measures in an object
#' 
#' Factorization under multiple rank values lead to measures stored in
#' a data frame inside a slot \code{measure}. In maximum likelihood using
#' \code{\link{factorize}}, this set of quality measures include 
#' dispersion and cophenetic
#' coeeficients for each rank. In Bayesian factorization using 
#' \code{\link{vb_factorize}},
#' log evidence for each rank is stored. \code{measure(object)} 
#' will return the data
#' frame. \code{measure(object) <- value}
#' can be used to modify it.
#' @param object Object of class \code{scNMFSet}.
#' @return Either \code{NULL} or a data frame containing measures.
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s,ranks=2:4)
#' measure(s)
#' @export
setGeneric('measure', function(object) standardGeneric('measure'))
#' Rank measure accessor
#' 
#' @param object Object containing measure
#' @return Data frame of measure
setMethod('measure', signature='scNMFSet', function(object) object@measure)

#' Quality measures under null hypotheses
#' 
#' Retrieve or set net quality measures in an object
#' 
#' Factorization using maximum likelihood using \code{\link{factorize}} 
#' can be followed 
#' by \code{\link{factorize_H0}}, which repeats inference for all rank 
#' values with
#' randomly permuted data. The differences in quality measures after 
#' subtracting
#' the null values are stored in \code{dmeasure} slot. 
#' \code{dmeasure(object)} 
#' will return the data frame. \code{dmeasure(object) <- value}
#' can be used to modify it. Not used when factorization is done with 
#' \code{\link{vb_factorize}}.
#' 
#' @param object Object of class \code{scNMFSet}.
#' @return Either \code{NULL} or a data frame containing difference measures.
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' s <- factorize(s,ranks=2:4)
#' s <- factorize_H0(s)
#' dmeasure(s)
#' @export
setGeneric('dmeasure', function(object) standardGeneric('dmeasure'))
#' Null rank measure accessor
#' 
#' @param object Object containing null measure
#' @return Data frame of null measure
setMethod('dmeasure', signature='scNMFSet', function(object) object@dmeasure)

#' Subsetting scNMFSet object
# 
# @name [
#' @aliases [,scNMFSet-method
#' @docType methods
#' @rdname subset-methods
#' @param x Object to be subsetted
#' @param i row index
#' @param j column index
setMethod('[', 'scNMFSet', function(x, i, j){
                    if(missing(i)) i=1:(dim(x)[1])
                    if(missing(j)) j=1:(dim(x)[2])
                    .count <- x@count[i,j]
                    if(is.null(dim(.count))){ # became vector
                    .count <- data.frame(.count)
                    rownames(.count) <- rownames(x@count)[i]
                    colnames(.count) <- colnames(x@count)[j]
                  }
                  .cells <- x@cells[j, ]
                  if(is.null(dim(.cells))){ # became vector
                    .cells <- data.frame(cell_id=.cells)
                    rownames(.cells) <- rownames(x@cells)[j]                    
                  }
                  .genes <- x@genes[i, ]
                  if(is.null(dim(.genes))){ # became vector
                    .genes <- data.frame(gene_id=.genes)
                    rownames(.genes) <- rownames(x@genes)[i]                    
                  }
                  if(length(x@basis)>0){
                    w <- x@basis
                    for(k in 1:length(w)) 
                     w[[k]] <- w[[k]][i,]
                  }
                  if(length(x@coeff)>0){
                    h <- x@coeff
                    for(k in 1:length(h))
                      h[[k]] <- h[[k]][,j]
                  }
               obj <- scNMFSet(count=.count, cells=.cells, genes=.genes)
               obj@measure <- x@measure
               if(length(x@basis)>0) obj@basis <- w
               if(length(x@coeff)>0) obj@coeff <- h
               obj@ranks <- x@ranks
                 
               obj
          })

setValidity('scNMFSet', function(object){
   msg <- NULL
   valid <- TRUE
   if(sum(dim(object@cells)==0 & dim(object@genes))==0)  # counts only
      return(TRUE)
   if(nrow(count(object)) != nrow(genes(object))){
     valid <- FALSE
     msg <- c(msg,'Number of genes and count data rows must match.')
   }
   if(ncol(count(object)) != nrow(cells(object))){
     valid <- FALSE
     msg <- c(msg, 'Number of cells and count data columns must match.')
   }
   if(!identical(rownames(count(object)), rownames(genes(object)))){
     valid <- FALSE
     msg <- c(msg, 'Count data rows and gene names must match.')
   }
   if(!identical(colnames(count(object)),rownames(cells(object)))){
    valid <- FALSE
    msg <- c(msg, 'Count data columns and cell names must match.')
   }
   if(valid) TRUE else msg
})

#' Modify count matrix
#' 
#' Can be used to access and modify count matrix 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value count matrix to be substituted
#' @return Input object with updated count matrix
#' @examples
#' s <- scNMFSet(matrix(0, 4,3))
#' count(s) <- matrix(rpois(n=12,lambda=3),4,3)
#' count(s)
#' @export
setGeneric('count<-', function(object,value) standardGeneric('count<-'))
#' Count matrix assignment
#' 
#' @param object Object containing count matrix
#' @param value New count matrix to be substituted
#' @return Object with count modified
setMethod('count<-', 'scNMFSet', 
       function(object, value){
          if(class(value)!='dgCMatrix')
            value <- as(as.matrix(value),'dgCMatrix') # coerce count
          object@count <- value
          if(validObject(object)) return(object)
})

#' Modify cell annotation
#'  
#' Can be used to access and modify cell annotation 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value cell annotation data frame to be substituted
#' @return Input object with updated cell annotation
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' cells(s) <- letters[1:3]
#' cells(s)
#' @return Input object with updated cell annotation
#' @export
setGeneric('cells<-', function(object,value) standardGeneric('cells<-'))
#' Cell annotation assignment
#' 
#' @param object Object containing cell annotation
#' @param value New cell annotation to be substituted
#' @return Object with modified cell annotation
#' @export
setMethod('cells<-','scNMFSet',
       function(object, value){
         if(class(value) != 'data.frame')
           value <- data.frame(value)
         object@cells<- value
         if(validObject(object)) return(object)
})

#' Modify gene annotation
#'  
#' Can be used to access and modify gene annotation 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Gene annotation data frame to be substituted
#' @return Input object with updated gene annotation
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' genes(s) <- letters[1:3]
#' genes(s)
#' @return Input object with updated gene annotation
#' @export
setGeneric('genes<-', function(object,value) standardGeneric('genes<-'))
#' Gene annotation assignment
#' 
#' @param object Object containing gene annotation
#' @param value New gene annotation to be substituted
#' @return Object with modified gene annotation
setMethod('genes<-','scNMFSet',
       function(object, value){
         if(class(value) != 'data.frame')
           value <- data.frame(value)
         object@genes <- value
         if(validObject(object)) return(object)
})

#' Modify ranks
#'  
#' Can be used to access and modify ranks 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Rank values (vector) to be substituted
#' @return Input object with updated ranks
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=2:3)
#' ranks(s) <- c('two','three')
#' ranks(s)
#' @return Input object with updated ranks
#' @export
setGeneric('ranks<-', function(object,value) standardGeneric('ranks<-'))
#' Rank values assignment
#' 
#' @param object Object containing ranks
#' @param value New rank vector to be substituted
#' @return Object with rank vector modified
setMethod('ranks<-','scNMFSet',
       function(object, value){
         object@ranks <- value
         if(validObject(object)) return(object)
})

#' Modify basis matrices
#'  
#' Can be used to access and modify basis matrices 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Basis matrix to be substituted
#' @return Input object with updated basis matrices
#' @examples
#' set.seed(1)
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=3)
#' basis(s)[[1]] <- apply(basis(s)[[1]],1:2,round,digits=3)
#' basis(s)
#' @return Input object with updated basis matrices
#' @export
setGeneric('basis<-', function(object,value) standardGeneric('basis<-'))
#' Basis matrix assignment
#' 
#' @param object Object containing basis matrix
#' @param value New basis matrices to be substituted
#' @return Object with basis matrix modified
setMethod('basis<-','scNMFSet',
       function(object, value){
         object@basis <- value
         if(validObject(object)) return(object)
})

#' Modify coefficient matrices
#'  
#' Can be used to access and modify coefficient matrices 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Coefficient matrix to be substituted
#' @return Input object with updated coefficient matrices
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=3)
#' coeff(s)[[1]] <- apply(coeff(s)[[1]],1:2,round,digits=2)
#' coeff(s)
#' @return Input object with updated coefficient matrices
#' @export
setGeneric('coeff<-', function(object,value) standardGeneric('coeff<-'))
#' Coefficient matrix assignment
#'
#' @param object Object containing coefficient matrices
#' @param value New coefficient matrices to be substituted
#' @return Object with coefficient matrix modified
setMethod('coeff<-','scNMFSet',
       function(object, value){
         object@coeff<- value
         if(validObject(object)) return(object)
})

#' Modify factorization measure
#'  
#' Can be used to access and modify factorization measure
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Measure to be substituted
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=3)
#' measure(s)[,-1] <- apply(measure(s)[,-1],1:2, round,digits=3)
#' measure(s)
#' @return Input object with updated measure
#' @export
setGeneric('measure<-', function(object,value) standardGeneric('measure<-'))
#' Rank measure assignment
#' 
#' @param object Object containing measure
#' @param value New measure (data frame) to be substituted
#' @return Object with modified measure 
setMethod('measure<-','scNMFSet',
       function(object, value){
         object@measure<- value
         if(validObject(object)) return(object)
})

#' Modify factorization null measures
#'  
#' Can be used to access and modify factorization null measures 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Null measure (data frame) to be substituted
#' @return Input object with updated \code{dmeasure}
#' @examples
#' s <- scNMFSet(matrix(rpois(n=12,lambda=3),4,3))
#' s <- factorize(s, ranks=2:4)
#' s <- factorize_H0(s)
#' dmeasure(s)[,-1] <- apply(dmeasure(s)[,-1],1:2, round,digits=3)
#' dmeasure(s)
#' @return Input object with updated dmeasure
#' @export
setGeneric('dmeasure<-', function(object,value) standardGeneric('dmeasure<-'))
#' Null rank measure assignment
#' 
#' @param object Object containing null measure
#' @param value New null measure (data frame) to be substituted
#' @return Object with null measure modified
setMethod('dmeasure<-','scNMFSet',
       function(object, value){
         object@dmeasure<- value
         if(validObject(object)) return(object)
})

#' @describeIn scNMFSet Plot measures of an object. For quality 
#' measures derived from
#' maximum likelihood inference, dispersion and cophenetic will 
#' be plotted separately.
#' 
#' For measure derived from Bayesian inference, log evidence as a function of 
#' rank values will be plotted. Parameter 
#' \code{measure = c('measure','dmeasure')} can be used
#' when slot \code{dmeasure} has been filled (maximum likelihood NMF)
#' to plot either H1 or H0 quality measures.
#' 
#' @param x Object containing measure
#' @param measure \code{c('measure','dmeasure')}
#' @return \code{NULL}
#' @export
#' @importFrom graphics plot
setMethod('plot',signature="scNMFSet",definition = 
            function(x, measure='measure'){
    
     bayes <- names(x@measure)[2]=='evidence'   
     if(measure=='dmeasure'){
       dflag <- TRUE
       mx <- x@dmeasure
       ylab <- c(expression(Delta*' likelihood'),
                 expression(Delta*' dispersion'),
                 expression(Delta*' cophenetic'))
     }
     else{
       dflag <- FALSE
       mx <- x@measure
       if(bayes) ylab <- c('log(Evidence)')
       else ylab <- c('Likelihood','Dispersion', 'Cophenetic')
     }
     if(sum(dim(mx))==0) stop('Quality measure empty.')
    
     if(bayes){
       par(mfrow=c(1,1))
       graphics::plot(x=mx$rank, y=mx$evidence, type='b',xlab='Rank', ylab=ylab,
            bty='n') 
     }
     if(!bayes){   
       par(mfrow=c(1,3))
       graphics::plot("",xlim=c(mx$rank[1], mx$rank[length(mx$rank)]),
          ylim=c(min(mx$likelihood),max(mx$likelihood)), xlab='Rank',
          ylab=ylab[1], bty='n')
       if(dflag)
         if(!is.na(mx$r_se[1])) error_bar(mx[,c('ranks','likelihood','r_se')])
       graphics::points(x=mx$rank,y=mx$likelihood,pch=21,bg='white',type='b')
     
       graphics::plot("",xlim=c(mx$rank[1], mx$rank[length(mx$rank)]),
          ylim=c(min(mx$dispersion),max(mx$dispersion)), xlab='Rank',
          ylab=ylab[2], bty='n')
       if(dflag)
         if(!is.na(mx$d_se[1])) error_bar(mx[,c('ranks','dispersion','d_se')])
       graphics::points(x=mx$rank,y=mx$dispersion,pch=21,bg='white',type='b')
     
       plot("",xlim=c(mx$rank[1], mx$rank[length(mx$rank)]),
          ylim=c(min(mx$cophenetic),max(mx$cophenetic)), xlab='Rank',
          ylab=ylab[3], bty='n')
       if(dflag)
         if(!is.na(mx$c_se[1])) error_bar(mx[,c('ranks','cophenetic','c_se')])
       graphics::points(x=mx$rank, y=mx$cophenetic, pch=21,bg='white',type='b')
     }
            
     return(invisible())
})

# Function to set error bars in plot.
error_bar <- function(dat){
  
  graphics::segments(x0=dat[,1], x1=dat[,1], y0=dat[,2]-dat[,3], 
           y1=dat[,2]+dat[,3], col='red')
  graphics::segments(x0=dat[,1]-0.2, x1=dat[,1]+0.2, y0=dat[,2]-dat[,3], 
           y1=dat[,2]-dat[,3], col='red')
  graphics::segments(x0=dat[,1]-0.2, x1=dat[,1]+0.2, y0=dat[,2]+dat[,3], 
           y1=dat[,2]+dat[,3], col='red')

}

# Remove rows or columns that are empty from an object
remove_zeros <- function(object){
  
  count <- object@count
  genes <- object@genes
  cells <- object@cells
  
  gene0 <- Matrix::rowSums(count)==0
  cell0 <- Matrix::colSums(count)==0
  ng0 <- sum(gene0)
  nc0 <- sum(cell0)
  
  if(ng0+nc0>0){
    object@count <- count[!gene0,!cell0]
    if(ng0>0){
      tmp <- as.data.frame(genes[!gene0,])
      names(tmp) <- names(genes)
      object@genes <- tmp
      cat(ng0,'empty genes removed\n')
    }
    if(nc0>0){
      tmp <- as.data.frame(cells[!cell0,])
      names(tmp) <- names(cells)
      object@cells <- tmp
      cat(nc0,'empty cells removed\n')
    }
  }
  object
}