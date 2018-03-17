#' Class \code{scNMFSet} for storing input data and results
#'
#' \code{S4} class derived from \link{\code{SingleCellExperiment}}
#' that can store single-cell count matrix, 
#' gene and cell annotation data frames, and factorization factors as well 
#' as quality measures for rank determination.
#' 
#' @slot assays Named list for count matrix \code{counts}. 
#' @slot rowData \code{DataFrame} for gene (feature) 
#'    names and annotations in columns.
#' @slot colData \code{DataFrame} for cell IDs and other 
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
#' 
#' Other slots inherited from \code{SingleCellExperiment} class are 
#' not explicitly used.
#' 
#' @examples
#' # toy matrix
#' ngenes <- 8; ncells <- 5
#' mat <- matrix(rpois(n=ngenes*ncells,lambda=3),ngenes,ncells)
#' 
#' abc <- letters[1:ngenes]
#' ABC <- LETTERS[1:ncells] 
#' genes <- DataFrame(gene_id=abc)
#' cells <- DataFrame(cell_id=ABC)
#' rownames(mat) <- rownames(genes) <- abc
#' colnames(mat) <- rownames(cells) <- ABC
#' 
#' # create scNMFSet object
#' s <- scNMFSet(count=mat,rowData=genes,colData=cells)
#' # alternative ways
#' s2 <- scNMFSet(count=mat)
#' s2 <- scNMFSet(assays=list(counts=mat))        
#'
#' # show dimensions
#' dim(s)
#' 
#' # show slots
#' rowData(s)
#' 
#' # modify slots
#' colData(s) <- DataFrame(cell_id=1:ncells,
#'               cell_type=c(rep('tissue1',2),
#'                           rep('tissue2',ncells-2)))
#' # slots be accessed directly
#' s@colData  
#'
#' @return Object of class \code{scNMFSet}
#' @export scNMFSet
# @import Matrix
#' @import methods
#' @import SingleCellExperiment
setClass('scNMFSet',
         slots=c(ranks='vector',
                 basis='list',
                 coeff='list',
                 measure='data.frame'),
         contains='SingleCellExperiment')

#' Create \code{scNMFSet} object
#' 
#' Object derived from \link{\code{SingleCellExperiment}}
#' 
#' @param count Count matrix
#' @param ... Other parameters \link{\code{SingleCellExperiment}}
#' @param remove.zeros Remove empty rows and columns
#' @return Object of class \code{scNMFSet}.
#' @export
scNMFSet <- function(count=NULL, ..., remove.zeros=TRUE){
           if(!is.null(count)) 
             Object <- SingleCellExperiment(assays=list(counts=count), ...)
           else
             Object <- SingleCellExperiment(...)
           if(min(counts(Object))<0) stop('Count data contains negative values.')
           if(remove.zeros) Object <- remove_zeros(Object)
           return(new('scNMFSet', Object))
}

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
             callNextMethod()
          })

#' Accessor for count matrix
#' 
#' @param Object Object containing count matrix
#' @examples
#' s <- scNMFSet(count = matrix(rpoise(n=12, 3,4)))
#' counts(s)
#' 
#' @export
#' @import BiocGenerics
setMethod('counts',signature='scNMFSet',
     function(object){
       assay(object, i='counts')
     }
)
# Assignment of count matrix
# 
# @param object Object containing count
# @param value Matrix-like object for replacement
# @examples
# 
setReplaceMethod('counts',c('scNMFSet','ANY'),
     function(object,value){
       assay(object, i='counts') <- value
       object
     }
)
#' Rank values in an Object
#' 
#' Retrieve or set the rank values in an object
#' 
#' Ranks for which factorization has been performed are stored 
#' in slot \code{ranks} of \code{scNMFSet} object.
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
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
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
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
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
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s,ranks=2:4)
#' measure(s)
#' @export
setGeneric('measure', function(object) standardGeneric('measure'))
#' Rank measure accessor
#' 
#' @param object Object containing measure
#' @return Data frame of measure
setMethod('measure', signature='scNMFSet', function(object) object@measure)

#' Subsetting scNMFSet object
# 
# @name [
#' @aliases [,scNMFSet-method
#' @docType methods
#' @rdname subset-methods
#' @param x Object to be subsetted
#' @param i row index
#' @param j column index
#' @return Subsetted object
setMethod('[', 'scNMFSet', function(x, i, j){
              if(missing(i)) i=1:nrow(x)
              if(missing(j)) j=1:ncol(x)
              x <- callNextMethod()
              if(length(basis(x)) > 0){
                w <- basis(x)
                  for(k in 1:length(w)) 
                    w[[k]] <- w[[k]][i,]
                basis(x) <- w
              }
              if(length(coeff(x)) > 0){
                h <- coeff(x)
                  for(k in 1:length(h))
                    h[[k]] <- h[[k]][,j]
                coeff(x) <- h
              }
              return(x)
          })

setValidity('scNMFSet', function(object){
    valid <- TRUE
    if(length(ranks(object))!=length(basis(object)) |
       length(ranks(object))!=length(coeff(object)) |
       length(basis(object))!=length(coeff(object))){
         valid <- FALSE
         msg <- 'rank, basis, or coeff data length do not match.'
    }
    if(valid) TRUE else msg
})

#' Gene annotation accessor
#' 
#' @param x Object containing data
#' @return DataFrame of row annotation
#' @export
setMethod('rowData','scNMFSet', 
          function(x){callNextMethod()})
#' Gene annotation assignment
#' 
#' @param x Object containing data
#' @param value DataFrame of row annotation to be substituted
#' @export
#' @import SummarizedExperiment
setMethod('rowData<-','scNMFSet', 
          function(x, value){callNextMethod()})

#' Cell annotation accessor
#' 
#' @param x Object containing cell annotation
#' @examples
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' colData(s) <- letters[1:3]
#' colData(s)
#' @export
setMethod('colData','scNMFSet',
          function(x){callNextMethod()})

#' Cell annotation assignment
#' 
#' @param x Object containing cell annotation
#' @param value DataFrame to be substituted
#' @examples
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' colData(s) <- letters[1:3]
#' colData(s)
#' @export
setReplaceMethod('colData','scNMFSet',
          function(x, value){x@colData <- value})

#' Modify ranks
#'  
#' Can be used to access and modify ranks 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Rank values (vector) to be substituted
#' @return Input object with updated ranks
#' @examples
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
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
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
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
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
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
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
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


#' @describeIn scNMFSet Plot measures of an object. 
#' For quality measures derived from maximum likelihood inference, 
#' dispersion and cophenetic will be plotted separately.
#' 
#' For measure derived from Bayesian inference, log evidence as a function of 
#' rank values will be plotted. 
#' 
#' @param x Object containing measure
#' @return \code{NULL}
#' @export
#' @importFrom graphics plot
setMethod('plot',signature="scNMFSet",definition = 
            function(x){
    
     bayes <- names(x@measure)[2]=='evidence'   
     dflag <- FALSE
     mx <- x@measure
     if(bayes) ylab <- c('log(Evidence)')
     else ylab <- c('Likelihood','Dispersion', 'Cophenetic')
     
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
       graphics::points(x=mx$rank,y=mx$likelihood,pch=21,bg='white',type='b')
     
       graphics::plot("",xlim=c(mx$rank[1], mx$rank[length(mx$rank)]),
          ylim=c(min(mx$dispersion),max(mx$dispersion)), xlab='Rank',
          ylab=ylab[2], bty='n')
       graphics::points(x=mx$rank,y=mx$dispersion,pch=21,bg='white',type='b')
     
       plot("",xlim=c(mx$rank[1], mx$rank[length(mx$rank)]),
          ylim=c(min(mx$cophenetic),max(mx$cophenetic)), xlab='Rank',
          ylab=ylab[3], bty='n')
       graphics::points(x=mx$rank, y=mx$cophenetic, pch=21,bg='white',type='b')
     }
            
     return(invisible())
})

# Remove rows or columns that are empty from an object
#' @export
remove_zeros <- function(object){
  
  count <- counts(object)
  genes <- rowData(object)
  cells <- colData(object)
  
  gene0 <- Matrix::rowSums(count)==0
  cell0 <- Matrix::colSums(count)==0
  ng0 <- sum(gene0)
  nc0 <- sum(cell0)
  
  if(ng0+nc0>0){
    object <- object[!gene0,!cell0]
    if(ng0>0)
      cat(ng0,'empty genes removed\n')
    if(nc0>0)
      cat(nc0,'empty cells removed\n')
  }
  object
}