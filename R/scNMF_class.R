#' Class \code{scNMFSet} for storing input data and results
#'
#' \code{S4} class derived from \code{\link{SingleCellExperiment}}
#' that can store single-cell count matrix, gene and cell annotation 
#' data frames, and factorization factors as well as quality measures 
#' for rank determination.
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
#' library(S4Vectors)
#' # toy matrix
#' ngenes <- 8 
#' ncells <- 5
#' mat <- matrix(rpois(n=ngenes*ncells,lambda=3),ngenes,ncells)
#' 
#' abc <- letters[seq_len(ngenes)]
#' ABC <- LETTERS[seq_len(ncells)] 
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
#' colData(s) <- DataFrame(cell_id=seq_len(ncells),
#'               cell_type=c(rep('tissue1',2),
#'                           rep('tissue2',ncells-2)))
#' colData(s)
#' @return Object of class \code{scNMFSet}
#' @export scNMFSet
#' @import SingleCellExperiment
#' @import S4Vectors
#' @import methods
setClass('scNMFSet',
         slots=c(ranks='vector',
                 basis='list', dbasis='list',
                 coeff='list', dcoeff='list',
                 measure='data.frame'),
         contains='SingleCellExperiment')

#' Create \code{scNMFSet} object
#'
#' Object derived from \code{\link{SingleCellExperiment}}
#' 
#' @param count Count matrix
#' @param ... Other parameters of \code{\link{SingleCellExperiment}}
#' @param remove.zeros Remove empty rows and columns
#' @return Object of class \code{scNMFSet}.
#' @examples
#' count <- matrix(rpois(n=12,lambda=2),4,3)
#' s <- scNMFSet(count=count)
#' s
#' @export
scNMFSet <- function(count=NULL, ..., remove.zeros=TRUE){
           if(!is.null(count))
             Object <- SingleCellExperiment(assays=list(counts=count), ...)
           else
             Object <- SingleCellExperiment(...)
           if(min(counts(Object))<0) 
             stop('Count data contains negative values.')
           if(remove.zeros) Object <- remove_zeros(Object)
           x <- new('scNMFSet', Object)
           return(x)
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
             cat('class:', class(object),'\n')
             cat('dim:', dim(object),'\n')
             cat('rownames: ')
             print(utils::head(rownames(object)))
             cat('colnames: ')
             print(utils::head(colnames(object)))
          })

#' Accessor for count matrix
#' 
#' @param object Object containing count matrix
#' @return Count matrix
#' @examples
#' s <- scNMFSet(count = matrix(rpois(n=12,lambda=3),3,4))
#' counts(s)
#' 
#' @export
setMethod('counts',signature='scNMFSet',
     function(object){
       assay(object, i='counts')
     }
)
#' Assignment of count matrix
#' 
#' Count matrix can be modified
#' 
#' @param object Object containing count
#' @param value Matrix-like object for replacement
#' @return Object with updated count
#' @export
#' @examples
#' mat <- matrix(rpois(n=12,lambda=3),3,4)
#' s <- scNMFSet(count = mat)
#' counts(s) <- mat^2
#' counts(s)
setMethod('counts<-','scNMFSet',
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
#' s <- vb_factorize(s,ranks=seq(2,4))
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
#' s <- vb_factorize(s,ranks=seq(2,4))
#' basis(s)[[1]]
#' @export
setGeneric('basis', function(object) standardGeneric('basis'))

#' Basis matrix accessor
#' 
#' @param object Object containing basis matrix
#' @return List of basis matrices
setMethod('basis', signature='scNMFSet', function(object) object@basis)

#' Basis SD matrix accessor
#' 
#' @param object Object containing dbasis matrix
#' @return List of dbasis matrices
#' @export
setGeneric('dbasis', function(object) standardGeneric('dbasis'))
#' Basis SD matrix accessor
#' 
#' @param object Object containing basis standard deviation (SD) matrix
#' @return List of dbasis matrices
setMethod('dbasis', signature='scNMFSet', function(object) object@dbasis)

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
#' s <- vb_factorize(s,ranks=seq(2,4))
#' coeff(s)[[1]]
#' @export
setGeneric('coeff', function(object) standardGeneric('coeff'))
#' Coefficient matrix accessor
#' 
#' @param object Object containing coefficient matrix
#' @return List of coefficient matrices
setMethod('coeff', signature='scNMFSet', function(object) object@coeff)

#' Coeff SD matrix accessor
#' 
#' @param object Object containing dcoeff matrix
#' @return List of dcoeff matrices
#' @export
setGeneric('dcoeff', function(object) standardGeneric('dcoeff'))

#' Coeffcient SD matrix accessor
#' 
#' @param object Object containing coeffient standard deviation (SD) matrix
#' @return List of dcoeff matrices
setMethod('dcoeff', signature='scNMFSet', function(object) object@dcoeff)

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
#' s <- vb_factorize(s,ranks=seq(2,4))
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
              if(missing(i)) i=seq_len(nrow(x))
              if(missing(j)) j=seq_len(ncol(x))
              x <- callNextMethod()
              if(length(basis(x)) > 0){
                w <- basis(x)
                dw <- dbasis(x)
                for(k in seq_len(length(w))){ 
                  w[[k]] <- w[[k]][i,]
                  dw[[k]] <- dw[[k]][i,]
                }
                basis(x) <- w
                dbasis(x) <- dw
              }
              if(length(coeff(x)) > 0){
                h <- coeff(x)
                dh <- dcoeff(x)
                for(k in seq_len(length(h))){
                  h[[k]] <- h[[k]][,j]
                  dh[[k]] <- dh[[k]][,j]
                }
                coeff(x) <- h
                dcoeff(x) <- dh
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

#' Feature annotation accessor
#' 
#' @param x Object containing data
#' @return DataFrame of feature annotation
#' @examples
#' x <- matrix(rpois(n=12,lambda=3),4,3)
#' rownames(x) <- seq_len(4)
#' colnames(x) <- seq_len(3)
#' s <- scNMFSet(count=x,rowData=seq_len(4),colData=seq_len(3))
#' rowData(s) 
#' @export
setMethod('rowData','scNMFSet', 
          function(x){
            callNextMethod()})
#' Gene annotation assignment
#' 
#' @param x Object containing data
#' @param value DataFrame of row annotation to be substituted
#' @return Row annotation DataFrame
#' @export
#' @import SummarizedExperiment
setMethod('rowData<-','scNMFSet', 
          function(x, value){
            callNextMethod()})

#' Sample annotation accessor
#' 
#' @param x Object containing sample annotation
#' @return Column annotation DataFrame 
#' @examples
#' library(S4Vectors)
#' x <- matrix(rpois(n=12,lambda=3),4,3)
#' rownames(x) <- seq_len(4)
#' colnames(x) <- c('a','b','c')
#' s <- scNMFSet(count=x,rowData=seq_len(4),colData=c('a','b','c'))
#' cols <- DataFrame(tissue=c('tissue1','tissue1','tissue2'))
#' rownames(cols) <- c('a','b','c')
#' colData(s) <- cols
#' s
#' @export
setMethod('colData','scNMFSet',
          function(x){callNextMethod()})

#' Cell annotation assignment
#' 
#' @param x Object containing cell annotation
#' @param value DataFrame to be substituted
#' @return Updated column annotation 
#' @examples
#' library(S4Vectors)
#' x <- matrix(rpois(n=12,lambda=3),4,3)
#' rownames(x) <- seq_len(4)
#' colnames(x) <- c('a','b','c')
#' s <- scNMFSet(count=x,rowData=seq_len(4),colData=c('a','b','c'))
#' cols <- DataFrame(tissue=c('tissue1','tissue1','tissue2'))
#' rownames(cols) <- c('a','b','c')
#' colData(s) <- cols
#' s
#' @export
setMethod('colData<-','scNMFSet', 
          function(x, value){
            callNextMethod()})

#' Generics for ranks assignment
#'
#' Replace \code{ranks} slot of \code{scNMFSet} object
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Rank values (vector) to be substituted
#' @return Input object with updated ranks
#' @examples
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=seq(2,3))
#' ranks(s) <- c('two','three')
#' ranks(s)
#' @export
setGeneric('ranks<-', function(object,value) standardGeneric('ranks<-'))
#' Modify ranks
#'
#' Replace \code{ranks} slot of \code{scNMFSet} object
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Rank values (vector) to be substituted
#' @return Input object with updated ranks
#' @examples
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=seq(2,3))
#' ranks(s) <- c('two','three')
#' ranks(s)
#' @export
setMethod('ranks<-','scNMFSet',
       function(object, value){
         object@ranks <- value
         if(validObject(object)) return(object)
})

#' Generics for basis matrix assignment
#'
#' Access and modify basis matrices 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Basis matrix to be substituted
#' @return Input object with updated basis matrices
#' @examples
#' set.seed(1)
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=3)
#' basis(s)[[1]] <- apply(basis(s)[[1]],seq(1,2),round,digits=3)
#' basis(s)
#' @export
setGeneric('basis<-', function(object,value) standardGeneric('basis<-'))
#' Modify basis matrices
#'
#' Access and modify basis matrices 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Basis matrix to be substituted
#' @return Input object with updated basis matrices
#' @examples
#' set.seed(1)
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=3)
#' basis(s)[[1]] <- apply(basis(s)[[1]],c(1,2),round,digits=3)
#' basis(s)
#' @export
setMethod('basis<-','scNMFSet',
       function(object, value){
         object@basis <- value
#        if(validObject(object)) return(object)
         return(object)
         })
#' Basis SD matrix assignment
#' 
#' @param object Object containing dbasis matrix
#' @param value List for assignment
#' @return Updated object
#' @export
setGeneric('dbasis<-', function(object,value) standardGeneric('dbasis<-'))
#' Modify dbasis matrices
#'
#' Access and modify dbasis matrices 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Basis SD matrix to be substituted
#' @return Modified object
setMethod('dbasis<-','scNMFSet',
          function(object, value){
            object@dbasis <- value
            return(object)
          })

#' Generics for coefficient matrix assignment
#'
#' Access and modify coefficient matrices 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Coefficient matrix to be substituted
#' @return Input object with updated coefficient matrices
#' @examples
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=3)
#' coeff(s)[[1]] <- apply(coeff(s)[[1]],c(1,2),round,digits=2)
#' coeff(s)
#' @export
setGeneric('coeff<-', function(object,value) standardGeneric('coeff<-'))
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
#' coeff(s)[[1]] <- apply(coeff(s)[[1]],c(1,2),round,digits=2)
#' coeff(s)
#' @export
setMethod('coeff<-','scNMFSet',
       function(object, value){
         object@coeff<- value
#        if(validObject(object)) return(object)
         return(object)
})
#' Coeff SD matrix assignment
#' 
#' @param object Object containing dcoeff matrix
#' @param value List for assignment
#' @return Updated object
setGeneric('dcoeff<-', function(object,value) standardGeneric('dcoeff<-'))
#' Modify dcoeff matrices
#'
#' Access and modify dcoeff matrices 
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Coeff SD matrix to be substituted
#' @return Updated object
setMethod('dcoeff<-','scNMFSet',
          function(object, value){
            object@dcoeff <- value
            return(object)
          })

#' Generics for factorization measure assignment
#'
#' Can be used to access and modify factorization measure
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Measure to be substituted
#' @examples
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=3)
#' measure(s)[,-1] <- apply(measure(s)[,-1], c(1,2), round,digits=3)
#' measure(s)
#' @return Input object with updated measure
#' @export
setGeneric('measure<-', function(object,value) standardGeneric('measure<-'))
#' Modify factorization measure
#'
#' Can be used to access and modify factorization measure
#' 
#' @param object Object of class \code{scNMFSet}
#' @param value Measure to be substituted
#' @examples
#' s <- scNMFSet(count=matrix(rpois(n=12,lambda=3),4,3))
#' s <- vb_factorize(s, ranks=3)
#' measure(s)[,-1] <- apply(measure(s)[,-1], c(1,2), round,digits=3)
#' measure(s)
#' @return Input object with updated measure
#' @export
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
    
     bayes <- names(x@measure)[2] %in% c('lml','E','evidence')   
     dflag <- FALSE
     mx <- x@measure
     if(bayes) ylab <- c('log ML')
     else ylab <- c('Likelihood','Dispersion', 'Cophenetic')
     
     if(sum(dim(mx))==0) stop('Quality measure empty.')
    
     if(bayes){
       graphics::plot(NULL,xlim=c(mx$rank[1],mx$rank[nrow(mx)]),
                      ylim=c(min(mx[,2]),max(mx[,2])),xlab='Rank', 
                      ylab=ylab, bty='n')
       graphics::points(x=mx$rank, y=mx[,2], type='b',xlab='Rank', 
                      ylab=ylab, bty='n',pch=21,bg='white')
     }
     else{   
       par(mfrow=c(1,3))
       graphics::plot("",xlim=c(mx$rank[1], mx$rank[length(mx$rank)]),
          ylim=c(min(mx$likelihood),max(mx$likelihood)), xlab='Rank',
          ylab=ylab[1], bty='n')
       graphics::points(x=mx$rank,y=mx$likelihood,pch=21,
                        bg='white',type='b')
     
       graphics::plot("",xlim=c(mx$rank[1], mx$rank[length(mx$rank)]),
          ylim=c(min(mx$dispersion),max(mx$dispersion)), xlab='Rank',
          ylab=ylab[2], bty='n')
       graphics::points(x=mx$rank,y=mx$dispersion,pch=21,bg='white',
                        type='b')
     
       plot("",xlim=c(mx$rank[1], mx$rank[length(mx$rank)]),
          ylim=c(min(mx$cophenetic),max(mx$cophenetic)), xlab='Rank',
          ylab=ylab[3], bty='n')
       graphics::points(x=mx$rank, y=mx$cophenetic, pch=21,bg='white',
                        type='b')
     }
            
     return(invisible())
})

#' Remove rows or columns that are empty from an object
#'
#' @param object Object containing data
#' @return Object with empty rows/columns removed
#' @examples
#' set.seed(1)
#' x <- matrix(rpois(n=100,lambda=0.1),10,10)
#' s <- scNMFSet(count=x,remove.zeros=FALSE)
#' s2 <- remove_zeros(s)
#' s2
#' @export
remove_zeros <- function(object){

  if(is(object,'matrix'))
    count <- object
  else
    count <- counts(object)

  gene0 <- Matrix::rowSums(count)==0
  cell0 <- Matrix::colSums(count)==0
  ng0 <- sum(gene0)
  nc0 <- sum(cell0)
  
  if(ng0+nc0>0){
    object <- object[!gene0,!cell0]
    if(ng0>0)
      cat(ng0,'empty rows removed\n')
    if(nc0>0)
        cat(nc0,'empty cells removed\n')
  }
  object
}
