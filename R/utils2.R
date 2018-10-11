#' Determine optimal rank
#' 
#' Takes as main argument \code{scNMFSet} object containing factorization output
#' and attempts to determine optimal rank.
#' 
#' @param object \code{scNMFSet} object containing factorization output.
#' @param df Degree of freedom for spline fit. Upper bound is the total 
#'        number of data (number of rank values scanned). 
#' @param type \code{c(1,2)}. Type 1 is where there is a clear maximum. 
#'        Type 2 is where marginal likelihood reaches a maximal level and 
#'        stays constant.
#' @param max.slope For Type 2, optimal rank is defined as the minimum rank 
#'        where the slope of log ML with respect to rank becomes smaller 
#'        than this parameter.
#' @return Optimal rank.
#' 
#' @details Spline fit is used for log marginal likelihood versus rank data. 
#'        For Type 1, the location of maximum is returned.
#'        For Type 2, lowest rank with slope smaller than specified value will be returned.
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50,ncol=100,rank=5)
#' s <- scNMFSet(x$x)
#' s <- vb_factorize(s,ranks=seq(2,8),nrun=5)
#' plot(s)
#' optimal_rank(s)
#' @export
optimal_rank <- function(object, df=10, type=1, max.slope=1e-4){
  
  if(is(object,'scNMFSet'))
    me <- measure(object)[,seq_len(2)]
  else if(is(object,'data.frame'))
    me <- object[,seq_len(2)]
  else stop('Inappropriate class of object')
  
  df <- min(df,nrow(me))
  fs <- stats::smooth.spline(x=me[,1],y=me[,2],df=df)
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
