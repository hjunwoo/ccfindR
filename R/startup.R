# Startup message
.onAttach <- function(lib, pkg){
  msg <- c(paste0(
"__________________________________________

  ccfindR (Cancer Clone findeR), v",
utils::packageVersion('ccfindR'),"  \n\n"
,
"  Institute for Health Informatics
  University of Minnesota                 
__________________________________________
"
  ))
  if(!interactive())
    msg[1] <- paste("Package 'ccfindR' version", 
                    utils::packageVersion('ccfindR'))
  packageStartupMessage(msg)
  invisible()
}