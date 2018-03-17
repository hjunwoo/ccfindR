## ------------------------------------------------------------------------
library(ccfindR)

## ---- eval=T,echo=T------------------------------------------------------
# A toy matrix for count data
set.seed(1)
mat <- matrix(rpois(n = 80, lambda = 2), nrow = 4, ncol = 20)
ABC <- LETTERS[1:4]
abc <- letters[1:20]

