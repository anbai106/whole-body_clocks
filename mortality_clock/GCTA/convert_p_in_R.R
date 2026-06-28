#!/usr/bin/env Rscript
compute_p <- function(lrt){
  p <- 0.5 * pchisq(lrt, df=1, lower.tail=FALSE)
  return(p)
}
