#' S3 plotting function for `CVIoutput`objects'
#'
#' @param x a CVIoutput object
#'
#' @export

plot.CVIoutput <- function(x, y=NULL, ...){

  plot(sapply(x$optimisation$ELBO, sum)[-1], type="l",
       xlab="Iterations", ylab="ELBO (log)")

}
