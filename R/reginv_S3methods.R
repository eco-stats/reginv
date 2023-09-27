#'@export
print.reginv <- function(x, ...) {
  cat("\n Call: ")
  print(x$call)
  cat("\n")
  print(x$theta)
  cat(" CI estimated using", x$method, "method\n")
}

#'@export
print.mle_cutt <- function(x, ...) {
  cat("\n Call: ")
  print(x$call)
  cat("\n  MLE:", x$theta[1],"\n")
  if(is.null(x$q)==FALSE)
  {
    cat("\n   CI:\n")
    print(x$ci)
  }
}
