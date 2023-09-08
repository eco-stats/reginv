#'@export
print.reginv <- function(x, ...) {
  cat("\n Call: ")
  print(x$call)
  cat("\n")
  print(x$theta)
}

#'@export
print.mle_fossil <- function(x, ...) {
  cat("\n Call: ")
  print(x$call)
  cat("\n  MLE:", x$theta[1],"\n")
  if(is.null(x$call$alpha)==FALSE)
  {
    cat("\n   CI:\n")
    print(x$ci)
  }
}
