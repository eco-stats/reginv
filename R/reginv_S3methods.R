#'@export
print.reginv <- function(x, ...) {
  cat("\n Call: ")
  print(x$call)
  cat("\n")
  print(x$theta)
  cat(" CI estimated using", x$method, "method\n")
}

#'@export
print.est_cutt <- function(x, ...) {
  cat("\n Call: ")
  print(x$call)
  if(is.null(x$q)==FALSE)
  {
    cat("\n   CI:\n")
    print(x$theta)
  }
  else
    cat("\n  MLE:", x$theta[1],"\n")
  librr
  cat(" CI estimated using", x$method, "method\n")
}

#'@export
qqenvelope.est_cutt <- function(x, ...) {
#insert function here
# also make qqenvelope a generic I guess...
}
