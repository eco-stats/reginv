#'@export

print.reginv <- function(x, ...) {
  cat("\n Call: ")
  print(x$call)
  cat("\n")
  print(x$theta)
}
