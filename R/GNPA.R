#' Générateur de nombres pseudo-aléatoire
#'
#' @description Les fonctions ne sont pas idiotproof, il faut rentrer les paramètres correctement.
#' @param n nombre d'uniformes
#' @param a multiplicateur (41358)
#' @param m modulo (2^31 - 1)
#' @param x0 valeur de départ
#' @export

GNPA <- function(n, a=41358, m=2147483647, x0){

  x <- numeric(n+1)
  x[1] <- x0

  for (i in 1:n)
    x[i+1] <- (a * x[i]) %% m

  x[-1] / m
}
#' @rdname GNPA
#' @export
