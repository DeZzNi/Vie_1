#' Rente discrète avec vecteur qx
#'
#' @description Il faut seulement entrer les paramètres et tout ce qu'Étienne demande normalement est calculé.
#' Les fonctions ne sont pas idiotproof, il faut rentrer les paramètres correctement.
#' @param x âge de l'individu (0 par défaut)
#' @param n durée de contrat (1)
#' @param m années différée(0)
#' @param g paiement de rente (1)
#' @param delta taux d'intérêt (0)
#' @param qx vecteur des qx
#' @param s valeurs pour évaluer Fz (0)
#' @param kappa kappa comme d'habitude (0)
#' @export

ren_dis_qx <- function(x=0, n=1, m=0, g=1, delta=0, qx, s=0, kappa=0){

  px <- 1 - qx
  v <- exp(-delta)

  #sx[i] = Pr(x > i)
  sx <- cumprod(px)

  #Pr(Tx > t)
  stx <- function(x,t)
    sx[x+t] / sx[x]

  #Fonction de valeur d'annuité discrète
  ann <- function(t)
    (1 - v^t) / (1 - v)

  k <- m:(m+n-1)

  #Vecteur de probabilités ordonné
  fz <- sapply(k, function(y) stx(x,y) - stx(x,y+1))
  fz[n] <- fz[n] + stx(x, m+n)
  fz <- c(1-stx(x,m), fz)
  fz_cum <- cumsum(fz)

  #Valeurs possibles ordonnées
  z <- g * v^m * sapply(1:n, function(y) ann(y))
  z <- c(0, z)

  ez <- sum(fz * z)
  ez2 <- sum(fz * z^2)
  va <- ez2 - (ez)^2

  Fzz <- sapply(s, function(y) sum(fz * (z <= y)))

  j <- sapply(kappa, function(y) which(fz_cum >= y)[1])

  var <- z[j]

  etronq <- sapply(j, function(y) sum(fz[(y+1):(n+1)] * z[(y+1):(n+1)]))
  autre <- var * (fz_cum[j] - kappa)
  tvar <- (etronq + autre) / (1-kappa)

  #Si un des j = n+1, var=tvar
  if(any(j == n+1)){
    i <- which(j == (n+1))
    tvar[i] <- var[i]
  }

  print(data.frame(z=z, fz=fz, Fz=fz_cum))
  print(data.frame(EZ=ez, var=va))
  print(data.frame(s=s, Fz=Fzz))
  print(data.frame(kappa=kappa, VaR=var, TVaR=tvar))
}

#' @rdname ren_dis_qx
#' @export

