#' Internal function of algofastflightcube
#' @noRd
jump <- function(X, pik) {
  EPS = 1e-11
  N = length(pik)
  p = round(length(X)/length(pik))
  X <- array(X, c(N, p))
  X1 = cbind(X, rep(0, times = N))

  kern <- svd(X1)$u[, p + 1]

  listek = abs(kern) > EPS

  buff1 <- (1 - pik[listek])/kern[listek]
  buff2 <- -pik[listek]/kern[listek]
  la1 <- min(c(buff1[(buff1 > 0)], buff2[(buff2 > 0)]))
  pik1 <- pik + la1 * kern

  buff1 <- -(1 - pik[listek])/kern[listek]
  buff2 <- pik[listek]/kern[listek]
  la2 <- min(c(buff1[(buff1 > 0)], buff2[(buff2 > 0)]))
  pik2 <- pik - la2 * kern

  q <- la2/(la1 + la2)
  if (stats::runif(1) < q)
    pikn <- pik1
  else pikn <- pik2
  pikn
}
