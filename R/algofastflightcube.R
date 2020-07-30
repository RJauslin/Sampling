#' Title
#'
#' @param X
#' @param pik
#'
#' @return
#' @export
#'
#' @examples
algofastflightcube2 <- function(X, pik) {

  N = length(pik)
  p = round(length(X)/length(pik))
  X <- array(X, c(N, p))
  A <- X/pik
  B <- A[1:(p + 1), ]


  psik <- pik[1:(p + 1)]
  ind <- seq(1, p + 1, 1)
  pp = p + 2
  B <- array(B, c(p + 1, p))

  ### HERE MODIFY B SUCH THAT IN FACT WE DO NOT USE p but loop until

  while (pp <= N) {
    tmp <- reduxB(B)
    tmp$B
    B[tmp$ind_row,tmp$ind_col]
    B_tmp <- tmp$B

    psik_tmp <- psik[tmp$ind_row]
    psik[tmp$ind_row]<- jump(B_tmp, psik_tmp)

    liste <- (psik > (1 - EPS) | psik < EPS)
    i <- 0
    while (i <= (p) & pp <= N) {
      i = i + 1
      if (liste[i] == TRUE) {
        pik[ind[i]] = psik[i]
        psik[i] = pik[pp]
        B[i, ] = A[pp, ]
        B = array(B, c(p + 1, p))
        ind[i] = pp
        pp = pp + 1
      }
    }
    # print(ind)
    # image(as(B,"sparseMatrix"))
  }

  if (length(pik[(pik > EPS & pik < (1 - EPS))]) == (p + 1)){
    psik <- jump(B, psik)
    pik[ind] = psik
    pik
  }

}





#' Title
#'
#' @param X
#' @param pik
#'
#' @return
#' @export
#'
#' @examples
algofastflightcube <- function(X, pik) {
  N = length(pik)
  p = round(length(X)/length(pik))
  X <- array(X, c(N, p))
  A <- X/pik
  B <- A[1:(p + 1), ]

  psik <- pik[1:(p + 1)]
  ind <- seq(1, p + 1, 1)
  pp = p + 2
  B <- array(B, c(p + 1, p))
  while (pp <= N) {
    psik <- jump(B, psik)
    liste <- (psik > (1 - EPS) | psik < EPS)
    i <- 0
    while (i <= (p) & pp <= N) {
      i = i + 1
      if (liste[i] == TRUE) {
        pik[ind[i]] = psik[i]
        psik[i] = pik[pp]
        B[i, ] = A[pp, ]
        B = array(B, c(p + 1, p))
        ind[i] = pp
        pp = pp + 1
      }
    }
  }
  if (length(pik[(pik > EPS & pik < (1 - EPS))]) == (p + 1)){
    psik <- jump(B, psik)
    pik[ind] = psik
    pik
  }

}
