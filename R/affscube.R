#' Internal function of fastflightcubeSPOT
#' @noRd
affscube <- function(X, pik,redux = TRUE) {
  
  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------
  
  EPS = 1e-11
  N = length(pik)
  p = ncol(X)
  X <- array(X, c(N, p))
  A <- X/pik
  B <- A[1:(p + 1), ]
  
  
  psik <- pik[1:(p + 1)]
  ind <- seq(1, p + 1, 1)
  pp = p + 2
  B <- array(B, c(p + 1, p))
  
  
  ##---------------------------------------------------------------
  ##                          Main loop                           -
  ##---------------------------------------------------------------
  
  while (pp <= N) {
    
    #---- if redux is needed
    if(redux == TRUE){
      tmp <- reduxB(B)
      B_tmp <- tmp$B
      psik_tmp <- psik[tmp$ind_row]
      psik[tmp$ind_row]<- jump(B_tmp, psik_tmp)
    }else{
      psik <- jump(B, psik)
    }
    
    #---- update
    
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
  
  
  ##---------------------------------------------------------------
  ##                          Last step                           -
  ##---------------------------------------------------------------
  
  if(length(pik[(pik > EPS & pik < (1 - EPS))]) == (p + 1)){
    psik <- jump(B, psik)
    pik[ind] = psik
    pik
  }
  return(pik)
}
