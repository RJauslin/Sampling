#' Internal function of algofastflightcube
#' @noRd
jump <- function(X, pik) {
  N = length(pik)
  p = ncol(X)
  
  # add a col of 0 to be sure that it is not empty
  X1 = cbind(X, rep(0, times = N))
  u = MASS::Null(X1)[,1]

  l1=min(pmax((1-pik)/u,-pik/u))
  l2=min(pmax((pik-1)/u,pik/u))
  if(runif(1)<l2/(l1+l2)){
    pik = pik+l1*u
  }else{
    pik = pik-l2*u
  }

  return(pik)

}
