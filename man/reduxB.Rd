% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/redux.R
\name{reduxB}
\alias{reduxB}
\title{Reduction of the matrix}
\usage{
reduxB(B)
}
\arguments{
\item{B}{a matrix of size (p+1 x p) sub-matrix of auxiliary matrix.}
}
\value{
a list
}
\description{
This function reduces the size of the matrix by removing alternatively columns and rows that have sum equal to 0.

In case where the number of auxiliary varibale is great (p very large), even if we use the fast implementation proposed by
(Chauvet and Tillé 2005) the problem is time consuming. If we have the chance that the matrix is strongly sparse,
we can then use the function to reduce the size of the matrix B by using this method. 

If the matrix is dense or every column have sum greater than 0, then nothing is changed.
}
\examples{
\dontrun{
#' rm(list = ls())
set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 50
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),5),
sampling::inclusionprobabilities(runif(N),10),
sampling::inclusionprobabilities(runif(N),15)),ncol = 3)
X <- PM(Pik)$PM
pik <- PM(Pik)$P
dim(X)
order = 2
EPS = 1e-11


p <- ncol(X)
A <-  X/pik
B <- A[1:(p + 1), ]

tmp <- reduxArma(B)
tmp2 <- reduxB(B)
B_redux <- tmp$B
B_redux[1:10,1:20]
B[tmp$ind_row[1:10],tmp$ind_col[1:20]]
}
}
\references{
Chauvet, G. and Tillé, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53–62.
}
