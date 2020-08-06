
rm(list = ls())
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




rm(list = ls())
set.seed(1)
B <- as.matrix(rsparsematrix(200,200,density = 0.01))
image(as(B,"sparseMatrix"))
system.time(test1 <- reduxArma(B))
dim(test1$B)
system.time(test2 <- reduxB(B))
dim(test2$B)



rm(list = ls())
set.seed(1)
B <- as.matrix(rsparsematrix(200,200,density = 0.001))
B <- cbind(rep(1,200),B)
image(as(B,"sparseMatrix"))
system.time(test1 <- reduxArma(B))
test1
dim(test1$B)
system.time(test2 <- reduxB(B))
test2
dim(test2$B)







rm(list = ls())
set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 50
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),5),
                sampling::inclusionprobabilities(runif(N),10),
                sampling::inclusionprobabilities(runif(N),15)),ncol = 3)
X <- PM(Pik)$PM
image(as(X,"sparseMatrix"))
pik <- PM(Pik)$P
dim(X)
order = 2
EPS = 1e-11

system.time(test1 <- reduxArma(X))
system.time(test2 <- reduxB(X))




rm(list = ls())
set.seed(1)
B <- as.matrix(rsparsematrix(4000,3000,density = 0.0001))
image(as(B,"sparseMatrix"))
system.time(test1 <- reduxArma(B))
dim(test1$B)
system.time(test2 <- reduxB(B))
dim(test2$B)
all(as.vector(test1$ind_row) == as.vector(test2$ind_row))