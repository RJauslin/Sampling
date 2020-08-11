library(BalancedSampling)
library(sampling)
library(WaveSampling)
library(pracma)
library(MASS)



rm(list=ls())




user <- c('setille\\switchdrive\\__PROJETS_DE_RECHERCHE\\',
          'eustachee\\switchdrive\\',
          'jauslinr\\switchdrive\\')[3]

# #---------LOADING PACKAGE---------#
# setwd(paste0('C:\\Users\\',user,'SpatioTempo\\SpotSampling\\'))
# devtools::load_all()
# 
# 



# #----------------DATA LIBELLULE---------------#


#-----Suisse
data      <- read.csv(file = paste0('C:\\Users\\',user,'Spatio-temporal\\simulations\\',"data_square_lib.csv"), sep = ',')
#-----Mitteland
data      <- read.csv(file = paste0('C:\\Users\\',user,'Spatio-temporal\\simulations\\',"mitteland_square_lib.csv"), sep = ',')





df       <- data
coord    <- as.matrix(df[,1:2])
interest <- df[,3]
N        <- nrow(df)


set.seed(1)
N_new <- 500
ind <- sample(seq(1,N,1),N_new)
N <- N_new
df <- df[ind,]
coord <- coord[ind,]





# #--------INCLUSION PROBABILITIES-----------#

t       <- 3
Pik     <- matrix(rep(0, t*N), ncol = t)
n1 <- N/3
n2 <- N/5
n3 <- N/7
Pik[,1] <- rep(n1/N,N)
Pik[,2] <- inclusionprobabilities(df[,4],n2)
Pik[,3] <- inclusionprobabilities(df[,3],n3)


system.time(test <- Orfs(Pik,arma = TRUE))
system.time(test <- Orfs(Pik,arma = FALSE))




#--------SIMULATIONs-----------#

nbSimu <- 5

RESULT.IB.SIMPLE  <- rep(0, nbSimu)
RESULT.IB.SPOT    <- rep(0, nbSimu)
RESULT.IB.XIAMEN  <- rep(0, nbSimu)

RESULT.sb.SIMPLE  <- rep(0, nbSimu)
RESULT.sb.SPOT    <- rep(0, nbSimu)
RESULT.sb.XIAMEN  <- rep(0, nbSimu)


for(j in 1:nbSimu){
  n <- colSums(Pik)
  
  #----------METHOD NEUCH---------------#
  
  #--Preselection
  # Pik.new      <- Preselection(Pik = Pik, coord = coord, L = 1)
  # Pik.remain   <- Pik.new[rowSums(Pik.new)>1e-7,]
  # coord.remain <- coord[rowSums(Pik.new)>1e-7,]
  
  #--Method SPOT
  A1                             <- Orfs(Pik.remain)
  A1.tot                         <- matrix(rep(0,N*t), nrow = N, ncol = t)
  A1.tot[rowSums(Pik.new)>1e-7,] <- A1
  
  #--Method ORFS
  A2                             <- Orsp(Pik.remain, coord.remain)
  A2.tot                         <- matrix(rep(0,N*t), nrow = N, ncol = t)
  A2.tot[rowSums(Pik.new)>1e-7,] <- A2
  
  #--Method ORSP
  A3                             <- Spot(Pik.remain, coord.remain)
  A3.tot                         <- matrix(rep(0,N*t), nrow = N, ncol = t)
  A3.tot[rowSums(Pik.new)>1e-7,] <- A3
  
  #----------METHOD XIAMEN---------------#
  setwd(paste0('C:\\Users\\',user,'Spatio-temporal\\simulations\\'))
  source("function.R")
  
  B     <- Xiamen(Pik, coord, 1)
  
  
  
  
  #-------PLOT
  SpreadPlot(A.tot, 1:4, IB, coord)
  SpreadPlot(rep(1,nrow(coord)), 1, IB, coord)
  
  SpreadPlot <- function(samples, time, criteria = c(IB,sb), coord){
    library(pracma)
    par(mfrow=c(ceil(length(time)/2),2))
    
    for(t in 1:length(time)){
      i <- time[t]
      plot(coord[,1], coord[,2], cex = 0.5)
      points(coord[samples[,i]==1,1],coord[samples[,i]==1,2], col = 'blue', pch = 19, cex = 0.5)
    }
  }
  
  
  
  #------- TEST ETALEMENT
  sb.simple <- rep(0,t)
  sb.Spot   <- rep(0,t)
  sb.Xiamen <- rep(0,t)
  
  IB.simple  <- rep(0,t)
  IB.Spot <- rep(0,t)
  IB.Xiamen <- rep(0,t)
  
  
  for(i in 1:t){
    #IB
    m.strat      <- wpik(X = coord, pik = Pik[,i])
    m.strat      <- m.strat-diag(diag(m.strat), nrow = nrow(m.strat), ncol = ncol(m.strat))
    IB.simple[i] <- IB(W = m.strat, s = A.tot[,i])
    IB.Spot[i]   <- IB(W = m.strat, s = srswor(n[i], N))
    IB.xiamen[i] <- IB(W = m.strat, s = B.tot[,i])
    
    #sb
    sb.simple[i] <- sb(Pik[,i], coord, (1:N)[A.tot[,i]==1])
    sb.Spot[i]   <- sb(Pik[,i], coord, (1:N)[srswor(n, N)==1])
    sb.xiamen[i] <- sb(Pik[,i], coord, (1:N)[B.tot[,i]==1])
  }
  
  RESULT.IB.SIMPLE[j]  <- mean(IB.simple)
  RESULT.IB.SPOT[j]    <- mean(IB.Spot)
  RESULT.IB.XIAMEN[j]  <- mean(IB.xiamen)
  
  RESULT.sb.SIMPLE[j]  <- mean(sb.simple)
  RESULT.sb.SPOT[j]    <- mean(sb.Spot)
  RESULT.sb.XIAMEN[j]  <- mean(sb.xiamen)
}






#--------------PLOT------------#

path <- paste0('C:\\Users\\',user,'Spatio-temporal\\simulations\\\\ggg_2019-LV03\\shp')
plotSwiss(X = data, Xs = data[A==1,], Region = c(2), Canton = TRUE, Commune = FALSE, path = path)
plotSwiss(X = data, Xs = NULL, Region = c(2), Canton = TRUE, Commune = FALSE, path = path)

