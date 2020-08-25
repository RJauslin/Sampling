library(BalancedSampling)
library(sampling)
library(WaveSampling)
library(pracma)
library(MASS)
library(gridExtra)


rm(list=ls())




user <- c('setille\\switchdrive\\__PROJETS_DE_RECHERCHE\\',
          'eustachee\\switchdrive\\',
          'jauslinr\\switchdrive\\')[3]

#---------LOADING PACKAGE---------#
# setwd(paste0('C:\\Users\\',user,'\\Spatio-temporal\\package\\SpotSampling\\'))
# devtools::load_all()

#---------LOADING OTHER FUNCTIONS---------#
# setwd(paste0('C:\\Users\\',user,'Spatio-temporal\\simulations\\'))
source("C:\\Users\\jauslinr\\switchdrive\\Spatio-temporal\\simulations\\function.R")

# #----------------DATA LIBELLULE---------------#
# setwd(paste0('C:\\Users\\',user,'Spatio-temporal\\simulations\\'))

#-----Suisse
# data      <- read.csv(file = "data_square_lib.csv", sep = ',')
#-----Mitteland
#data      <- read.csv(file = "mitteland_square_lib.csv", sep = ',')
#-----Neuch Vaud Gruyere
data <- read.csv(file = "C:\\Users\\jauslinr\\switchdrive\\Spatio-temporal\\simulations\\Neuch_Vaud_Gruyere_2.csv", sep = ',')

df       <- data[sample(1:1407,200),]
coord    <- as.matrix(df[,1:2])
N        <- nrow(df)







# #--------INCLUSION PROBABILITIES-----------#
t       <- 3
pik     <- matrix(rep(0, t*N), ncol = t)
pik[,1] <- rep(20/N,N)
pik[,2] <- rep(20/N,N)
pik[,3] <- rep(20/N,N)
#pik[,2] <- inclusionprobabilities(df[,4],120)
# pik[,1] <- inclusionprobabilities(df[,3],16)
# pik[,2] <- inclusionprobabilities(df[,3],18)
# pik[,3] <- inclusionprobabilities(df[,3],20)


# S <- matrix(rep(0,N*t), nrow = N, ncol = t)
#----------METHOD NEUCH---------------#
#--Preselection
# EPS <- 1e-7
# pik_new      <- Preselection(pik, coord, L = 1)
# TEST         <- rowSums(pik_new)>EPS
# pik_remain   <- pik_new[TEST,]
# coord_remain <- coord[TEST,]
#
# system.time(res  <- Spot(pik_remain, coord_remain, comment = TRUE)) # -0.2715581 -0.3267983 -0.3212758 (IB)
# S[TEST,]  <- res
#
# system.time(res2 <- Orfs(pik, comment = TRUE))
# system.time(res3 <- Orsp(pik, coord, comment = TRUE)) # -0.2698679 -0.3459062 -0.3413718
# system.time(res4 <- Xiamen(pik, coord)) # -0.2905476 -0.2422639 -0.2519039








#--------SIMULATIONs-----------#

nbSimu <- 10
# method <- c('Spot', 'Orsp', 'Xiamen', 'Orfs')
# 
# 
# for(i in 1:4){
#   m <- method[i]
#   simu(pik, coord, nbSimu, m, EPS = 1e-7)
# }

# simu(pik, coord, nbSimu, m = 'Orfs', EPS = 1e-7)
# t      <- ncol(pik)
# N      <- nrow(pik)
# IB <- matrix(rep(0, nbSimu*t), ncol = t)
# sb <- matrix(rep(0, nbSimu*t), ncol = t)


set.seed(8)
EPS = 1e-7
n <- colSums(pik)
pik_new      <- Preselection(pik, coord, L = 1)
TEST         <- rowSums(pik_new)>EPS
pik_remain   <- pik_new[TEST,]
coord_remain <- coord[TEST,]
S <- Orfs(pik = pik_remain,comment = TRUE)
S
