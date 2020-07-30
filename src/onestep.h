#ifndef onestep_H
#define onestep_H

#include <RcppArmadillo.h>

arma::vec onestep(arma::mat B,arma::vec pik,double EPS=0.0000001);

#endif