#pragma once 
#include <vector>

#include "unif.h"
#include "mhstate.h"

class Estimate {
public:
  double mu_hat, se, var_hat, lambda;
  Estimate(const std::vector<double> & , const std::vector<double> &, double , double);
};
