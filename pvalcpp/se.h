#pragma once 
#include <cmath>
#include <cstdint>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <iostream>
#include <cassert>
#include <queue>


#include "unif.h"
#include "mhstate.h"

class Estimate {
public:
  double mu_hat, se, var_hat, lambda;
  Estimate(const std::vector<double> & trajectory, const MHstate & );
};
