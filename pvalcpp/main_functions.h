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

#include "se.h"
#include "wl.h"
#include "unif.h"
#include "peptide.h"
#include "psm.h"


std::vector<double> mh_weighted (const std::vector<double> , 
                      const std::vector<std::vector<double>> &,
                      const std::vector<std::pair<unsigned, unsigned>> & ,
                      std::vector<double> , double , const WLsimulator & , int ) 
;
void hit_run(const std::vector<double> , 
                      const std::vector<std::vector<double>> &,
                      const std::vector<std::pair<unsigned, unsigned>> &, 
                      std::vector<double> mass, double ,
                      const WLsimulator & , int , int ,  double , double );