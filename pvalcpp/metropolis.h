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

class Metropolis {
private:
	std::vector<double> exp_spectrum_;
	std::vector<std::vector<double>> mat_;
	std::vector<std::pair<unsigned, unsigned>> rule_;
	double peptide_mass_;
	WLsimulator wl_;
	std::vector<double> trajectory;

public:
	Metropolis(const std::vector<double> exp_spectrum, 
               const std::vector<std::vector<double>> &mat,
               const std::vector<std::pair<unsigned, unsigned>> &rule,
               double peptide_mass, const WLsimulator & wl)
	:
	exp_spectrum_(exp_spectrum),
	mat_(mat),
	rule_(rule),
    peptide_mass_(peptide_mass),
    wl_(wl)
    {}

	std::vector<double>  mh_weighted (std::vector<double> &, int );
	void hit_run(std::vector<double> , int , double , double );

};