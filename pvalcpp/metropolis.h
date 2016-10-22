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
#include "mhstate.h"
#include "unif.h"
#include "peptide.h"
#include "psm.h"

class Metropolis {
private:
	double peptide_mass_;
	MHstate mh_;
	std::vector<double> trajectory;

public:
	Metropolis(const MHstate & mh)
	:
	mh_(mh)
	{}

	std::vector<double>  mh_weighted (int);
	void hit_run( int , double , double );

};