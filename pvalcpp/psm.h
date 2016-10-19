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
#include "peptide.h"

class Psm {
public:
	double score_;
	void score_peak (const std::vector<double> &, 
		const Peptide & peptide ,  bool );
};