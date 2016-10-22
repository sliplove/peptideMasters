
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

#include "peptide.h"
#include "unif.h"
#include "mhstate.h"

class WLsimulator {
private:
	std::vector<double> weights_;
	std::vector<int> hist_;
	double phi_begin_, phi_end_;
	int thr_;

public:
	WLsimulator(double phi_begin,double phi_end, int thr) 
	:
	phi_begin_(phi_begin),
	phi_end_(phi_end),
	thr_(thr) 
	{}

	std::vector<double> get_weights() const { return weights_; }

	bool hist_flatness() const {
		double mean_h, sum_h = 0;
		for (auto h : hist_)
			sum_h += h;
		mean_h = sum_h/hist_.size();
		for (const auto& h : hist_) {
			if ((h < 0.6 * mean_h) || (h < 20)) 	
				return false;
		}
		return true;
	}

	void wl_step(MHstate & mh, double , bool );
	std::vector<double>  wl_full(MHstate & , bool trace);
	
	void print();
};




