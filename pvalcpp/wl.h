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
#include "metropolis.h"

class WLsimulator {
private:
	std::vector<double> weights_;
	std::vector<int> hist_;
	double phi_begin_, phi_end_;
	int thr_;
	Metropolis mh_;

public:
	WLsimulator(const Metropolis &mh, double phi_begin,
						double phi_end, int thr) 
	:
	phi_begin_(phi_begin),
	phi_end_(phi_end),
	thr_(thr),
	mh_(mh)
	{
		weights_.resize(mh_.get_range());
		hist_.resize(mh_.get_range());
	}

	std::vector<double> get_weights() const { return weights_; }
	double get_single_weight(int idx) const {return weights_[idx];}

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

	void wl_step(double , bool );
	std::vector<double>  wl_full(bool);
	
	
	void print();
};




