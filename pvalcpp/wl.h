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

class WLsimulator {
private:
	std::vector<double> weights;
	std::vector<int> hist;
	double min_score_, max_score_, phi_begin_, phi_end_;
	int thr_;

public:
	WLsimulator(double min_score, double max_score, double phi_begin, 
					double phi_end, int thr) : min_score_{min_score}, max_score_{max_score},
							phi_begin_{phi_begin}, phi_end_{phi_end}, thr_{thr} {
		weights.resize(int(max_score_ - min_score_ + 1));
		hist.resize(int(max_score_ - min_score_ + 1));
		}

	std::vector<double> get_weights() const { return weights; }
	double get_single_weight(double score) const {return weights[int(score - min_score_)]; }

	bool hist_flatness() {
		double mean_h, sum_h = 0;
		for (auto h : hist)
    		sum_h += h;
    	mean_h = sum_h/hist.size();
		for (const auto& h : hist) {
			if ((h < 0.6 * mean_h) || (h < 20)) 	
					return false;
		}
		return true;
	}

	void wl_step(const std::vector<double> , 
		const std::vector<std::vector<double>> &,
        const std::vector<std::pair<unsigned, unsigned>> &,
        std::vector<double> , double , double , bool );

	std::vector<double>  wl_full(const std::vector<double> exp_spectrum,
	 const std::vector<std::vector<double>> &mat,
	 const std::vector<std::pair<unsigned, unsigned>> &rule,
     double peptide_mass,  bool trace);
	
	void print();
};




