#pragma once 
#include <vector>

#include "peptide.h"

class Scorer {
private:
	std::vector<double> exp_spectrum_;
	double score_;
	double product_ion_thresh_;
public:
	Scorer();
	Scorer(const std::vector<double> &exp_spectrum, 
		double product_ion_thresh)
	:
	exp_spectrum_(exp_spectrum),
	product_ion_thresh_(product_ion_thresh) {}

	double get_score_() const {return score_;}
	std::vector<double>  get_exp_spectrum_() const {
		return exp_spectrum_;
	}
	void score_peak (const Peptide &,  bool);
};