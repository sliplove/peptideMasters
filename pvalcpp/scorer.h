#pragma once 
#include <vector>

#include "peptide.h"

class Scorer {
private:
	std::vector<double> exp_spectrum_;
	double score_;
public:
	Scorer(const std::vector<double> &exp_spectrum)
	:
	exp_spectrum_(exp_spectrum) {}

	double get_score_() const {return score_;}
	std::vector<double>  get_exp_spectrum_() const {return exp_spectrum_;}
	void score_peak (const Peptide &,  bool );
};