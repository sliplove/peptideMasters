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
#include "peptide.h"
#include "scorer.h"


class Metropolis {
private:
	std::vector<std::vector<double>> mat_;
	std::vector<std::pair<unsigned, unsigned>> rule_;
	double peptide_mass_;
	MHstate state_;
	std::vector<double> trajectory;
	Scorer scorer_;
	double min_score_, max_score_;

public:
	Metropolis() {}
	Metropolis(std::vector<std::vector<double>> mat,
		std::vector<std::pair<unsigned, unsigned>> rule,
		double peptide_mass, double min_score, double max_score,
		const Scorer & scorer)
	:
	mat_(mat),
	rule_(rule),
	peptide_mass_(peptide_mass),
	min_score_(min_score),
	max_score_(max_score),	
	scorer_(scorer),
	state_(mat[0].size(), peptide_mass) {}

	const double & get_min_score_() const { return min_score_; }
	const double & get_max_score_() const { return max_score_; }
	double get_range() const { 
		return (int)(max_score_ - min_score_ + 1);
	}
	const MHstate & get_state_() const { return state_; }
	MHstate &get_state_() { return state_; }
	
	int step(const  std::vector<double> &);
	std::vector<double>  mh_weighted (int, const std::vector<double> &);
	void hit_run( int , double , double , const std::vector<double> &);



};