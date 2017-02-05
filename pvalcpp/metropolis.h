#pragma once 
#include <vector>
#include <fstream>

#include "mhstate.h"
#include "peptide.h"
#include "scorer.h"
#include "../include/pcg_random.hpp"

using namespace scoring;

class Metropolis {
private:
	std::vector<std::vector<double>> mat_;
	std::vector<std::pair<unsigned, unsigned>> rule_;
	double peptide_mass_;
	MHstate state_;
	std::vector<double> trajectory;
	SPCScorer scorer_;
	Spectrum spectrum_;
	double min_score_, max_score_;
	Peptide peptide;

public:
	Metropolis(std::vector<std::vector<double>> mat,
		std::vector<std::pair<unsigned, unsigned>> rule,
		double peptide_mass, double min_score, double max_score,
		const MHstate &state,
		const Peptide &peptide,
		const Spectrum & spectrum, const SPCScorer & scorer)
	:
	mat_(mat),
	rule_(rule),
	peptide_mass_(peptide_mass),
	min_score_(min_score),
	max_score_(max_score),	
	state_(state),
	peptide(peptide),
	spectrum_(spectrum), 
	scorer_(scorer) {}

	const double & get_min_score_() const { return min_score_; }
	const double & get_max_score_() const { return max_score_; }
	double get_range() const { 
		return (int)(max_score_ - min_score_ + 1);
	}

	const MHstate & get_state_() const { return state_; }
	
	int step(pcg_extras::seed_seq_from<std::random_device> & ,
	 const  std::vector<double> &);
	std::vector<double>  mh_weighted (pcg_extras::seed_seq_from<std::random_device> & , int, const std::vector<double> &);
	void hit_run(pcg_extras::seed_seq_from<std::random_device> &,
	 int , double , double , const std::vector<double> &);
	std::vector<double>  update_spectrum_by_new_mass(pcg_extras::seed_seq_from<std::random_device> &,
	 const std::vector<double> &);

};