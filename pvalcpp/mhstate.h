
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
#include "psm.h"

class MHstate {
private:	
	std::vector<double> exp_spectrum_;
	std::vector<std::vector<double>> mat_;
	std::vector<std::pair<unsigned, unsigned>> rule_;
	double peptide_mass_;
	std::vector<double> weights_;
	std::vector<double> current_state_;
	double min_score_, max_score_;

public:
	MHstate(const std::vector<double> &exp_spectrum,
     const std::vector<std::vector<double>> &mat,
     const std::vector<std::pair<unsigned, unsigned>> &rule,
     double peptide_mass, double min_score, double max_score)
	:
	exp_spectrum_(exp_spectrum),
	mat_(mat),
	rule_(rule),
    peptide_mass_(peptide_mass),
    min_score_(min_score),
    max_score_(max_score)
    {
    	weights_.resize(int(max_score_ - min_score_ + 1));
    	current_state_ = get_start_mass(mat_[0].size(), peptide_mass_);
    }
    void set_weights(std::vector<double> weights) {
    	for (int i = 0; i < weights_.size(); ++i) {
    		weights_[i] = weights[i];
    	}

    }
    void drop_current_state_() {
    	current_state_ = get_start_mass(mat_[0].size(), peptide_mass_);
    }
    std::vector<double> get_weights_() const { return weights_; }
    double get_single_weight(double score) const {
    	return weights_[int(score - min_score_)]; 
    }
    std::vector<double> get_current_state_() const { return current_state_; }
    double get_min_score_() const {return min_score_;}
    double get_max_score_() const {return max_score_;}

    int step();
    void print_current_state_() {
        for (auto & entry : current_state_) {
            std::cout << entry << " ";
        }
        std::cout << std::endl;
    }
    void print_weights_() {
        for (auto & entry : weights_) {
            std::cout << entry << " ";
        }
        std::cout << std::endl;
    }


};

