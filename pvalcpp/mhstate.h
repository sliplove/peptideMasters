#pragma once 
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "peptide.h"


class MHstate {
private:	
	// std::vector<double> weights_;
	std::vector<double> current_state_;

public:
    MHstate() { }
	MHstate(int length, double peptide_mass)
    {
    	// weights_.resize(int(max_score - min_score + 1));
    	current_state_ = get_start_mass(length, peptide_mass);
    }
    
    void set_current_state_(std::vector<double> state) {
        current_state_ = state; 
    }

    void drop_current_state_() {
        double l = current_state_.size();
        double mass = 0;
        for (int i = 0; i < l; ++i) {
            mass += current_state_[i];
        }
    	current_state_ = get_start_mass(l, mass);
    }
    
    const std::vector<double> & get_current_state_() const{
        return current_state_;
    }

    void print_current_state_() const{
        for (auto & entry : current_state_) {
            std::cout << entry << " ";
        }
        std::cout << std::endl;
    }

};

