#pragma once 
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include "peptide.h"


class MHstate {
private:	
	std::vector<double> current_state_;

public:
	MHstate(int length, double peptide_mass, pcg_extras::seed_seq_from<std::random_device> & rd) {
    	current_state_ = get_start_mass(rd, length, peptide_mass);
    }

    MHstate(const std::vector<double> & mass) {
        current_state_ = mass;
    }
    
    void set_current_state_(std::vector<double> state) {
        current_state_ = state; 
    }
   
    const std::vector<double> & get_current_state_() const {
        return current_state_;
     }

    void print_current_state_() const {
        for (auto & entry : current_state_) {
            printf("%e ", entry) ;
        }
        std::cout << std::endl;
    }

    void print_sum() const {
        double sum = 0;
        for (int i = 0; i < current_state_.size(); ++i)
        {
            sum += current_state_[i];
        }
        printf("sum=%e\n", sum );
    }

};

