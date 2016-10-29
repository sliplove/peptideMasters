#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "wl.h"

void WLsimulator::wl_step(double phi, bool trace) {
	
	std::cout << "wl weights" <<  std::endl;
	for (auto & entry : weights_) {
		std::cout << entry << " ";
	}

	std::cout << std::endl;
	
	for (auto &h : hist_) { h = 0; }
	
	double lphi =  log(phi);
	int i = 0, idx;
	
	do {
		idx = mh_.step(weights_);
		hist_[idx]++ ;
		weights_[idx] -= lphi;
		i++;
		
		if (i % 1000 == 0) {
			if (trace) {
				std::cout << "iteration " << i << std::endl;
				for (const auto& w : weights_) {	
					std::cout << w << " ";
				}
				std::cout << std::endl;

				for (const auto& h : hist_) {
					std::cout << h << " ";
				}

				std::cout << std::endl;
			}		
		}

	} while (!hist_flatness() && (i < thr_)) ;

}

std::vector<double>  WLsimulator::wl_full(bool trace) {
	double phi = phi_begin_;
	if (trace) 
		std::cout << "wl step: phi = " << phi << std::endl;
	
	wl_step(phi, true);
	
	for (auto & entry : weights_) {
		std::cout << entry;
	}
	std::cout << std::endl;


	auto maxw =  std::max_element(weights_.begin(), weights_.end());

	for (int i = 0; i < weights_.size(); ++i) {
		weights_[i] -= *maxw;
	}
	std::cout << std::endl;

	while(phi > phi_end_) {
		phi = sqrt(phi);
  		wl_step(phi, true);


		auto maxw =  std::max_element(weights_.begin(), weights_.end());

		for (int i = 0; i < weights_.size(); ++i) {
			weights_[i] -= *maxw;
		}
		std::cout << std::endl;

	}

	return weights_;
}

void WLsimulator::print() {
	std::cout << std::endl; 
	std::cout << "------- Wang-Landau parameters -----------";
	std::cout << std::endl; 
	std::cout << "Initial phi value: " << phi_begin_ << std::endl << "Сrowning phi value:" << phi_end_ << std::endl;
	std::cout << "Number of iterations in wl step: " << thr_ << std::endl;

	std::cout << std::endl; 
	std::cout << "Weights: " << std::endl;

	for (auto & entry : weights_) 
		std::cout << entry << " "; 
	std::cout <<  std::endl;

	std::cout << std::endl; 
	std::cout << "Histogram: " << std::endl;	
	
	for (auto & entry : hist_) 
		std::cout << entry << " "; 
	std::cout <<  std::endl;

	mh_.get_state_().print_current_state_();
	
}


