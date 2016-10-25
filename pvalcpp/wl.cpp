#include "wl.h"

void WLsimulator::wl_step(MHstate & mh, double phi, bool trace) {
	weights_ = mh.get_weights_();

	std::cout << "wl weights" <<  std::endl;
	for (auto & entry : weights_) {
		std::cout << entry << " ";
	}

	std::cout << std::endl;
	hist_.resize(weights_.size());
	
	for (auto & h : hist_) { h = 0; }

	double lphi =  log(phi);
	int i = 0, idx;
	
	do {
		idx = mh.step();
			// std::cout << "idx = " << idx << std::endl;
		hist_[idx] = hist_[idx] + 1;
		weights_[idx] = weights_[idx] - lphi;
		i++;
		
		mh.set_weights(weights_);

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

std::vector<double>  WLsimulator::wl_full(MHstate & mh, bool trace) {
	double phi = phi_begin_;
	if (trace) 
		std::cout << "wl step: phi = " << phi << std::endl;
	
	mh.print_current_state_();
	mh.drop_current_state_();
	std::cout << "drop current state" << std::endl; 
	mh.print_current_state_();
	wl_step(mh, phi, true);
	std::cout << "mh weights_" << std::endl; 
	mh.print_weights_();
	
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
  		// std::cout << "phi = " << phi << std::endl;
		mh.drop_current_state_();
  		// mh.print_current_state_(	);
		mh.print_weights_();
		wl_step(mh, phi, true);


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
}


