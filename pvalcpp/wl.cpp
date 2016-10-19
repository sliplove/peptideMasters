#include "wl.h"
#include "psm.h"


void WLsimulator::wl_step(const std::vector<double> exp_spectrum, 
					const std::vector<std::vector<double>> &mat,
                      const std::vector<std::pair<unsigned, unsigned>> &rule,
                      std::vector<double> mass, double peptide_mass, double phi, bool trace) {
		Peptide peptide(mat, rule, mass, peptide_mass);		
		Peptide npeptide(mat, rule, mass, peptide_mass);

		std::vector<double> nmass;

		Psm psm, npsm;

		double lphi =  log(phi);
		double alpha, sscore;
		int i = 0, idx;
		
		for (auto & h : hist) { h = 0; }
		
		do {
			npeptide.clear(mass);
			// Peptide npeptide(mat, rule, mass, peptide_mass);
			nmass = update_spectrum_by_new_mass(mass, rule, npeptide);

			psm.score_peak(exp_spectrum, peptide, false);
			npsm.score_peak(exp_spectrum, npeptide, false);
			
			double score_old =  std::min(psm.score_, MAX_SCORE);
			double score_new =  std::min(npsm.score_, MAX_SCORE);
			
			double w_old =  get_single_weight(score_old);
			double w_new =  get_single_weight(score_new);

			// std::cout << "score_old = " << score_old << " score_new = " << score_new << std::endl;
			// std::cout << "w_old = " << w_old << " w_new = " << w_new << std::endl;

			alpha = std::min(1.0, exp(w_new - w_old));	
			// std::cout << "alpha = " << alpha << std::endl;

			if (unif_rand() < alpha){
				mass = nmass;
				peptide.copy_spectrum_(npeptide);
				sscore = std::min(npsm.score_, MAX_SCORE);
			} else {
				sscore = std::min(psm.score_, MAX_SCORE);
			}
				
			idx = (int)sscore - min_score_;

    		hist[idx] = hist[idx] + 1;
    		weights[idx] = weights[idx] - lphi;
    		i++;

			// print stuff 
	    	if (i % 1000 == 0) {
	    		if (trace) {
	    			std::cout << "iteration " << i << std::endl;
	    			for (const auto& w : weights) {	
	    				std::cout << w << " ";
	    			}
	    			std::cout << std::endl;

	    			for (const auto& h : hist) {
	    				std::cout << h << " ";
	    			}

	    			std::cout << std::endl;
	    		}		
	     	}
		} while (!hist_flatness() && (i < thr_)) ;
}

std::vector<double>  WLsimulator::wl_full (const std::vector<double> exp_spectrum,
	 const std::vector<std::vector<double>> &mat,
	 const std::vector<std::pair<unsigned, unsigned>> &rule,
                      double peptide_mass,  bool trace) {

	int nrow = mat.size() - 1;
    int ncol = mat[0].size();

	double phi = phi_begin_;
	if (trace) 
		std::cout << "wl step: phi = " << phi << std::endl;
	
	
	std::vector<double> st = get_start_mass(ncol, peptide_mass);
	// peptide.print();
	wl_step(exp_spectrum, mat, rule, st, peptide_mass, phi, true);

	auto maxw =  std::max_element(std::begin(weights), std::end(weights));

	for (int i = 0; i < weights.size(); ++i) {
		weights[i] -= *maxw;
	}
	std::cout << std::endl;

	while(phi > phi_end_) {
  		phi = sqrt(phi);
  		// std::cout << "phi = " << phi << std::endl;
  		std::vector<double> st = get_start_mass(ncol, peptide_mass);
		wl_step(exp_spectrum, mat, rule, st, peptide_mass, phi, true);

  		auto maxw =  std::max_element(weights.begin(), weights.end());

		for (int i = 0; i < weights.size(); ++i) {
			weights[i] -= *maxw;
		}
		std::cout << std::endl;

	}

	return weights;
}

void WLsimulator::print() {
	std::cout << std::endl; 
	std::cout << "------- Wang-Landau parameters -----------";
    std::cout << std::endl; 
	    std::cout << "Min score: " <<  min_score_ << std::endl << "Max score: " <<  max_score_<< std::endl;
	    std::cout << "Initial phi value: " << phi_begin_ << std::endl << "Сrowning phi value:" << phi_end_ << std::endl;
	    std::cout << "Number of iterations in wl step: " << thr_ << std::endl;

		std::cout << std::endl; 
	    std::cout << "Weights: " << std::endl;

	    for (auto & entry : weights) 
	        std::cout << entry << " "; 
	        std::cout <<  std::endl;

	    std::cout << std::endl; 
	    std::cout << "Histogram: " << std::endl;	
	    
	    for (auto & entry : hist) 
	        std::cout << entry << " "; 
	        std::cout <<  std::endl;
}


	