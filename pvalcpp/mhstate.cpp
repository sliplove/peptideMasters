#include "mhstate.h"

int MHstate::step() {
	Peptide peptide(mat_, rule_, current_state_, peptide_mass_);		
	Peptide npeptide(mat_, rule_, current_state_, peptide_mass_);

	std::vector<double> nmass;

	Psm psm, npsm;

	double alpha, sscore;
	int i = 0, idx;
	
	nmass = update_spectrum_by_new_mass(current_state_, rule_, npeptide);

	psm.score_peak(exp_spectrum_, peptide, false);
	npsm.score_peak(exp_spectrum_, npeptide, false);
	
	double score_old =  std::min(psm.score_, max_score_);
	double score_new =  std::min(npsm.score_, max_score_);
	
	double w_old =  get_single_weight(score_old);
	double w_new =  get_single_weight(score_new);

		// std::cout << "score_old = " << score_old << " score_new = " << score_new << std::endl;
		// std::cout << "w_old = " << w_old << " w_new = " << w_new << std::endl;

	alpha = std::min(1.0, exp(w_new - w_old));	
		// std::cout << "alpha = " << alpha << std::endl;
	
	if (unif_rand() < alpha){
		current_state_ = nmass;
		peptide.copy_spectrum_(npeptide);
		sscore = std::min(npsm.score_, max_score_);
	} else {
		sscore = std::min(psm.score_, max_score_);
	}
	
	idx = (int)sscore - min_score_;
	
	
	return idx;
} 