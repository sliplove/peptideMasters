#include "metropolis.h"
#include "unif.h"

int Metropolis::step(const  std::vector<double> &weights) {
  Peptide peptide(mat_, rule_, state_.get_current_state_(), peptide_mass_);    
  Peptide npeptide(mat_, rule_, state_.get_current_state_(), peptide_mass_);

  std::vector<double> nmass;

  double alpha, sscore;
  int i = 0, idx;
  
  nmass = update_spectrum_by_new_mass(state_.get_current_state_(), rule_, npeptide);

  scorer_.score_peak(peptide, false);
  double score_old =  std::min(scorer_.get_score_(), max_score_);
  scorer_.score_peak(npeptide, false);
  double score_new =  std::min(scorer_.get_score_(), max_score_);
  
  double w_old =  weights[score_old + min_score_];
  double w_new =  weights[score_new + min_score_];

  std::cout << "score_old = " << score_old << " score_new = " << score_new << std::endl;
  std::cout << "w_old = " << w_old << " w_new = " << w_new << std::endl;

  alpha = std::min(1.0, exp(w_new - w_old));  
    std::cout << "alpha = " << alpha << std::endl;
  
  if (unif_rand() < alpha){
    state_.set_current_state_(nmass);
    peptide.set_spectrum_(npeptide);
    sscore = score_new;
  } else {
    sscore = score_old;
  }
  
  idx = (int)sscore - min_score_;
  std::cout << "idx = " << idx << std::endl;

  return idx;
} 

std::vector<double>  Metropolis::mh_weighted (int N, const std::vector<double> &weights) {
  int i = 0;
  std::vector<double> v(N);
  double idx;

  while (i < N) {
    idx = step(weights);
    v[i] = idx + min_score_;        
    i++;
  }

  return v;
}



void Metropolis::hit_run(int step, double eps_, double level,
 const std::vector<double> & weights) {

  double z = 1.96;
  std::vector<double> trajectory = mh_weighted(step, weights);
  std::vector<double> temp;
  double w, eps;

  do {
    Estimate est(trajectory, weights, min_score_, max_score_);
    w = 2*z*est.se;
    eps = w/sqrt(est.lambda);

    std::cout << "length = " << trajectory.size() << " mean = " << est.mu_hat << " sdev = " << est.se << 
    " w = " << w << " eps = " << eps << " lambda = "  << est.lambda << std::endl; 

    
    temp = mh_weighted(step, weights);  
    trajectory.insert(std::end(trajectory), std::begin(temp), std::end(temp));

  } while (eps > eps_);
}





