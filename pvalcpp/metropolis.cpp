#include "metropolis.h"
#include "unif.h"


std::vector<double>  Metropolis::update_spectrum_by_new_mass(const std::vector<double> &mass_) {
    
    unsigned id = floor(unif_rand() * rule_.size());
    unsigned beg = rule_[id].first;
    unsigned end = rule_[id].second;
    
    double beg_mass = mass_[beg], end_mass = mass_[end];
    double delta = unif(-beg_mass, end_mass);

    std::vector<double> mass (mass_);
    
    mass[beg] = beg_mass + delta;
    mass[end] = end_mass - delta;
    
    auto it_moved_plus = peptide.moved_plus_begin();
    auto it_moved_minus = peptide.moved_minus_begin();
    auto it_unmoved = peptide.unmoved_begin();

    for (size_t i : peptide.get_sorting_permutation_()) {
        if (mat_[i][beg] == mat_[i][end]) {
            *it_unmoved++ = i;
        } else if (mat_[i][beg]) {
            *it_moved_plus++ = i;
            peptide.increment_spectrum_(i, delta);
        } else {
            *it_moved_minus++ = i;
            peptide.increment_spectrum_(i, (-1)*delta);
        }
    }

    int sz = peptide.get_spectrum_().size();

    // Add sentineles
    *it_moved_plus = sz - 1;
    *it_moved_minus = sz - 1;
    *it_unmoved = sz - 1;

    // Rewind
    it_moved_plus = peptide.moved_plus_begin();
    it_moved_minus = peptide.moved_minus_begin();
    it_unmoved = peptide.unmoved_begin();

    for (size_t i = 0; i < peptide.get_sorting_permutation_().size(); ++i) {
        if ((peptide.get_spectrum_())[*it_unmoved] < (peptide.get_spectrum_())[*it_moved_plus] && 
            (peptide.get_spectrum_())[*it_unmoved] < (peptide.get_spectrum_())[*it_moved_minus]) {
            peptide.set_sorting_permutation_(i, *it_unmoved);
            it_unmoved++;
    } else if ((peptide.get_spectrum_())[*it_moved_minus] < (peptide.get_spectrum_())[*it_moved_plus]) {
        peptide.set_sorting_permutation_(i, *it_moved_minus);
        it_moved_minus++;
    } else {
        peptide.set_sorting_permutation_(i, *it_moved_plus);
        it_moved_plus++;
    }
}

return mass;

}

int Metropolis::step(const  std::vector<double> &weights) {
  peptide.clear(state_.get_current_state_());
  Peptide ppeptide(peptide);

  std::vector<double> nmass;

  double alpha, sscore;
  int i = 0, idx;
  
  nmass = update_spectrum_by_new_mass(state_.get_current_state_());

  scorer_.score_peak(ppeptide, false);
  double score_old =  std::min(scorer_.get_score_(), max_score_);
  scorer_.score_peak(peptide, false);
  double score_new =  std::min(scorer_.get_score_(), max_score_);
  
  double w_old =  weights[score_old + min_score_];
  double w_new =  weights[score_new + min_score_];


  alpha = std::min(1.0, exp(w_new - w_old));  
  
  if (unif_rand() < alpha){
    state_.set_current_state_(nmass);
    sscore = score_new;
  } else {
    peptide.set_spectrum_(peptide);
    sscore = score_old;
  }
  
  idx = (int)sscore - min_score_;

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