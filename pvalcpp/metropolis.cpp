#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "metropolis.h"
#include "unif.h"
#include "se.h"

std::vector<double>  Metropolis::update_spectrum_by_new_mass(pcg_extras::seed_seq_from<std::random_device>& rd, const std::vector<double> &mass_) {
    
    // std::cout << "rule_size = " << rule_.size() << std::endl;
    unsigned id = unif_int(rd, 0,  rule_.size() - 1);
    unsigned beg = rule_[id].first;
    unsigned end = rule_[id].second;
    
    double beg_mass = mass_[beg], end_mass = mass_[end];
    double delta = unif_real(rd, -beg_mass, end_mass);

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
    peptide.set_sorted_spectrum();

return mass;

}


int Metropolis::step(pcg_extras::seed_seq_from<std::random_device>& rd,
                     const  std::vector<double> &weights) {
   peptide.clear(state_.get_current_state_());
  // state_.print_current_state_();

  // std::cout << "First Peptide" << std::endl;
  // peptide.print(); 
  Peptide ppeptide(peptide);

  std::vector<double> nmass;

  double alpha, sscore;
  int i = 0, idx;
  

  nmass = update_spectrum_by_new_mass(rd, state_.get_current_state_());

  // double score_old = std::min(scorer_.score_peak(spectrum_, ppeptide, false), max_score_);
  // double score_new = std::min(scorer_.score_peak(spectrum_, peptide, false), max_score_);

  double score_old = std::min(scorer_.score(spectrum_, ppeptide, false), max_score_);
  double score_new = std::min(scorer_.score(spectrum_, peptide, false), max_score_);

  // std::cout << "old score : " << score_old << " : " << score_new << " new score : " << nscore_old << " : " << nscore_new << std::endl;
  double w_old =  weights[score_old - min_score_];
  double w_new =  weights[score_new - min_score_];
      

  alpha = std::min(1.0, exp(w_new - w_old));  
  double temp = unif_real(rd, 0, 1);
  if (temp < alpha){
    state_.set_current_state_(nmass);
    sscore = score_new;
    // std::cout << "accepted" << temp << ":" << alpha << std::endl;

  } else {
    // peptide.set_spectrum_(ppeptide);
    sscore = score_old;
  }
  
  idx = (int)sscore - min_score_;

  // state_.print_current_state_();
  // state_.print_sum();
  // peptide.clear(state_.get_current_state_());
  // peptide.print();

  // std::cout << "score=" << sscore << std::endl;

  return idx;
} 


std::vector<double>  Metropolis::mh_weighted (pcg_extras::seed_seq_from<std::random_device>& rd, int N, const std::vector<double> &weights) {
  int i = 0;
  std::vector<double> v(N);
  double idx;
  // state_.drop_current_state_();

  while (i < N) {
    idx = step(rd, weights);
    v[i] = idx + min_score_;        
    i++;
  }

  return v;
}


void Metropolis::hit_run(pcg_extras::seed_seq_from<std::random_device> & rd,
               int step, double eps_, double level, const std::vector<double> & weights) {

  double z = 1.96;
  std::vector<double> trajectory = mh_weighted(rd, step, weights);
  std::vector<double> temp;
  double w, eps;
  // Estimate est(trajectory, weights, min_score_, max_score_);  
  // w = 2*z*est.se;
  // eps = w/sqrt(est.lambda);

    
  // std::cout << "length = " << trajectory.size() << " mean = " << est.mu_hat << " sdev = " << est.se << 
  // " w = " << w << " eps = " << eps << " lambda = "  << est.lambda <<
  //  " lower = " << est.mu_hat - w/2 << " upper = " << est.mu_hat + w/2  << std::endl;

  do {
    Estimate est(trajectory, weights, min_score_, max_score_);
    w = 2*z*est.se;
    eps = w/sqrt(est.lambda);

    
    std::cout << "length = " << trajectory.size() << " mean = " << est.mu_hat << " sdev = " << est.se << 
    " w = " << w << " eps = " << eps << " lambda = "  << est.lambda <<
     " lower = " << est.mu_hat - w/2 << " upper = " << est.mu_hat + w/2  << std::endl;
    
    state_.print_current_state_();
    
    temp = mh_weighted(rd, step, weights);  
    trajectory.insert(std::end(trajectory), std::begin(temp), std::end(temp));

  } while (eps > eps_);
  // for (int i = 0; i < trajectory.size(); ++i)
  // {
  //   std::cout << trajectory[i] << " ";
  // }
} 