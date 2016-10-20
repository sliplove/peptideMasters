#include "metropolis.h"

std::vector<double>  Metropolis::mh_weighted (std::vector<double> & mass, int N) {
    int i = 0;
    std::vector<double> v(N);
    double sscore, alpha;


    Peptide peptide(mat_, rule_, mass, peptide_mass_);
    Peptide npeptide(mat_, rule_, mass, peptide_mass_);
    
    std::vector<double> nmass;
    Psm psm, npsm;

    while (i < N) {
      npeptide.clear(mass);
      nmass = update_spectrum_by_new_mass(mass, rule_, npeptide);
      psm.score_peak(exp_spectrum_, peptide, false);
      npsm.score_peak(exp_spectrum_, npeptide, false);
      
      double score_old =  std::min(psm.score_, wl_.get_max_score_());
      double score_new =  std::min(npsm.score_, wl_.get_max_score_());
      
      double w_old =  wl_.get_single_weight(score_old);
      double w_new =  wl_.get_single_weight(score_new);

      // std::cout << "score_old = " << score_old << " score_new = " << score_new << std::endl;
      // std::cout << "w_old = " << w_old << " w_new = " << w_new << std::endl;

      alpha = std::min(1.0, exp(w_new - w_old));  
      // std::cout << "alpha = " << alpha << std::endl;

      if (unif_rand() < alpha){
        peptide.copy_spectrum_(npeptide);
        mass = nmass;
        // std::cout << "accepted" << std::endl;
        sscore = std::min(npsm.score_, wl_.get_max_score_());
      } else {
        sscore = std::min(psm.score_,  wl_.get_max_score_());
      }
      v[i] = sscore;        
      i++;

    }
      return v;
  }

    

void Metropolis::hit_run(std::vector<double> mass, int step, double eps_, double level) {

    double z = 1.96;
    std::vector<double> trajectory = mh_weighted(mass, step);
    std::vector<double> temp;
    double w, eps;

    do {
      Estimate est(trajectory, wl_);
      w = 2*z*est.se;
      eps = w/sqrt(est.lambda);

      std::cout << "length = " << trajectory.size() << " mean = " << est.mu_hat << " sdev = " << est.se << 
                         " w = " << w << " eps = " << eps << " lambda = "  << est.lambda << std::endl; 

  
      temp = mh_weighted(mass, step);  
      trajectory.insert(std::end(trajectory), std::begin(temp), std::end(temp));

    } while (eps > eps_);
}


  


