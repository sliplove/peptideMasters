#include "main_functions.h"

std::vector<double> mh_weighted (const std::vector<double> exp_spectrum, 
                        const std::vector<std::vector<double>> &mat,
                      const std::vector<std::pair<unsigned, unsigned>> &rule,
                      std::vector<double> mass, double peptide_mass,const WLsimulator & wl, int N) {
    int i = 0;
    std::vector<double> v(N);
    double sscore, alpha;


    Peptide peptide(mat, rule, mass, peptide_mass);
    std::vector<double> nmass = update_mass_single_(mass, rule);    
    Peptide npeptide(mat, rule, nmass, peptide_mass);

    Psm psm, npsm;

    while (i < N) {
      nmass = update_mass_single_(mass, rule);
      peptide.set_spectrum(mass);
      npeptide.set_spectrum(nmass);
      psm.score_peak(exp_spectrum, peptide.get_spectrum_(), false);
      npsm.score_peak(exp_spectrum, npeptide.get_spectrum_(), false);
      
      double score_old =  std::min(psm.score_, MAX_SCORE);
      double score_new =  std::min(npsm.score_, MAX_SCORE);
      
      double w_old =  wl.get_single_weight(score_old);
      double w_new =  wl.get_single_weight(score_new);

      // std::cout << "score_old = " << score_old << " score_new = " << score_new << std::endl;
      // std::cout << "w_old = " << w_old << " w_new = " << w_new << std::endl;

      alpha = std::min(1.0, exp(w_new - w_old));  
      // std::cout << "alpha = " << alpha << std::endl;

      if (unif_rand() < alpha){
        peptide = npeptide;
        mass = nmass;
        // std::cout << "accepted" << std::endl;
        sscore = std::min(npsm.score_, MAX_SCORE);
      } else {
        sscore = std::min(psm.score_, MAX_SCORE);
      }
      v[i] = sscore;        
      i++;

    }
      return v;
  }

    

void hit_run(const std::vector<double> exp_spectrum, 
                        const std::vector<std::vector<double>> &mat,
                      const std::vector<std::pair<unsigned, unsigned>> &rule,
                      std::vector<double> mass, double peptide_mass,
                      const WLsimulator & wl, int step, int max_n, double eps_, double level) {

    double z = 1.96;
    std::vector<double> trajectory = mh_weighted(exp_spectrum, mat, rule, mass, peptide_mass, wl, step);
    std::vector<double> temp;
    double w, eps;

    do {
      Estimate est(trajectory, wl);
      w = 2*z*est.se;
      eps = w/sqrt(est.lambda);

      std::cout << "length = " << trajectory.size() << " mean = " << est.mu_hat << " sdev = " << est.se << 
                         " w = " << w << " eps = " << eps << " lambda = "  << est.lambda << std::endl; 

  
      temp = mh_weighted(exp_spectrum, mat, rule, mass, peptide_mass, wl, step);  
      trajectory.insert(std::end(trajectory), std::begin(temp), std::end(temp));

    } while ((eps > eps_) || (trajectory.size() < max_n) );
}


  


