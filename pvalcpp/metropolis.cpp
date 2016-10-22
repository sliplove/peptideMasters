#include "metropolis.h"

std::vector<double>  Metropolis::mh_weighted (int N) {
  int i = 0;
  std::vector<double> v(N);
  double idx;

  while (i < N) {
    idx = mh_.step();
    v[i] = idx + mh_.get_min_score_();        
    i++;
  }

  return v;
}



void Metropolis::hit_run(int step, double eps_, double level) {

  double z = 1.96;
  std::vector<double> trajectory = mh_weighted(step);
  std::vector<double> temp;
  double w, eps;

  do {
    Estimate est(trajectory, mh_);
    w = 2*z*est.se;
    eps = w/sqrt(est.lambda);

    std::cout << "length = " << trajectory.size() << " mean = " << est.mu_hat << " sdev = " << est.se << 
    " w = " << w << " eps = " << eps << " lambda = "  << est.lambda << std::endl; 

    
    temp = mh_weighted(step);  
    trajectory.insert(std::end(trajectory), std::begin(temp), std::end(temp));

  } while (eps > eps_);
}





