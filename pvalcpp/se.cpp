#include "se.h"

  Estimate::Estimate(const std::vector<double> & trajectory, const WLsimulator & wl) {
  int n = trajectory.size();
  int b = floor(sqrt(n));
  int a = n - b + 1;
  double prob_const = 0;

  std::vector<double> elements(n);
  std::vector<double> cumsum(n + 1);
  std::vector<double> part_means(a);
  
  cumsum[0] = 0;

  for (int i = 0; i < n; ++i) {
    prob_const += exp((-1) *  wl.get_single_weight(trajectory[i]));;
  }
  
  prob_const = prob_const / trajectory.size();

  // std::cout << "prob_const = " << prob_const << std::endl;

  for (int i = 0; i < n; ++i) {
    elements[i] = (trajectory[i] >=  wl.get_max_score_())*exp((-1) * wl.get_single_weight(trajectory[i])) / prob_const;
    cumsum[i + 1] = cumsum[i] + elements[i];
  }

  // std::cout << "sum of elements = " << cumsum[n] << std::endl;

  for (int i = 0; i < a; ++i) {
    part_means[i] = (cumsum[i + b] - cumsum[i])/b;
  }

  // std::cout << std::endl;

  mu_hat = cumsum[n]/n;

  double tmp_sum = 0, tmp_sum1 = 0;
  for (int i = 0; i < a; ++i) {
    tmp_sum += (part_means[i] - mu_hat)*(part_means[i] - mu_hat);
  } 

  // std::cout << n*b*tmp_sum << std::endl;


  var_hat = n * b * tmp_sum / (a - 1) / a ;
  se = sqrt(var_hat / n);
  
  for (int i = 0; i < n; ++i) {
    tmp_sum1 += (elements[i] - mu_hat) * (elements[i] - mu_hat);
  }

  lambda = tmp_sum1 / n; 
}





