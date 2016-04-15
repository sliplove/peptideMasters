#include <Rcpp.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace Rcpp;


const double nlp_mass = 1007.65;
const double MASS_PROTON = 1.00728;
const bool keep_both = false;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double score_peak(DoubleVector spectrum, DoubleVector peak_masses) {
  std::sort(peak_masses.begin(), peak_masses.end(), [](double a, double b) { return a < b; });

  double score = 0;
  double product_ion_thresh = 0.5;
  auto pmb = peak_masses.begin();
  auto pme = std::reverse_iterator<decltype(peak_masses.end())>(peak_masses.end());
  auto rpme = std::reverse_iterator<decltype(peak_masses.begin())>(peak_masses.begin());

  for (const auto& rp: spectrum) {
    double thr = rp - product_ion_thresh - MASS_PROTON;
            
    while (pmb != peak_masses.end() && *pmb <= thr) {
      pmb++;
      // std::cout << *pmb << std::endl;
    }

    if (pmb != peak_masses.end()) {
      if (std::abs(*pmb + MASS_PROTON - rp) < product_ion_thresh) {
        score += 1;
        // std::cout << *pmb << std::endl;
        continue;
      }
    }

    if (keep_both) {
      if (pmb == peak_masses.end())
        break;
      continue;
    }

    while (pme != rpme && nlp_mass - *pme <= thr) {
      pme++;
      // std::cout << nlp_mass - *pme + MASS_PROTON << std::endl;
    }

    if (pme != rpme) {
      if (std::abs(nlp_mass - *pme + MASS_PROTON - rp) < product_ion_thresh) {
        score += 1;
        // std::cout << nlp_mass - *pme + MASS_PROTON << std::endl;
        continue;
      }
    }

    if (pmb == peak_masses.end() && pme == rpme) {
      break;
    }
  }

  return score;
}
