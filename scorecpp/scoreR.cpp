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


double score_peak(std::vector<double> &spectrum,std::vector<double> &peak_masses){
    double score = 0;
    double product_ion_thresh = 0.5;
    auto pmb = peak_masses.begin();
    auto pme = peak_masses.rbegin();

    // for (int i = 0; i < peak_masses.size(); ++i)
    // 	std::cout << peak_masses[i] << std::endl;
    
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

        while (pme != peak_masses.rend() && nlp_mass - *pme <= thr){
            pme++;
         	// std::cout << nlp_mass - *pme + MASS_PROTON << std::endl;
        }

        if (pme != peak_masses.rend()) {
            if (std::abs(nlp_mass - *pme + MASS_PROTON - rp) < product_ion_thresh) {
                score += 1;
                // std::cout << nlp_mass - *pme + MASS_PROTON << std::endl;
                continue;
            }
        }

        if (pmb == peak_masses.end() && pme == peak_masses.rend()) {
      	      break;
        }
    }

    return score;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
double get_score(DoubleVector s, DoubleVector t) {
  std::vector<double> spectrum = Rcpp::as<std::vector<double>>(s);
  std::vector<double> theoretic = Rcpp::as<std::vector<double>>(t);
  return score_peak(spectrum, theoretic);
}

