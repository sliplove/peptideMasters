#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <iostream>
#include <cassert>
#include <queue>

#include "unif.h" 
#include "scorer.h"

const double MASS_PROTON = 1.00728;

void Scorer::score_peak(const Peptide & peptide , bool keep_both) {
  
  std::vector<double> spect = peptide.get_spectrum_();
  std::vector<size_t> sort_perm = peptide.get_sorting_permutation_();
  double nlp = peptide.get_peptide_mass_();


  std::vector<double> spectrum_(sort_perm.size());

  for (int i = 0; i < sort_perm.size(); ++i) {
    spectrum_[i]  = spect[sort_perm[i]];
  }
  
  score_ = 0;

  auto pmb = spectrum_.begin();
  auto pme = std::reverse_iterator<decltype(spectrum_.end())>(spectrum_.end());
  auto rpme = std::reverse_iterator<decltype(spectrum_.begin())>(spectrum_.begin());

  for (const auto& rp: exp_spectrum_) {
    double thr = rp - product_ion_thresh_ - MASS_PROTON;

    while (pmb != spectrum_.end() && *pmb <= thr) {
      pmb++;
    }
    
    if (pmb !=spectrum_.end()) {
      if (std::abs(*pmb + MASS_PROTON - rp) < product_ion_thresh_) {
        score_ += 1;
        // std::cout << *pmb << " " << rp << std::endl;
        continue;
      }
    }

    if (keep_both) {
      if (pmb == spectrum_.end())
        break;
      continue;
    }

    while (pme != rpme && nlp- *pme <= thr) {
      pme++;
    }

    if (pme != rpme) {
      if (std::abs(nlp - *pme + MASS_PROTON - rp) < product_ion_thresh_) {
        score_ += 1;
        // std::cout << nlp - *pme + MASS_PROTON << " " << rp << std::endl;
        continue;
      }
    }

    if (pmb == spectrum_.end() && pme == rpme) {
      break;
    }
    
  }
}
