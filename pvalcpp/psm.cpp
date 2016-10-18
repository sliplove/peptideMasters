#include "psm.h"

void Psm::score_peak(const std::vector<double> & exp_spectrum_, 
		const std::vector<double> & spectrum_, bool keep_both) {
  score_ = 0;

  auto pmb = spectrum_.begin();
  auto pme = std::reverse_iterator<decltype(spectrum_.end())>(spectrum_.end());
  auto rpme = std::reverse_iterator<decltype(spectrum_.begin())>(spectrum_.begin());

  for (const auto& rp: exp_spectrum_) {
    double thr = rp - PRODUCT_ION_THRESH - MASS_PROTON;

    while (pmb != spectrum_.end() && *pmb <= thr) {
      pmb++;
    }
    
    if (pmb !=spectrum_.end()) {
      if (std::abs(*pmb + MASS_PROTON - rp) < PRODUCT_ION_THRESH) {
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

    while (pme != rpme && NLP_MASS - *pme <= thr) {
      pme++;
    }

    if (pme != rpme) {
      if (std::abs(NLP_MASS - *pme + MASS_PROTON - rp) < PRODUCT_ION_THRESH) {
        score_ += 1;
        // std::cout << nlp_mass - *pme + MASS_PROTON << " " << rp << std::endl;
        continue;
      }
    }

    if (pmb == spectrum_.end() && pme == rpme) {
      break;
    }
    // score_ = std::min(score_, MAX_SCORE);
  }
}
