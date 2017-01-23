#include "utils.h"

#include <cfloat>
#include <cassert>
#include "scorer.h"

const double MASS_PROTON = 1.00728;

using namespace scoring;
template <class It1, class It2>
static double score_match_absolute(It1 pe, It1 e_end,
                                   It2 pt, It2 t_end,
                                   double delta,
                                   double threshold = 0.) {
    unsigned score = 0;

    // assert(std::is_sorted(pe, e_end));
    // assert(std::is_sorted(pt, t_end));

    while (pe != e_end && pt != t_end) {
        // std::cout << *pe - MASS_PROTON - *pt << "\n";
        if (std::abs(*pe - MASS_PROTON - *pt) < delta) {
            // std::cout << *pe << " === " << *pt << "\n";
            score += 1;

            auto prev_e_match = *pe;
            auto prev_t_match = *pt;
            ++pt; ++pe;

            if (threshold) {
                while (pe != e_end && *pe - prev_e_match < threshold) ++pe;  // < for consistency with threshold = 0 case
                while (pt != t_end && *pt - prev_t_match < threshold) ++pt;
            }
        } else if (*pe - MASS_PROTON < *pt) {
            ++pe;
        } else {
            ++pt;
        }
    }

    return score;
}

template <class It1, class It2>
static double score_match_relative(It1 pe, It1 e_end,
                                   It2 pt, It2 t_end,
                                   double delta,
                                   double threshold = 0.) {
    unsigned score = 0;

    // assert(std::is_sorted(pe, e_end));
    // assert(std::is_sorted(pt, t_end));

    // Rewind negative peaks
    while (pe != e_end && *pe - MASS_PROTON <= 0) ++pe;
    while (pt != t_end && *pt <= 0) ++pt;

    while (pe != e_end && pt != t_end) {
        if (*pe - MASS_PROTON < (1 + delta) * *pt && *pt < (1 + delta) * (*pe - MASS_PROTON)) {
            score += 1;
            //std::cout << *pe << " === " << *pt << "\n";
            auto prev_e_match = *pe;
            auto prev_t_match = *pt;
            ++pt; ++pe;

            if (threshold) {
                while (pe != e_end && (*pe - MASS_PROTON) < (prev_e_match - MASS_PROTON) * (1 + threshold)) ++pe;  // < for consistency with threshold = 0 case
                while (pt != t_end && *pt < prev_t_match * (1 + threshold)) ++pt;
            }
        } else if (*pe - MASS_PROTON < *pt) {
            ++pe;
        } else {
            ++pt;
        }
    }

    return score;
}

template<class It>
class spectrum_adaptor : public llvm::iterator_adaptor_base<spectrum_adaptor<It>,
                                                            It,
                                                            typename std::iterator_traits<It>::iterator_category,
                                                            decltype(std::declval<It>()->mass)> {
    using mass_type = decltype(std::declval<It>()->mass);
  public:
    // This is really hack, but we simply provide some 'view', so should be relatively safe.
    typedef mass_type reference;

    spectrum_adaptor(It it)
            : spectrum_adaptor::iterator_adaptor_base(it) {}

    mass_type operator*() const { return this->I->mass; }
};

template<class It>
spectrum_adaptor<It> make_spectrum(It it) {
    return spectrum_adaptor<It>(it);
}

template <class T>
class assign_iterator : public std::iterator<std::output_iterator_tag,
                                             void, void, void, void> {
  protected:
    T& container_;
    size_t idx_;

  public:
    typedef T container_type;

    explicit assign_iterator(T& x, size_t sz)
            : container_(x), idx_(0) {
        container_.clear();
        container_.resize(sz);
    }
    assign_iterator& operator=(const typename T::value_type& value) {
        container_[idx_] = value;
        return *this;
    }
    assign_iterator& operator=(typename T::value_type&& value) {
        container_[idx_] = std::move(value);
        return *this;
    }
    typename T::value_type& operator*() { return container_[idx_]; }
    assign_iterator& operator++()    { idx_ += 1; return *this; }
    assign_iterator  operator++(int) { assign_iterator tmp = *this; idx_ +=1; return tmp; }
};

template<class T>
assign_iterator<T> make_assigner(T &container, size_t sz) {
    return assign_iterator<T>(container, sz);
}

double SPCScorer::score(Spectrum &spectrum, Peptide &nlp, bool keep_both) const {
    unsigned max_charge = spectrum.charge_;

    // bool keep_both = nlp.keep_both;
    double res = 0;
   
    // Expand spectrum, if necessary
    tspectrum_.clear();
    tspectrumkb_.clear();

    std::vector<double> & otspectrum = nlp.sorted_spectrum_;

    
    // for (auto & entry :  otspectrum) {
    //         printf("%e\n", entry) ;
    // }

    tspectrumkb_.reserve(otspectrum.size() * (keep_both ? 1 : 2));
    tspectrum_.reserve(max_charge * tspectrumkb_.capacity());


    if (!keep_both) {
        add_right_tail_to_spectrum_sorted(otspectrum.begin(),
                                          otspectrum.end(),
                                          make_assigner(tspectrumkb_, otspectrum.size() * 2),
                                          nlp.get_peptide_mass_());
    } 
    // std::cout << "after add_right_tail" << std::endl;
    
    // for (auto & entry :  tspectrum_) {
    //         printf("%e\n", entry) ;
    // }
    // std::cout << std::endl;

    add_charged_peaks_to_spectrum_sorted(tspectrumkb_.begin(), tspectrumkb_.end(),
                                         make_assigner(tspectrum_, max_charge * tspectrumkb_.size()),
                                         max_charge);
    
   // std::cout << "after add_peaks" << std::endl;
   // for (auto & entry :  tspectrum_) {
   //         printf("%e\n", entry) ;
   // }
   // std::cout << std::endl;

    const auto& raw_peaks = spectrum.exp_spectrum_;

    // std::cout << "exp_spectrum_" << std::endl;
    // for (auto & entry : raw_peaks) {
    //         printf("%e\n", entry) ;
    // }
    // std::cout << std::endl;
    
    // std::sort(tspectrum_.begin(), tspectrum_.end());
    
    if (ppm_threshold_)
        res = score_match_relative(raw_peaks.begin(), raw_peaks.end(),
                                   tspectrum_.cbegin(), tspectrum_.cend(),
                                   product_ion_thresh_);
    else
        res = score_match_absolute(raw_peaks.begin(), raw_peaks.end(),
                                   tspectrum_.cbegin(), tspectrum_.cend(),
                                   product_ion_thresh_);

    return res;

}

double SPCScorer::score_peak(Spectrum & spectrum, Peptide & peptide , bool keep_both) {
  double nlp = peptide.peptide_mass_;

  std::vector<double> spectrum_ = peptide.sorted_spectrum_;
  // std::vector<double> spectrum_(sort_perm.size());

  // for (int i = 0; i < sort_perm.size(); ++i) {
  //   spectrum_[i]  = spect[sort_perm[i]];
  // }

  // for (int i = 0; i < spectrum_.size(); ++i)
  // {
  //   std::cout << spectrum_[i] << " : " << peptide.sorted_spectrum_[i];
  //   // std::cout << spectrum_[i] << " " ;
  // }


  // for (int i = 0; i < peptide.sorted_spectrum_.size(); ++i)
  // {
  //   std::cout << spectrum_[i] << " : " << peptide.sorted_spectrum_[i] << std::endl;
  //   // std::cout << peptide.sorted_spectrum_[i] << " " ;
  // }

  
  
  double score_ = 0;

  auto pmb = spectrum_.begin();
  auto pme = std::reverse_iterator<decltype(spectrum_.end())>(spectrum_.end());
  auto rpme = std::reverse_iterator<decltype(spectrum_.begin())>(spectrum_.begin());

  for (const auto& rp: spectrum.exp_spectrum_) {
    double thr = rp - product_ion_thresh_ - MASS_PROTON;

    while (pmb != spectrum_.end() && *pmb <= thr) {
      pmb++;
    }
    
    if (pmb !=spectrum_.end()) {
      if (std::abs(*pmb + MASS_PROTON - rp) <= product_ion_thresh_) {
        score_ += 1;
        // printf("bion %e : %e\n", *pmb + MASS_PROTON, rp);
        // printf("real = %e : %e\n",  *pmb , rp);
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
      if (std::abs(nlp - *pme + MASS_PROTON - rp) <= product_ion_thresh_) {
        score_ += 1;
        // printf("yion %e : %e\n", nlp - *pme + MASS_PROTON, rp);
        // printf("real = %e : %e\n",  *pme ,rp);
        continue;
      }
    }

    if (pmb == spectrum_.end() && pme == rpme) {
      break;
    }
    
  }
  return score_;
  // std::cout << "+++++++++" << std::endl;
}

