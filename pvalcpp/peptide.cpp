#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <iostream>

#include "unif.h"
#include "peptide.h"

std::vector<double> get_start_mass(unsigned N, double peptide_mass_) {
    std::vector<double> mass_(N);

    double sum = 0;
    for (size_t i = 0; i < N; ++i) {
        double val = -log(unif_rand());
        sum += val;
        mass_[i] = val;
    }

    double coef = peptide_mass_ / sum;
    for (auto &entry : mass_)
        entry *= coef;

    return mass_;
}


Peptide::Peptide(const std::vector<std::vector<double>> &mat,
  const std::vector<std::pair<unsigned, unsigned>> &rule,
  const std::vector<double> & mass, 
  double peptide_mass) : rule_(rule), mat_(mat), peptide_mass_(peptide_mass) {
    
    int nrow = mat_.size() - 1;
    int ncol = mat_[0].size();


    spectrum_.resize(nrow);

    for (size_t i = 0; i < nrow; ++i) {
        double res = 0;
        for (size_t j = 0; j < ncol; ++j)
            res += mass[j] * mat_[i][j];
        spectrum_[i] = res;
    }

    // Add sentinel
    spectrum_.push_back(std::numeric_limits<double>::max());

    sorting_permutation_.resize(spectrum_.size() - 1);
    std::iota(sorting_permutation_.begin(),
      sorting_permutation_.end(),
      0);
    std::sort(sorting_permutation_.begin(),
      sorting_permutation_.end(),
      [this](size_t i, size_t j) -> bool { return this->spectrum_[i] < this->spectrum_[j]; });

    unmoved_.resize(spectrum_.size());
    moved_plus_.resize(spectrum_.size());
    moved_minus_.resize(spectrum_.size());

    // std::sort(spectrum_.begin(), spectrum_.end());
}



void Peptide::set_spectrum_(const Peptide & peptide) {
    spectrum_ = peptide.get_spectrum_();
    sorting_permutation_ = peptide.get_sorting_permutation_();  
}

void Peptide::print() {

    std::cout << "Spectrum of peptide: " << std::endl;

    for (auto & entry : spectrum_) 
        std::cout << entry << ", "; 
    std::cout <<  std::endl;

    std::cout << "Sorted Spectrum of peptide: " << std::endl;

    for (auto & entry : sorting_permutation_) 
        std::cout << spectrum_[entry] << ", "; 
    std::cout <<  std::endl;


    std::cout << "Sorting permutation: " << std::endl;

    for (auto & entry : sorting_permutation_) 
        std::cout << entry << " "; 
    std::cout <<  std::endl;
    

}

void Peptide::clear(const std::vector<double> & mass) {
    int nrow = mat_.size() - 1;
    int ncol = mat_[0].size();

    spectrum_.resize(nrow);

    for (size_t i = 0; i < nrow; ++i) {
        double res = 0;
        for (size_t j = 0; j < ncol; ++j)
            res += mass[j] * mat_[i][j];
        spectrum_[i] = res;
    }

    // Add sentinel
    spectrum_.push_back(std::numeric_limits<double>::max());

    sorting_permutation_.resize(spectrum_.size() - 1);
    std::iota(sorting_permutation_.begin(),
      sorting_permutation_.end(),
      0);
    std::sort(sorting_permutation_.begin(),
      sorting_permutation_.end(),
      [this](size_t i, size_t j) -> bool { return this->spectrum_[i] < this->spectrum_[j]; });

    unmoved_.resize(spectrum_.size());
    moved_plus_.resize(spectrum_.size());
    moved_minus_.resize(spectrum_.size());

}
