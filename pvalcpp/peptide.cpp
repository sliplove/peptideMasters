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

std::vector<double>  update_spectrum_by_new_mass(const std::vector<double> &  mass_, 
    const std::vector<std::pair<unsigned, unsigned>> & rule_, Peptide & peptide) {
    
    unsigned id = floor(unif_rand() * rule_.size());
    unsigned beg = rule_[id].first;
    unsigned end = rule_[id].second;
    
    double beg_mass = mass_[beg], end_mass = mass_[end];
    double delta = unif(-beg_mass, end_mass);

    std::vector<double> mass (mass_);
    
    mass[beg] = beg_mass + delta;
    mass[end] = end_mass - delta;
    
    auto it_moved_plus = peptide.moved_plus_.begin();
    auto it_moved_minus = peptide.moved_minus_.begin();
    auto it_unmoved = peptide.unmoved_.begin();

    for (size_t i : peptide.sorting_permutation_) {
        if (peptide.mat_[i][beg] == peptide.mat_[i][end]) {
            *it_unmoved++ = i;
        } else if (peptide.mat_[i][beg]) {
            *it_moved_plus++ = i;
            peptide.spectrum_[i] += delta;
        } else {
            *it_moved_minus++ = i;
            peptide.spectrum_[i] -= delta;
        }
    }

    peptide.unmoved_.resize(peptide.spectrum_.size());
    peptide.moved_plus_.resize(peptide.spectrum_.size());
    peptide.moved_minus_.resize(peptide.spectrum_.size());

    // Add sentineles
    *it_moved_plus = peptide.spectrum_.size() - 1;
    *it_moved_minus = peptide.spectrum_.size() - 1;
    *it_unmoved = peptide.spectrum_.size() - 1;

    // Rewind
    it_moved_plus = peptide.moved_plus_.begin();
    it_moved_minus = peptide.moved_minus_.begin();
    it_unmoved = peptide.unmoved_.begin();
    
    // 3-merge
    for (size_t i = 0; i < peptide.sorting_permutation_.size(); ++i) {
        if (peptide.spectrum_[*it_unmoved] < peptide.spectrum_[*it_moved_plus] && 
            peptide.spectrum_[*it_unmoved] < peptide.spectrum_[*it_moved_minus]) {
            peptide.sorting_permutation_[i] = *it_unmoved++;
    } else if (peptide.spectrum_[*it_moved_minus] < peptide.spectrum_[*it_moved_plus]) {
        peptide.sorting_permutation_[i] = *it_moved_minus++;
    } else {
        peptide.sorting_permutation_[i] = *it_moved_plus++;
    }
}

    // std::sort(peptide.spectrum_.begin(), peptide.spectrum_.end());

return mass;
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



void Peptide::copy_spectrum_(const Peptide & peptide) {
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
