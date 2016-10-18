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

std::vector<double> update_mass_single_(const std::vector<double> &  mass_, 
            const std::vector<std::pair<unsigned, unsigned>> & rule_) {
    
    unsigned id = floor(unif_rand() * rule_.size());
    unsigned beg = rule_[id].first;
    unsigned end = rule_[id].second;
    
    double beg_mass = mass_[beg], end_mass = mass_[end];
    double delta = unif(-beg_mass, end_mass);

    std::vector<double> mass (mass_);
    
    mass[beg] = beg_mass + delta;
    mass[end] = end_mass - delta;

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


    std::sort(spectrum_.begin(), spectrum_.end());
}


void Peptide::set_spectrum(const std::vector<double> & mass) {
    
    int nrow = mat_.size() - 1;
    int ncol = mat_[0].size();

     for (size_t i = 0; i < nrow; ++i) {
        double res = 0;
        for (size_t j = 0; j < ncol; ++j)
            res += mass[j] * mat_[i][j];
        spectrum_[i] = res;
    }


    std::sort(spectrum_.begin(), spectrum_.end());
}



   
void Peptide::print() {

    std::cout << "Spectrum of peptide: " << std::endl;

    for (auto & entry : spectrum_) 
        std::cout << entry << " "; 
        std::cout <<  std::endl;


}

