#pragma once 
#include <vector>


std::vector<double> get_start_mass(unsigned , double);

class Peptide {
private:
    std::vector<std::pair<unsigned, unsigned>> rule_;
    std::vector<std::vector<double>> mat_;
    double peptide_mass_;
    std::vector<double> spectrum_;
    std::vector<size_t> sorting_permutation_;
    std::vector<size_t> unmoved_, moved_plus_, moved_minus_;

public:
    Peptide() {};
    Peptide(const std::vector<std::vector<double>> & ,
       const std::vector<std::pair<unsigned, unsigned>> & ,
       const std::vector<double> &, double);
    Peptide (const Peptide & peptide)
    :
    rule_(peptide.rule_),
    mat_(peptide.mat_),
    spectrum_(peptide.spectrum_),
    sorting_permutation_(peptide.sorting_permutation_),
    unmoved_(peptide.unmoved_),
    moved_minus_(peptide.moved_minus_),
    peptide_mass_(peptide.peptide_mass_)
    {}
    
    void set_spectrum_(const Peptide & );
    
    double get_peptide_mass_() const { return peptide_mass_; }
    std::vector<double> get_spectrum_() const { return spectrum_; }
    
    double get_size() const { return spectrum_.size(); }
   
    std::vector<size_t> get_sorting_permutation_() const {
        return sorting_permutation_; 
    }

    void increment_spectrum_(unsigned i, double delta) {
        spectrum_[i] += delta;
    }

    void set_sorting_permutation_(unsigned i, unsigned value) {
        sorting_permutation_[i] = value;
    }
 
    void print();
    void clear(const std::vector<double> &);

    auto moved_minus_begin() {return moved_minus_.begin();}
    auto moved_plus_begin() {return moved_plus_.begin();}
    auto unmoved_begin() {return unmoved_.begin();}
    

};
