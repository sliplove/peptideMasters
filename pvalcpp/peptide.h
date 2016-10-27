#pragma once 
#include <cmath>
#include <cstdint>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <iostream>
#include <cassert>
#include <queue>
#include "unif.h"

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
    std::vector<size_t> get_sorting_permutation_() const {
        return sorting_permutation_; 
    }

    void print();
    friend std::vector<double> update_spectrum_by_new_mass(const std::vector<double> & ,
        const std::vector<std::pair<unsigned, unsigned>> &,
        Peptide & );
    void clear(const std::vector<double> &);

};
