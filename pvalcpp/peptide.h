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

std::vector<double> update_mass_single_(const std::vector<double> & ,
            const std::vector<std::pair<unsigned, unsigned>> &);

class Peptide {
private:
    std::vector<std::pair<unsigned, unsigned>> rule_;
    std::vector<std::vector<double>> mat_;
    double peptide_mass_;
    std::vector<double> spectrum_;
    // double score_;

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
                peptide_mass_(peptide.peptide_mass_)
                {}
	
    // double get_score() { return score_; }
    void set_spectrum(const std::vector<double> & );
    std::vector<double> get_spectrum_() {return spectrum_; }
	void print();

};
