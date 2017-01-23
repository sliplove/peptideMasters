#pragma once 
#ifndef _SCORING_H
#define _SCORING_H

#include <cstring>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>

#include "peptide.h"
#include "spectrum.h"

namespace scoring {

class SPCScorer {
protected:
    bool ppm_threshold_;
    double product_ion_thresh_;
private:
    // Some temporary buffers in order to save memory allocations
    mutable std::vector<double> tspectrum_, tspectrumkb_;

public:
    SPCScorer(double product_ion_thresh, bool ppm_threshold = false)
    :
    product_ion_thresh_(product_ion_thresh),
    ppm_threshold_(ppm_threshold) {}
    SPCScorer(const SPCScorer & scorer)
    :
    product_ion_thresh_(scorer.product_ion_thresh_),
    ppm_threshold_(scorer.ppm_threshold_) {}

    double score(Spectrum &spectrum, Peptide &nlp, bool) const;
    double score_peak(Spectrum & , Peptide &,  bool);
};
}

#endif
