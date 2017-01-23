#pragma once 
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include "../include/pcg_random.hpp"

static double unif_real(pcg_extras::seed_seq_from<std::random_device> & seed_source, double a, double b) {
    pcg32 rng(seed_source);
    std::uniform_real_distribution<> dis(a, b);
    double rval = dis(rng);   
    return rval;
}


static int unif_int(pcg_extras::seed_seq_from<std::random_device> & seed_source, int a, int b) {
    pcg32 rng(seed_source);
    std::uniform_int_distribution<> dis(a, b);
    double rval = dis(rng);   
    return rval;
}
