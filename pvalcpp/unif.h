#pragma once 

const double NLP_MASS = 1007.65;
const double MASS_PROTON = 1.00728;
const double PRODUCT_ION_THRESH = 0.5;
const double MIN_SCORE = 0;
const double MAX_SCORE = 14;
const double PHI_B = exp(0.6);
const double PHI_E = exp(0.0003);
const int STEP_LENGTH = 10000;

static double unif_rand() {
    return double(rand()) / RAND_MAX;
}


static double unif(double a, double b) {
    return unif_rand() * (b - a) + a;
}
