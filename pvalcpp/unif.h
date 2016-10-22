#pragma once 

const double MASS_PROTON = 1.00728;
const double PRODUCT_ION_THRESH = 0.5;

static double unif_rand() {
	return double(rand()) / RAND_MAX;
}


static double unif(double a, double b) {
	return unif_rand() * (b - a) + a;
}
