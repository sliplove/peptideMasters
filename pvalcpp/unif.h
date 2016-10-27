#pragma once 

static double unif_rand() {
	return double(rand()) / RAND_MAX;
}


static double unif(double a, double b) {
	return unif_rand() * (b - a) + a;
}
