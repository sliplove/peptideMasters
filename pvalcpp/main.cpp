#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "peptide.h"
#include "wl.h"
#include "metropolis.h"
#include "unif.h"


int main(int argc, char *argv[])
{
	srand (13);

	// read input files

	std::ifstream file_mat(argv[1]);
	std::ifstream file_rule(argv[2]);
	std::ifstream file_spectrum(argv[3]);

	double NLP_MASS = atof(argv[4]);
	double MIN_SCORE = atoi(argv[5]);
	double MAX_SCORE = atoi(argv[6]);
	double PHI_B = atof(argv[7]);
	double PHI_E = atof(argv[8]);	
	int STEP_LENGTH = atoi(argv[9]);
	int MIN_STEPS_RUN = atoi(argv[10]);
	double EPS = atof(argv[11]);
	double LEVEL = atof(argv[12]);



	std::vector<std::vector<double> > mat;

    double elem;
	std::string line;
	
	while(!file_mat.eof()) {
		getline(file_mat, line);
		std::istringstream iss(line);
		std::vector<double> row;
		while (iss >> elem) {
			row.push_back(elem);
		}
		mat.push_back(row);
	}
	
	int nrow = mat.size() - 1;
	int ncol = mat[0].size();

	// std::cout << ncol  << " " << nrow << std::endl;

    std::vector<std::pair<unsigned, unsigned> > rule (ncol);
	double elem1, elem2;
	int i = 0;
	while(!file_rule.eof()) {
		getline(file_rule, line);
		std::istringstream iss(line);
		while (iss >> elem1 >> elem2) {
			rule[i].first = elem1;
			rule[i].second = elem2;
		}
		++i;
	}

	std::vector<double> exp_spectrum;

	while(!file_spectrum.eof()) {
		getline(file_spectrum, line);
		std::istringstream iss(line);
		while (iss >> elem) {
			exp_spectrum.push_back(elem);
		}
	}

	file_mat.close();
	file_rule.close();
	file_spectrum.close();
	std::cout << MIN_SCORE << MAX_SCORE << PHI_B << PHI_E << STEP_LENGTH;

	// get weights	
	WLsimulator wl(MIN_SCORE, MAX_SCORE, PHI_B, PHI_E, STEP_LENGTH);
	wl.print();

	std::vector<double> weights;
	weights = wl.wl_full(exp_spectrum, mat, rule, NLP_MASS, true);

	std::cout << "Wang-Landau weights" << std::endl;
	for (auto & w: weights) {
		std::cout << w << " " ;
	}

	std::cout << std::endl << std::endl;

	// mh step
	Metropolis run(exp_spectrum, mat, rule, NLP_MASS, wl);

	std::vector<double> start_mass = get_start_mass(ncol, NLP_MASS);
    run.hit_run(start_mass, MIN_STEPS_RUN , EPS, LEVEL);

	return 0;
}