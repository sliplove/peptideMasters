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
#include "main_functions.h"
#include "unif.h"



int main(int argc, char const *argv[])
{
	srand (13);

	// read input files
	std::vector<std::vector<double> > mat;
    std::ifstream file_mat("../input/matrix.txt");
	std::ifstream file_rule("../input/rule_graph.txt");
	std::ifstream file_spectrum("../input/spectrum.txt");
	

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

	std::cout << ncol  << " " << nrow << std::endl;

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
	
	// get weights	
	WLsimulator wl(MIN_SCORE, MAX_SCORE, PHI_B, PHI_E, STEP_LENGTH);

	std::vector<double> weights;
	weights = wl.wl_full(exp_spectrum, mat, rule, NLP_MASS, true);

	std::cout << "Wang-Landau weights" << std::endl;
	for (auto & w: weights) {
		std::cout << w << " " ;
	}

	std::cout << std::endl << std::endl;

	// mh step
	std::vector<double> start_mass = get_start_mass(ncol, NLP_MASS);

    hit_run(exp_spectrum, mat, rule, start_mass, NLP_MASS,  wl, 50000 , 500000, 0.02, 0.95);

	return 0;
}