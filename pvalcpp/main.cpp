#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "cxxopts.hpp"
#include "peptide.h"
#include "wl.h"
#include "metropolis.h"
#include "mhstate.h"
#include "scorer.h"

int main(int argc, char *argv[])
{
	srand (13);
	cxxopts::Options options(argv[0], "Markov Chain Monte Carlo method");
	options.add_options()
	("h,help", "Print help")
	("s,spectrum", "Input spectrum", cxxopts::value<std::string>(), "FILE")
	("m,matrix", "Fragmentation Marix", cxxopts::value<std::string>(), "FILE")
	("r,rule", "Rule graph", cxxopts::value<std::string>(), "FILE")
	("precursor", "Precursor mass", cxxopts::value<double>(), "FLOAT")
	("min", "Min score", cxxopts::value<double>()->default_value("0"), "FLOAT")
	("max", "Max score", cxxopts::value<double>(), "FLOAT")
	("phi_begin", "Initial phi value (WL option)", cxxopts::value<double>()->default_value("1.8"), "FLOAT")
	("phi_end", "Final phi value (WL option)", cxxopts::value<double>()->default_value("1"), "FLOAT")
	("step", "Number of Wang-Landau iterations", cxxopts::value<unsigned>()->default_value("10000"), "N")  
	("run_iter", "Number of Monte-Carlo iterations", cxxopts::value<unsigned>()->default_value("50000"), "N")  
	("eps", "Accuracy", cxxopts::value<double>()->default_value("0.02"), "FLOAT")
	("level", "Quantile level for confident interval", cxxopts::value<double>()->default_value("0.95"), "FLOAT")
	("product_ion_thresh", "Score parameter", cxxopts::value<double>()->default_value("0.5"), "FLOAT");
     
    // options.parse_positional(std::vector<std::string>({"spectrum", "matrix", "rule"}));
	options.parse(argc, argv);

	if (options.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

	std::ifstream file_mat(options["matrix"].as<std::string>());
	std::ifstream file_rule(options["rule"].as<std::string>());
	std::ifstream file_spectrum(options["spectrum"].as<std::string>());
	double NLP_MASS = options["precursor"].as<double>();
	double MIN_SCORE = options["min"].as<double>();
	double MAX_SCORE = options["max"].as<double>();
	double PHI_B = options["phi_begin"].as<double>();
	double PHI_E = options["phi_end"].as<double>();	
	unsigned STEP_LENGTH = options["step"].as<unsigned>();
	unsigned MIN_STEPS_RUN = options["run_iter"].as<unsigned>();
	double EPS = options["eps"].as<double>();	
	double LEVEL = options["level"].as<double>();
	double PRODUCT_ION_THRESH = options["product_ion_thresh"].as<double>();


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

	std::cout << MIN_SCORE << " " << MAX_SCORE << " " << PHI_B << " " <<
									 PHI_E << " " << STEP_LENGTH << std::endl;

	// set scorer and metropolis parameters 
	Scorer scorer(exp_spectrum, PRODUCT_ION_THRESH);
	MHstate state(ncol, NLP_MASS);
	Peptide peptide(mat, rule, state.get_current_state_(), NLP_MASS);

	Metropolis mh(mat, rule, NLP_MASS, MIN_SCORE, MAX_SCORE, state, peptide, scorer);
	
	// get weights	
	WLsimulator wl(mh, PHI_B, PHI_E, STEP_LENGTH);
	wl.print();

	// wl.wl_step(PHI_B, true);
	std::vector<double> weights;
	weights = wl.wl_full(true);

	std::cout << "Wang-Landau weights" << std::endl;
	for (auto & w: weights) {
		std::cout << w << " " ;
	}

	std::cout << std::endl << std::endl;

	// mh step
	mh.hit_run(MIN_STEPS_RUN, EPS, LEVEL, weights);

	return 0;
}