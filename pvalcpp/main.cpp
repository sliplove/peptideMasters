#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>

#include "cxxopts.hpp"
#include "peptide.h"
#include "wl.h"
#include "metropolis.h"
#include "mhstate.h"
#include "scorer.h"
#include "unif.h"
#include "../include/pcg_random.hpp"


double st_m[] = { 29.498327, 42.155009, 7.926439, 7.030549, 118.557788,  
	12.898671, 134.758404,  66.869701,  41.517182 }; 	

int main(int argc, char *argv[])
{
	srand (42);
	cxxopts::Options options(argv[0], "Markov Chain Monte Carlo method");
	options.add_options("General")
		("h,help", "Print help")
		("s,spectrum", "Input spectrum", cxxopts::value<std::string>(), "FILE")
		("m,matrix", "Fragmentation Marix", cxxopts::value<std::string>(), "FILE")
		("r,rule", "Rule graph", cxxopts::value<std::string>(), "FILE")
		("precursor", "Precursor mass", cxxopts::value<double>(), "FLOAT")
		("min", "Min score", cxxopts::value<double>()->default_value("0"), "FLOAT")
		("max", "Max score", cxxopts::value<double>(), "FLOAT")
		("charge", "Spectrum parameter", cxxopts::value<unsigned>()->default_value("1"), "FLOAT");
	
	options.add_options("Advanced")
		("phi_begin", "Initial phi value (WL option)", cxxopts::value<double>()->default_value("1.822"), "FLOAT")
		("phi_end", "Final phi value (WL option)", cxxopts::value<double>()->default_value("1"), "FLOAT")
		("step", "Length of Wang-Landau iteration", cxxopts::value<unsigned>()->default_value("10000"), "N")  
		("run_iter", "Number of Monte-Carlo iterations", cxxopts::value<unsigned>()->default_value("50000"), "N")  
		("eps", "Accuracy", cxxopts::value<double>()->default_value("0.02"), "FLOAT")
		("level", "Quantile level for confident interval", cxxopts::value<double>()->default_value("0.95"), "FLOAT")
		("product_ion_thresh", "Score parameter", cxxopts::value<double>()->default_value("0.5"), "FLOAT");
	    

    const std::vector<std::string> all_groups({"General", "Advanced"});
 
    // options.parse_positional(std::vector<std::string>({"spectrum", "matrix", "rule"}));
	options.parse(argc, argv);
    
	if (options.count("help")) {
        std::cout << options.help(all_groups) << std::endl;
        exit(0);
    }

	std::ifstream file_mat(options["matrix"].as<std::string>());
	std::ifstream file_rule(options["rule"].as<std::string>());
	std::ifstream file_spectrum(options["spectrum"].as<std::string>());
	double NLP_MASS = options["precursor"].as<double>();
	double MIN_SCORE = options["min"].as<double>();
	double MAX_SCORE = options["max"].as<double>();
	unsigned CHARGE = options["charge"].as<unsigned>();
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

	std::cout << nrow << " " << ncol << std::endl;
	
	std::vector<std::pair<unsigned, unsigned> > rule;
	double elem1, elem2;
	int i = 0;
	std::istringstream iss(line);

	while(file_rule >> elem1 >> elem2) {
		rule.push_back(std::make_pair(elem1, elem2));
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
	
	pcg_extras::seed_seq_from<std::random_device> rd;   
	// for (int i = 0; i < nrow; ++i) {
	// for (int j = 0; j < ncol; ++j) {
	// 	std::cout << mat[i][j] << " ";
	// 	}
	// 	std::cout << std::endl;
	// }

	// std::cout << " ----------------- " << std::endl;
	// for (int i = 0; i < rule.size(); ++i) {
	// 	std::cout << rule[i].first << " " << rule[i].second << std::endl;
	// }

	// std::cout << " ----------------- " << std::endl;


	// std::cout << MIN_SCORE << " " << MAX_SCORE << " " << PHI_B << " " <<
	// 								 PHI_E << " " << STEP_LENGTH << " " << NLP_MASS << std::endl;

	
	// set scorer and metropolis parameters 
	Spectrum spectrum(exp_spectrum, CHARGE);
	SPCScorer scorer(PRODUCT_ION_THRESH);
	// std::vector<double> start_mass(st_m, st_m + sizeof(st_m) / sizeof(st_m[0]));
	// MHstate state(start_mass);

	MHstate state(ncol, NLP_MASS, rd);

	Peptide peptide(mat, rule, state.get_current_state_(), NLP_MASS);
	Metropolis mh(mat, rule, NLP_MASS, MIN_SCORE, MAX_SCORE, state, peptide, spectrum, scorer);

	// get weights	
	WLsimulator wl(mh, PHI_B, PHI_E, STEP_LENGTH);
	// wl.print();

	// wl.wl_step(rd, PHI_B, true);
	 
	std::vector<double> weights;
	weights = wl.wl_full(rd, false);

	// std::cout << "Wang-Landau weights" << std::endl;
	// for (auto & w: weights) {
	// 	std::cout << w << " " ;
	// }

	// std::cout << std::endl << std::endl;

	// // mh step
	mh.hit_run(rd, MIN_STEPS_RUN, EPS, LEVEL, weights);

	return 0;
}