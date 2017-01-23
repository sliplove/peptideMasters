#pragma once 
class Spectrum {
public:
	std::vector<double> exp_spectrum_;
	int charge_; 

	Spectrum();
	Spectrum(const std::vector<double> &exp_spectrum, int charge)
	:
	exp_spectrum_(exp_spectrum), 
	charge_(charge)
	{}
	void print() {
		for (int i = 0; i < exp_spectrum_.size(); ++i)
		{
			std::cout << exp_spectrum_[i] << " ";
		}
		std::cout << std::endl;
	}
};
