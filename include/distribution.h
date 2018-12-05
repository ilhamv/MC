#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H


struct DistributionWatt
{
	// Vector of parameters corresponding to thermal(< 1eV), 1 MeV, and 14 MeV
	std::vector<double> vec_a;    
	std::vector<double> vec_b;
	std::vector<double> vec_g;
	
	DistributionWatt(const std::vector<double> p1, const std::vector<double> p2,
                   const std::string label = "");
	double sample(const double E = 0.0);
};


#endif // DISTRIBUTION_H
