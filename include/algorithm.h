#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <vector>

#include "type.h"
#include "particle.h"

ull binary_search(const double x, const std::vector<double>& vec);
double calculate_XS(const ull idx, const double E, 
		                const std::vector<std::vector<double>>& XS, 
										const int header);
double interpolate(const double x, const double x1, const double x2,
                   const double y1, const double y2 );
Point scatter_direction(const Point dir_i, const double mu0);
double watt_spectrum(const double E, const std::vector<double>& vec_a, 
		                 const std::vector<double>& vec_b, 
										 const std::vector<double>& vec_g);

#endif // ALGORITHM_H
