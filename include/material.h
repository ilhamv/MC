#ifndef _MATERIAL_H
#define _MATERIAL_H

#include <vector>
#include <memory>
#include <cstring>

#include "type.h"

//==============================================================================
// Nuclide
//==============================================================================

struct Nuclide
{
  const std::string name;
	const double      A;
	const bool        fissionable;

	// Watt fission spectrum parameters corresponding to 
	// thermal(< 1eV), 1 MeV, and 14 MeV
  std::vector<double> a;
  std::vector<double> b;
  std::vector<double> g;

  // micro XS data in Energy
	// 0: energy (eV)
	// 1: total
	// 2: capture
	// 3: scattering
	// 4: fission
	// 5: nu
	const std::vector<std::vector<double>> XS;
	
	Nuclide(const std::string n, const double AA, const bool f, 
			    const std::vector<double> wa, const std::vector<double> wb, 
			    const std::vector<double> wg,
					const std::vector<std::vector<double>> xs): 
		name(n), A(AA), a(wa), b(wb), g(wg), fissionable(f), XS(xs) {};
};

//==============================================================================
// Material
//==============================================================================

struct Material
{
	const std::string name;
	const ull         ID;

	// Pair of nuclide and its nuclide density
	const std::vector<std::pair<std::shared_ptr<Nuclide>,double>> nuclides;

	Material(const std::string n, const ull i,
			     const std::vector<std::pair<std::shared_ptr<Nuclide>,double>> nd):
		name(n), ID(i), nuclides(nd) {};
};


#endif // MATERIAL_H
