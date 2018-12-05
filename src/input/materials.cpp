#include <fstream>

#include "input.h"
#include "constant.h"

void setup_materials(const pugi::xml_node input_materials, 
		                 std::vector<std::shared_ptr<Nuclide>>& nuclides,
  							     std::vector<std::shared_ptr<Material>>& materials,
										 const std::string xs_directory)
{
	//===========================================================================
	// Nuclides: Sweep through material and define the nuclides
	//===========================================================================
	
	for(const auto& node_material : input_materials.children("material")){
		for(const auto& node_nuclide : node_material.children("nuclide")){
			std::shared_ptr<Nuclide> the_nuclide;

			// Inputs
			const std::string name = node_nuclide.attribute("name").value();
			const double      ao   = node_nuclide.attribute("ao").as_double();

			// Skip if nuclide is already defined
			if(find_by_name(nuclides,name)) continue;

			// XS data
      const std::string xs_file = xs_directory + name + ".txt";
      std::ifstream     xs_data(xs_file);
			std::string       line;
			// A
			std::getline(xs_data, line);
			const double A = parse_vector<double>(line)[0];

			// Watt spectrum
      std::vector<double> a(3);
      std::vector<double> b(3);
      std::vector<double> g(3);
			for(int i = 0; i < 3; i++){
				std::getline(xs_data, line);
				std::vector<double> tmp = parse_vector<double>(line);
				a[i] = tmp[0];
				b[i] = tmp[1];
        const double C = (1.0 + a[i]*b[i]/8.0);
				g[i] = std::sqrt(C*C - 1.0) + C;
			}
			// Fissionable
			bool fissionable = false;
			for(int i = 0; i < 3; i++)
				if(a[i] != 0.0 || b[i] != 0.0) fissionable = true;

			// XS
			// 0: energy (eV)
			// 1: total
			// 2: capture
			// 3: scattering
			// 4: fission
			// 5: nu-fission
			std::vector<std::vector<double>> XS(6);
			while(std::getline(xs_data, line)){
				std::vector<double> tmp = parse_vector<double>(line);
				const double energy     = tmp[0];
				const double scatter    = tmp[1];
				const double capture    = tmp[2];
				const double fission    = tmp[3];
				const double nu         = tmp[4];
				const double total      = scatter + capture + fission;
				XS[XS_ENERGY].push_back(energy);
				XS[XS_TOTAL].push_back(total);
				XS[XS_CAPTURE].push_back(capture);
				XS[XS_SCATTERING].push_back(scatter);
				XS[XS_FISSION].push_back(fission);
				XS[XS_NU].push_back(nu);
			}
			
			// Create object
			the_nuclide = std::make_shared<Nuclide>(name, A, fissionable, a, b, g, XS);

			// Push to system
			nuclides.push_back(the_nuclide);
		}
	}
	
	//===========================================================================
	// Materials
	//===========================================================================
	
	for(const auto& node_material : input_materials.children("material")){
    std::shared_ptr<Material> the_material;

		// Inputs
		const std::string name    = node_material.attribute("name").value();
		const ull         ID      = node_material.attribute("id").as_ullong();
		const double      density = node_material.attribute("density").as_double();

		// Nuclides
		std::vector<std::shared_ptr<Nuclide>> nuclide_v;
		std::vector<double>                   ao_v;
		for(const auto& node_nuclide : node_material.children("nuclide")){
			// Inputs
			const std::string nuclide_name = node_nuclide.attribute("name").value();
			const double      nuclide_ao   = node_nuclide.attribute("ao").as_double();
			std::shared_ptr<Nuclide> the_nuclide = find_by_name(nuclides,nuclide_name);
			nuclide_v.push_back(the_nuclide);
			ao_v.push_back(nuclide_ao);
		}

		// Normalize
		double ao_total = 0.0;
		for(ull i = 0; i < nuclide_v.size(); i++) ao_total += ao_v[i];
		for(ull i = 0; i < nuclide_v.size(); i++) ao_v[i]  /= ao_total;

		// Density
		std::vector<double> density_v;
		double denom = 0.0;
		for(ull i = 0; i < nuclide_v.size(); i++) 
			denom += ao_v[i] * nuclide_v[i]->A;
		for(ull i = 0; i < nuclide_v.size(); i++)
			density_v.push_back(N_AVOGADRO * density / denom * ao_v[i] * 1.E-24);

		// The pair
		std::vector<std::pair<std::shared_ptr<Nuclide>,double>> nd;
		for(ull i = 0; i < nuclide_v.size(); i++)
			nd.push_back(std::make_pair(nuclide_v[i],density_v[i]));
		
		// Create object
		the_material = std::make_shared<Material>(name, ID, nd);

		// Push to system
		materials.push_back(the_material);
	}
}
