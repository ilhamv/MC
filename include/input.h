#ifndef INPUT_H
#define INPUT_H

#include <vector>
#include <sstream>
#include "pugixml.hpp"

#include "geometry.h"
#include "material.h"
#include "particle.h"

void xml_handler(const std::string io_dir, ull& particles, ull& batches, 
						     ull& inactives,
  							 std::vector<std::shared_ptr<Nuclide>>&  nuclides,
  							 std::vector<std::shared_ptr<Material>>& materials,
		  					 std::vector<std::shared_ptr<Surface>>&  surfaces,
  							 std::vector<std::shared_ptr<Cell>>&     cells,
								 std::vector<Particle>&                  source_bank);

void setup_simulation(const pugi::xml_node input_simulation, 
		                  ull& particles, ull& batches, ull& inactives);
void setup_materials(const pugi::xml_node input_materials, 
		                 std::vector<std::shared_ptr<Nuclide>>& nuclides,
  							     std::vector<std::shared_ptr<Material>>& materials,
										 const std::string xs_directory);
void setup_geometry(const pugi::xml_node input_geometry, 
		                std::vector<std::shared_ptr<Surface>>& surfaces,
  							    std::vector<std::shared_ptr<Cell>>& cells,
  							    const std::vector<std::shared_ptr<Material>>& materials);
void setup_source(const pugi::xml_node input_source, 
		              std::vector<Particle>& source_bank, const ull particles,
									const std::vector<std::shared_ptr<Cell>>& cells);


//=============================================================================
// Some functions for search and parsing
//=============================================================================

template<typename T>
std::shared_ptr<T> find_by_ID(
		const std::vector<std::shared_ptr<T>>& vec, const ull ID)
{
	for(auto& v : vec){
	if(v->ID == ID) return v;
  }
  return nullptr;
}

template<typename T>
std::shared_ptr<T> find_by_name(
	const std::vector<std::shared_ptr<T>>& vec, const std::string name)
{
	for(auto& v : vec){
	if(v->name == name) return v;
	}
	return nullptr;
}

template<class T>
std::vector<T> parse_vector(std::string const& pointLine)
{
  std::istringstream iss(pointLine);

  return std::vector<T>{std::istream_iterator<T>(iss), 
		                    std::istream_iterator<T>()};
}

#endif // INPUT_H
