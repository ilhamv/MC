#include "input.h"


//=============================================================================
// Handler
//=============================================================================

void xml_handler(const std::string io_dir, ull& particles, ull& batches, 
 						     ull& inactives,
  							 std::vector<std::shared_ptr<Nuclide>>&  nuclides,
  							 std::vector<std::shared_ptr<Material>>& materials,
	 	  					 std::vector<std::shared_ptr<Surface>>&  surfaces,
  							 std::vector<std::shared_ptr<Cell>>&     cells,
								 std::vector<Particle>&                  source_bank)
{
  // XML input file
  const std::string input_name = io_dir + "input.xml";
  pugi::xml_document input_file;
  input_file.load_file(input_name.c_str());

	// Nodes
	pugi::xml_node input_simulation = input_file.child("simulation");
	pugi::xml_node input_materials  = input_file.child("materials");
	pugi::xml_node input_geometry   = input_file.child("geometry");
	pugi::xml_node input_source     = input_file.child("source");
	pugi::xml_node input_xs         = input_file.child("xs");

	// XS file
	const std::string xs_directory = input_xs.attribute("directory").value();

	// Set-up
	setup_simulation(input_simulation, particles, batches, inactives);
	setup_materials(input_materials, nuclides, materials, xs_directory);
	setup_geometry(input_geometry, surfaces, cells, materials);
	setup_source(input_source, source_bank, particles, cells);
}
