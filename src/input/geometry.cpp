#include <iostream>

#include "input.h"
#include "error.h"


// Get sign
template <typename T> int get_sign(T val) { return (T(0) < val)-(val < T(0)); }

void setup_geometry(const pugi::xml_node input_geometry, 
		                std::vector<std::shared_ptr<Surface>>& surfaces,
  							    std::vector<std::shared_ptr<Cell>>& cells,
  							    const std::vector<std::shared_ptr<Material>>& materials)
{
	//===========================================================================
	// Surfaces
	//===========================================================================

	for(const auto& node_surface : input_geometry.children("surface")){
    std::shared_ptr<Surface> the_surface;

		// Inputs
		const std::string name          = node_surface.attribute("name").value();
		const ull         ID            = node_surface.attribute("id").as_ullong();
		const std::string type          = node_surface.attribute("type").value();
		const std::string coeffs_string = node_surface.attribute("coeffs").value();
		const std::string BC_string     = node_surface.attribute("boundary")
			                                            .value();

		// Type: Sphere
		if(type == "sphere"){
			// BC
			int BC;
			if(BC_string == "vacuum") BC = 0;
			else error("Surface boundary \""+BC_string+"\" is not supported.");

			// Coeffs
			std::vector<double> coeffs = parse_vector<double>(coeffs_string);
			if(coeffs.size() != 4)
				error("Sphere surface requires coeffs of \"x0 y0 z0 R\"");
			const double x0 = coeffs[0];
			const double y0 = coeffs[1];
			const double z0 = coeffs[2];
			const double r  = coeffs[3];

			// Create object
			the_surface = std::make_shared<SurfaceSphere>(name,ID,BC,x0,y0,z0,r);
		}
		else error("Surface type \""+type+"\" is not supported.");

		// Push to system
		surfaces.push_back(the_surface);
	}
	
	//===========================================================================
	// Cells
	//===========================================================================

	for(const auto& node_cell : input_geometry.children("cell")){
    std::shared_ptr<Cell> the_cell;

		// Inputs
		const std::string name            = node_cell.attribute("name").value();
		const ull         ID              = node_cell.attribute("id").as_ullong();
		const std::string region_string   = node_cell.attribute("region").value();
		const std::string material_string = node_cell.attribute("material").value();

		// Cell surfaces
  	std::vector<std::pair<std::shared_ptr<Surface>,int>> cell_surfaces;
		std::vector<lol> region = parse_vector<lol>(region_string);
		for(const auto reg : region){
			const int sign       = get_sign(reg);
			const ull surface_ID = sign * reg;
			std::shared_ptr<Surface> the_surface = find_by_ID(surfaces,surface_ID);
			if(!the_surface)
				error("Surface with id \""+std::to_string(surface_ID)
						  +"\" is not defined.");
			cell_surfaces.push_back(std::make_pair(the_surface,sign));
		}

		// Cell material
		std::shared_ptr<Material> cell_material = find_by_name(materials,
				                                                   material_string);
		if(!cell_material)
			error("Material \""+material_string+"\" is not defined.");

		// Cell index in cells
		const ull index = cells.size();

		// Create object
		the_cell = std::make_shared<Cell>(name,ID,cell_surfaces,cell_material,
				                              index);

		// Push to system
		cells.push_back(the_cell);
	}
	
}
