#include "input.h"
#include "error.h"

void setup_source(const pugi::xml_node input_source, 
		              std::vector<Particle>& source_bank, const ull particles, 
									const std::vector<std::shared_ptr<Cell>>& cells)
{
	// Inputs
	const std::string type          = input_source.attribute("type").value();
	const std::string params_string = input_source.attribute("params").value();

	if(type == "delta"){
		// Params
		std::vector<double> params = parse_vector<double>(params_string);
		if(params.size() != 7)
			error("Delta source requires params of \"x y z u v w E\"");
		const double x = params[0];
		const double y = params[1];
		const double z = params[2];
		const double u = params[3];
		const double v = params[4];
		const double w = params[5];
		const double E = params[6];

		// Sample source
		const Point pos                  = Point(x,y,z);
		const Point dir                  = Point(u,v,w);
		const std::shared_ptr<Cell> cell = search_cell(pos, cells);
		Particle P(pos, dir, E, cell);

		source_bank.resize(particles,P);
	}
	else error("Source type \""+type+"\" is not supported.");
}
