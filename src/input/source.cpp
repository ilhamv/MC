#include <mpi.h>

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

		//===========================================================================
		// Local source
		//===========================================================================
		
		// Get # of processes and id
		int np; MPI_Comm_size(MPI_COMM_WORLD, &np);
		int id; MPI_Comm_rank(MPI_COMM_WORLD, &id);

		// Local particles buffer
		ull particles_local = 1 + (particles-1) / np;

		// Process with excessive source_bank size
		ull excess = particles % np;
		if(excess > 0){
			const bool excessive_bank = id >= excess;
			if(excessive_bank) particles_local--;
		}

		source_bank.resize(particles_local,P);
	}
	else error("Source type \""+type+"\" is not supported.");
}
