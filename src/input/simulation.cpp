#include "input.h"

void setup_simulation(const pugi::xml_node input_simulation, 
		                  ull& particles, ull& batches, ull& inactives)
{
	// Inputs
	particles = input_simulation.attribute("particles").as_ullong();
	batches   = input_simulation.attribute("batches").as_ullong();
	inactives = input_simulation.attribute("inactives").as_ullong();
}
