#ifndef PHYSICS_H
#define PHYSICS_H

#include <vector>
#include <memory>

#include "particle.h"
#include "material.h"

void sample_reaction_scattering(Particle& P, const std::shared_ptr<Nuclide> N);
void sample_reaction_fission(Particle& P, const double nu_bar,
		                         const std::shared_ptr<Nuclide> N, 
		                         std::vector<Site>& fission_bank);

#endif // PHYSICS_H
