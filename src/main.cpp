#include <iostream>
#include <fstream>
#include <cstring> 
#include <memory>
#include <vector>
#include <stack>
#include <cmath>

#include "input.h"
#include "error.h"
#include "random.h"
#include "constant.h"
#include "algorithm.h"
#include "physics.h"


int main( int argc, char* argv[] )
{
	//===========================================================================
  // Input set-up
	//===========================================================================
	
  // I/O Directory
  if(argc == 1) error("Please provide input.xml directory...");
  const std::string io_dir = std::string(argv[1]) + "/"; 

  // Simulation parameters
	ull particles, batches, inactives;
	
  // Particle banks
  std::vector<Particle> fission_bank, source_bank;
	std::stack<Particle>  secondary_bank;

  // System objects
  std::vector<std::shared_ptr<Nuclide>>  nuclides;
  std::vector<std::shared_ptr<Material>> materials;
  std::vector<std::shared_ptr<Surface>>  surfaces;
  std::vector<std::shared_ptr<Cell>>     cells;

	// Set-up
	xml_handler(io_dir, particles, batches, inactives, nuclides, materials,
			        surfaces, cells, source_bank);

	// Report setup
	std::cout<<"\nRunning with:\n"
		       <<"  Generations = "<<batches<<"\n"
		       <<"  Histories   = "<<particles<<"\n\n"
					 <<"Gen# - k\n";

	// k-eigenvalue
	std::vector<double> kgen(batches);
	double k = 1.0; // Current k
	double k_tally; // Tallying next k

	// RN seed stride for fission bank sampling
	const ull fbank_stride = 1 + (particles-1)/152917;

	//===========================================================================
	// LOOP 1: Generation
	//===========================================================================
	
	for(ull igen = 0; igen < batches; igen++){
		// Reset k tally
		k_tally = 0.0;

		//=========================================================================
		// LOOP 2: History
		//=========================================================================
	
		for(ull ihist = 0; ihist < particles; ihist++){
			// Initialize RN
			ull particle_number = particles * igen + ihist;
			RN_init_particle(&particle_number);

			// Get particle from source bank and put into secondary bank
			secondary_bank.push(source_bank[ihist]);

			//=======================================================================
			// LOOP 3: Secondaries
			//=======================================================================

			while(!secondary_bank.empty()){
				// Get particle from secondary as the current particle
				Particle P = secondary_bank.top();
				secondary_bank.pop();

				// XS cache of the current particle
				std::shared_ptr<Material> xs_material = nullptr;
				double                    xs_energy   = 0.0;
				double                    SigmaT, nuSigmaF;
				std::vector<double>       sigmaT, sigmaS, sigmaC, sigmaF, nu;

				//=====================================================================
				// LOOP 4: Particle random walk
				//=====================================================================
				
				while(P.alive){
					//===================================================================
					// XS data
					//===================================================================

					// Different material, resize vectors
					if(!xs_material || xs_material->ID != P.cell->material->ID){
						sigmaT.resize(P.cell->material->nuclides.size());
						sigmaS.resize(P.cell->material->nuclides.size());
						sigmaC.resize(P.cell->material->nuclides.size());
						sigmaF.resize(P.cell->material->nuclides.size());
						nu.resize(P.cell->material->nuclides.size());
						goto different_material;
					}

					// Different energy, binary search data
					if(xs_energy != P.energy){
						different_material:

						// Get current material
						std::shared_ptr<Material> M = P.cell->material;
	
						// Iterate over nuclides for microXS
						for(ull i = 0; i < M->nuclides.size(); i++){
							// Binary search index
							ull idx = binary_search(P.energy, 
									                    M->nuclides[i].first->XS[XS_ENERGY]);
							// MicroXS
							sigmaT[i] = calculate_XS(idx, P.energy, M->nuclides[i].first->XS,
									                     XS_TOTAL);
							sigmaS[i] = calculate_XS(idx, P.energy, M->nuclides[i].first->XS,
									                     XS_SCATTERING);
							sigmaC[i] = calculate_XS(idx, P.energy, M->nuclides[i].first->XS,
									                     XS_CAPTURE);
							sigmaF[i] = calculate_XS(idx, P.energy, M->nuclides[i].first->XS,
									                     XS_FISSION);
							nu[i]     = calculate_XS(idx, P.energy, M->nuclides[i].first->XS,
									                     XS_NU);
						}

						// Calculate macroXS
						SigmaT = 0.0; nuSigmaF = 0.0;
						for(ull i = 0; i < M->nuclides.size(); i++){
							SigmaT   +=  sigmaT[i]        * M->nuclides[i].second;
							nuSigmaF += (sigmaF[i]*nu[i]) * M->nuclides[i].second;
						}
					}

					//===================================================================
					// Move
					//===================================================================

					// Determine nearest surface and its distance
					std::shared_ptr<Surface> S;
					double                   dsurf;
					dsurf = MAX_float;
					for(const auto& s : P.cell->surfaces){
						double d = s.first->distance(P);
						if(d < dsurf){ 
							S     = s.first;
							dsurf = d;
						}
					}

					// Determine collision distance
					double dcoll = -std::log(Urand()) / SigmaT;

					// Move particle
					double dmove = std::min(dcoll,dsurf);
					P.pos.x += P.dir.x * dmove;
					P.pos.y += P.dir.y * dmove;
					P.pos.z += P.dir.z * dmove;

					// Tally k
					k_tally += nuSigmaF * dmove;
					

					//===================================================================
					// Surface hit
					//===================================================================
					
					if(dcoll > dsurf){
						if(S->BC == BC_VACUUM) P.alive = false;
					}
					
					//===================================================================
					// Collision
					//===================================================================
					
					else{
						// Get current material
						std::shared_ptr<Material> M = P.cell->material;
	
						// Sample nuclide
						double xi  = SigmaT * Urand();
						double sum = 0.0;
						ull    nuc_i;
						for(nuc_i = 0; nuc_i < M->nuclides.size(); nuc_i++){
							sum += sigmaT[nuc_i] * M->nuclides[nuc_i].second;
							if(sum > xi) break;
						}
						std::shared_ptr<Nuclide> N = nuclides[nuc_i];

						// Sample reaction
						xi  = sigmaT[nuc_i] * Urand();
						sum = 0.0;

						// Capture: kill particle
						sum += sigmaC[nuc_i];
						if(sum > xi){
							P.alive = false;
							continue;
						}

						// Scattering
						sum += sigmaS[nuc_i];
						if(sum > xi){
							sample_reaction_scattering(P, N);
							continue;
						}

						// Fission
						sum += sigmaF[nuc_i];
						if(sum > xi){
							sample_reaction_fission(P, nu[nuc_i]/k, N, fission_bank);
						}
					}
				} // LOOP 4: Particle random walk
			}	// LOOP 3: Secondaries
		} // LOOP 2: History

		// Calculate new k estimate
		kgen[igen] = k_tally/particles;

		// Update current k
		k = kgen[igen];

		// Report current k
		std::cout<<igen<<"  "<<kgen[igen]<<"\n";

		// Sample next source bank from fission bank
		ull sampling_seed = fbank_stride * igen;
		RN_init_particle(&sampling_seed);
		for(ull i = 0; i < particles; i++){
			ull idx = std::floor(Urand() * fission_bank.size());
			source_bank[i] = fission_bank[idx];
		}

		// Clear fission bank
		fission_bank.clear();

	} // LOOP 1: Generation
	
	// Output file
	std::ofstream ofile(io_dir+"k.txt");
	for(ull i = 0; i < batches; i++) ofile<<kgen[i]<<"\n";
}
