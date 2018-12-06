#include <iostream>
#include <fstream>
#include <cstring> 
#include <memory>
#include <vector>
#include <stack>
#include <cmath>
#include <mpi.h>

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
	
  // MPI initialization
	MPI_Init(&argc, &argv);
  MPI_Status status;
	  
	// Get # of processes and id
	int np; MPI_Comm_size(MPI_COMM_WORLD, &np);
	int id; MPI_Comm_rank(MPI_COMM_WORLD, &id);

  // I/O Directory
	if(id == 0)
		if(argc == 1) error("Please provide input.xml directory...");
	MPI_Barrier(MPI_COMM_WORLD);
	const std::string io_dir = std::string(argv[1]) + "/";

  // Simulation parameters
	ull particles, batches, inactives;
	
  // Particle banks
  std::vector<Particle> source_bank;
	std::stack<Particle>  secondary_bank;
  std::vector<Site>     fission_bank;
  std::vector<Site>     sample_bank;

  // System objects
  std::vector<std::shared_ptr<Nuclide>>  nuclides;
  std::vector<std::shared_ptr<Material>> materials;
  std::vector<std::shared_ptr<Surface>>  surfaces;
  std::vector<std::shared_ptr<Cell>>     cells;

	// Set-up
	xml_handler(io_dir, particles, batches, inactives, nuclides, materials,
			        surfaces, cells, source_bank);

	// Report setup
	if(id == 0)
		std::cout<<"\nRunning with:\n"
						 <<"  Generations = "<<batches<<"\n"
						 <<"  Histories   = "<<particles<<"\n\n"
						 <<"Gen# - k\n";

	// k-eigenvalue
	std::vector<double> kgen(batches);
	double k = 1.0;  // Current k
	double k_local;  // Tallying next k
	double k_global; // root only

	// RN seed stride for fission bank sampling
	const ull fbank_stride = 1 + (particles-1)/152917;
	

	//===========================================================================
	// Parallel preparation
	//===========================================================================
	
	// Status and request
  MPI_Status  status_left, status_right;
  MPI_Request request_left, request_right;

	// Neighbors
	int left  = id-1;
	int right = id+1;

	// Local particles per generation
	const ull particles_local = source_bank.size();

	// Starting and ending indexes of source, fission, and sample bank
	ull source_start, source_end, fission_start, fission_end, sample_start, 
			sample_end;

	// Source bank indexes
	MPI_Exscan(&particles_local, &source_start, 1, MPI_UNSIGNED_LONG_LONG, 
			       MPI_SUM, MPI_COMM_WORLD);
	source_end = source_start + particles_local - 1;

  // Set up MPI Point type
  MPI_Datatype MPI_POINT;
	MPI_Datatype point_type[1]     = {MPI_DOUBLE};
	int          point_blocklen[1] = {3};
	MPI_Aint     point_disp[1]     = {0};
	MPI_Type_create_struct(1, point_blocklen, point_disp, point_type, 
			                   &MPI_POINT);
	MPI_Type_commit(&MPI_POINT);

  // Set up MPI Site type
  MPI_Datatype MPI_SITE;
	MPI_Datatype site_type[3]     = {MPI_POINT, MPI_DOUBLE, 
		                               MPI_UNSIGNED_LONG_LONG};
	int          site_blocklen[3] = {2,1,1};
	MPI_Aint     site_disp[3], extent, dummy;
	site_disp[0] = 0;
	MPI_Type_get_extent(MPI_POINT, &dummy, &extent);
	site_disp[1] = extent * 2;
	MPI_Type_get_extent(MPI_DOUBLE, &dummy, &extent);
	site_disp[2] = extent * 1;
	MPI_Type_create_struct(3, site_blocklen, site_disp, site_type, &MPI_SITE);
	MPI_Type_commit(&MPI_SITE);


	//===========================================================================
	// LOOP 1: Generation
	//===========================================================================

	for(ull igen = 0; igen < batches; igen++){
		// Reset k tally
		k_local = 0.0;

		//=========================================================================
		// LOOP 2: History
		//=========================================================================
	
		for(ull ihist = 0; ihist < particles_local; ihist++){
			// Initialize RN
			ull particle_number = particles * igen + source_start + ihist;
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
					k_local += nuSigmaF * dmove;
					

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

		//=========================================================================
		// Update k
		//=========================================================================
	
		MPI_Barrier(MPI_COMM_WORLD);

		// Calculate new k estimate
  	MPI_Allreduce(&k_local, &k_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		kgen[igen] = k_global/particles;

		// Update current k
		k = kgen[igen];

		// Report current k
		if(id == 0)
			std::cout<<igen<<"  "<<kgen[igen]<<"\n";
		MPI_Barrier(MPI_COMM_WORLD);

		//=========================================================================
		// Fission bank operations
		//=========================================================================

		// Fission bank indexes
		const ull fissions_local = fission_bank.size();
		MPI_Exscan(&fissions_local, &fission_start, 1, MPI_UNSIGNED_LONG_LONG, 
				       MPI_SUM, MPI_COMM_WORLD);
		fission_end = fission_start + fissions_local - 1;

		// Locally sample the fission bank
		ull sampling_seed = fbank_stride * igen;
		RN_init_particle(&sampling_seed);
		for(ull i = 0; i < particles; i++){
			ull idx = std::floor(Urand() * fission_bank.size());
			if(fission_start <= idx && idx <= fission_end)
				sample_bank.push_back(fission_bank[idx]);
		}

		// Sample bank index
		ull       sample_start, sample_end;
		const ull samples_local = sample_bank.size();

		MPI_Exscan(&samples_local, &sample_start, 1, MPI_UNSIGNED_LONG_LONG, 
				       MPI_SUM, MPI_COMM_WORLD);
		sample_end = sample_start + samples_local - 1;

		// Send/recv sampled fission bank if necessary
		bool more_left, less_left, more_right, less_right;
		more_left  = sample_start < source_start;
		less_left  = sample_start > source_start;
		more_right = sample_end   > source_end;
		less_right = sample_end   < source_end;

		ull post_samples      = particles_local;
		ull post_sample_start = 0;
		ull post_sample_shift = 0;

		if(more_left){
			ull n = source_start - sample_start;
      MPI_Isend(&sample_bank[0], n, MPI_SITE, left, 0, MPI_COMM_WORLD, 
					      &request_left);
			post_sample_shift = n;
		}
		if(less_left){
			ull n = sample_start - source_start;
      MPI_Recv(&source_bank[0], n, MPI_SITE, right, 0, MPI_COMM_WORLD, 
					     &status_left);
			post_samples -= n;
			post_sample_start = n;
		}
		if(more_right){
			ull n = sample_end - source_end;
      MPI_Isend(&sample_bank[samples_local-n], n, MPI_SITE, right, 0, 
					      MPI_COMM_WORLD, &request_right);
		}
		if(less_right){
			ull n = source_end - sample_end;
      MPI_Recv(&source_bank[particles_local-n], n, MPI_SITE, right, 0, 
					      MPI_COMM_WORLD, &status_right);
			post_samples -= n;
		}

    if(more_left  || less_left)  MPI_Wait(&request_left, &status_left);
    if(more_right || less_right) MPI_Wait(&request_right, &status_right);

		// Accordingly copy post sample bank to source bank
		for(int i = post_sample_start; i < post_samples; i++){
			source_bank[i].pos    = sample_bank[i+post_sample_shift].pos;
			source_bank[i].dir    = sample_bank[i+post_sample_shift].dir;
			source_bank[i].energy = sample_bank[i+post_sample_shift].energy;
			source_bank[i].cell   = cells[sample_bank[i+post_sample_shift]
				                                          .cell_index];
		}

		// Clear fission and sample bank
		fission_bank.clear();
		sample_bank.clear();

	} // LOOP 1: Generation
	
	// Output file
	std::ofstream ofile(io_dir+"k.txt");
	for(ull i = 0; i < batches; i++) ofile<<kgen[i]<<"\n";
	
	MPI_Finalize();
	return 0;
}
