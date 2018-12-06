#include <cmath>

#include "physics.h"
#include "random.h"
#include "constant.h"
#include "algorithm.h"
#include "geometry.h"

//=============================================================================
// Scattering Reaction
//=============================================================================
// Scatter with sampled scattering angle mu0, with nucleus mass A
// Scattering is trated in Center of mass (COM) frame
// Current model: Free gas scattering with constant cross section

void sample_reaction_scattering(Particle& P, const std::shared_ptr<Nuclide> N)
{
  const double mu0 = 2.0 * Urand() - 1.0; // Isotropic scattering

	const double A = N->A;
	
	//===========================================================================
  // Sampling nuclide velocity
	//===========================================================================
	
  double V_tilda;  // Nuclide speed candidate
  double mu_tilda; // Nuclide-neutron polar cosine candidate

  const double beta = std::sqrt(2.0659834e-11 * A); // Eq. 19
  // note, the constant above is 
  // (1.674927471e-27 kg) / (1.38064852e-19 cm^2 kg s^-2 K^-1) / (293.6 K)/2
  // 	 293.6 comes from room temperature
	
	// Particle speed
	double P_speed = 13831.5926439 * std::sqrt(P.energy) * 100.0;
  // constant: 2.0 * (1.60217662e-19 J/eV) / (1.674927471e-27 kg)
	
  const double y = beta * P_speed; // Eq. 32
	
  // Sample candidate V_tilda and mu_tilda
  do{
		double x;
		if(Urand() < 2.0/(2.0 + PI_sqrt*y)) // Eq. 37
	    // w2 -> sample g2
	    x = std::sqrt(-std::log(Urand()*Urand())); // Eq. 39
		else{
	    // w1 --> sample g1
	    const double cos_val = std::cos(PI_half*Urand());
	    x = std::sqrt(-std::log(Urand()) - std::log(Urand())*cos_val*cos_val); 
			// Eq. 38
		}
		V_tilda  = x / beta; // Eq. 32
		mu_tilda = 2.0*Urand() - 1.0; // Eq. 40
  }	
  // Accept candidate V_tilda and mu_tilda?
  while(Urand() > std::sqrt(P_speed*P_speed + V_tilda*V_tilda 
				                    - 2.0 * P_speed * V_tilda * mu_tilda)
                  /(P_speed + V_tilda)); // Eq. 41
  // Nuclide direction, Eq. 42
  Point nuclide_dir = scatter_direction(P.dir, mu_tilda); 
  // Nuclide velocity - LAB, Eq. 43
  Point V_lab(nuclide_dir.x*V_tilda, nuclide_dir.y*V_tilda, 
			        nuclide_dir.z*V_tilda); 


	//===========================================================================
  // COM Kinematics
	//===========================================================================

  // Particle velocity - LAB
  Point v_lab(P_speed*P.dir.x, P_speed*P.dir.y, P_speed*P.dir.z);
		
  // COM velocity	
  const Point u((v_lab.x + A*V_lab.x)/(1.0+A), 
                (v_lab.y + A*V_lab.y)/(1.0+A), 
                (v_lab.z + A*V_lab.z)/(1.0+A)); // Eq. 6
	
  // Particle velocity - COM
  Point v_c(v_lab.x-u.x, v_lab.y-u.y, v_lab.z-u.z);
	
  // Particle speed - COM
  const double speed_c = std::sqrt(v_c.x*v_c.x + v_c.y*v_c.y + v_c.z*v_c.z);

  // Particle initial direction - COM
  const Point dir_c(v_c.x/speed_c, v_c.y/speed_c, v_c.z/speed_c);
	
  // Scattering the direction in COM
  Point dir_cNew = scatter_direction(dir_c, mu0); // Final direction - COM

  // Final velocity - COM
  v_c.x = speed_c * dir_cNew.x;
  v_c.y = speed_c * dir_cNew.y;
  v_c.z = speed_c * dir_cNew.z;
	
	//===========================================================================
  // Convert to LAB
	//===========================================================================

  // Final velocity - LAB
  v_lab.x = v_c.x + u.x;
  v_lab.y = v_c.y + u.y;
  v_lab.z = v_c.z + u.z;

  // Final speed - LAB
  P_speed  = std::sqrt(v_lab.x*v_lab.x + v_lab.y*v_lab.y + v_lab.z*v_lab.z);
  P.energy = 5.2270376e-13 * P_speed * P_speed;
  // constant: 0.5 / (1.60217662e-19 J/eV) * (1.674927471e-27 kg) 
  //           / (10000 cm^2/m^2)

  // Final direction - LAB
  P.dir = Point(v_lab.x/P_speed, v_lab.y/P_speed, v_lab.z/P_speed);
}


//=============================================================================
// Fission Reaction
//=============================================================================
void sample_reaction_fission(Particle& P, const double nu_bar,
		                         const std::shared_ptr<Nuclide> N, 
		                         std::vector<Site>& fission_bank)
{
	// Absorb/kill particle
	P.alive = false;

	// Sample number of fission neutrons
	const int nu = std::floor(nu_bar + Urand());
		
	// Sample fission neutron
	for(int i = 0; i < nu; i++){
		// Energy
		double E = watt_spectrum(P.energy, N->a, N->b, N->g);

		// Direction
    const double mu  = 2.0 * Urand() - 1.0;
    const double azi = PI_2 * Urand();
    const double c   = std::sqrt(1.0 - mu * mu);
		Point dir(mu, std::cos(azi) * c, std::sin(azi) * c);

		// Push fission neutron to fission bank
		Site P_fission(P.pos, dir, E, P.cell->index);
		fission_bank.push_back(P_fission);
	}
}
