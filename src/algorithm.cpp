#include "algorithm.h"
#include "constant.h"
#include "random.h"


//=============================================================================
// Binary search
//=============================================================================
// Note:
// 	value < lowest  grid --> -1
// 	value > highest grid --> vector.size - 1 (or number of bins)
// 	value = grid points  --> location of bin whose upper bound is the value
// 	                         (-1 if value = lowest grid)

ull binary_search(const double x, const std::vector<double>& vec)
{
	int left  = 0;
  int right = vec.size()- 1;
  int mid;

  while(left <= right){
    mid = (left + right)/2;        
    if(vec[mid] < x)left  = mid + 1;
    else             right = mid - 1;
  }
	
  return right; 
}


//=============================================================================
// Calculate XS
//=============================================================================
double calculate_XS(const ull idx, const double E, 
		                const std::vector<std::vector<double>>& XS, 
										const int header)
{
  // If incident energy exceed the table
  //   use the last data corresponding to the highest energy provided
  if(idx == XS[0].size()- 1)return XS[header].back();
    	
  // Similarly for energy below the lowest energy in the table
  else if(idx == -1)return XS[header][0];
    
	// Interpolate with the given bin index
  else{
    const double E1  = XS[XS_ENERGY][idx]; 
		const double E2  = XS[XS_ENERGY][idx+1];
		const double xs1 = XS[header][idx];    
		const double xs2 = XS[header][idx+1];
    return interpolate(E,E1,E2,xs1,xs2);
  }
}

//=============================================================================
// Interpolation
//=============================================================================

double interpolate(const double x, const double x1, const double x2,
                   const double y1, const double y2)
{ return (x - x2)/(x1 - x2)*y1 + (x - x1)/(x2 - x1)*y2; }


//=============================================================================
// Scatter direction
//=============================================================================

Point scatter_direction(const Point dir_i, const double mu0)
{
  // Sample azimuthal direction
  const double     azi = PI_2*Urand();
  const double cos_azi = std::cos(azi);
  const double sin_azi = std::sin(azi);
  const double      Ac = std::sqrt(1.0 - mu0*mu0);
  Point      dir_f; // Final direction

  if(dir_i.z != 1.0){
    const double       B = std::sqrt(1.0 - dir_i.z*dir_i.z);
    const double       C = Ac / B;
	
    dir_f.x = dir_i.x*mu0 
              + (dir_i.x*dir_i.z*cos_azi - dir_i.y*sin_azi)* C;
    dir_f.y = dir_i.y*mu0 
              + (dir_i.y*dir_i.z*cos_azi + dir_i.x*sin_azi)* C;
    dir_f.z = dir_i.z*mu0 - cos_azi*Ac*B;
  }
	
  // If dir_i = 0i + 0j + k, interchange z and y in the scattering formula
  else{
    const double B = std::sqrt(1.0 - dir_i.y*dir_i.y);
    const double C = Ac / B;
	
    Point q; // to store new direction point
    
    dir_f.x = dir_i.x*mu0 
              + (dir_i.x*dir_i.y*cos_azi - dir_i.z*sin_azi)* C;
    dir_f.z = dir_i.z*mu0 
              + (dir_i.z*dir_i.y*cos_azi + dir_i.x*sin_azi)* C;
    dir_f.y = dir_i.y*mu0 - cos_azi*Ac*B;
  }
  return dir_f;
}


//=============================================================================
// Watt spectrum
//=============================================================================

double watt_spectrum(const double E, const std::vector<double>& vec_a, 
		                 const std::vector<double>& vec_b, 
										 const std::vector<double>& vec_g)
{
  double a;
  double b;
  double g;
  double xi;   // xi_1 in formula
  double C;    // Acceptance parameter
  double Eout;

  // Binary search is not employed as there are only three grid points
  // E <= 1 eV (thermal)
  if(E <= 1.0){
		a = vec_a[0];
		b = vec_b[0];
		g = vec_g[0];
  }
  // 1 eV < E <= 1 MeV
  else if(E <= 1.0e6){
		a = interpolate(E, 1.0 , 1.0e6, vec_a[0], vec_a[1]);
		b = interpolate(E, 1.0 , 1.0e6, vec_b[0], vec_b[1]);
		g = interpolate(E, 1.0 , 1.0e6, vec_g[0], vec_g[1]);
  }
  // E >= 1 MeV, note for E > 14 MeV the values are extrapolated
  else{
		a = interpolate(E, 1.0e6, 14.0e6, vec_a[1], vec_a[2]);
		b = interpolate(E, 1.0e6, 14.0e6, vec_b[1], vec_b[2]);
		g = interpolate(E, 1.0e6, 14.0e6, vec_g[1], vec_g[2]);
  }
  
  do{
    	xi   = Urand();
			Eout = -a*g * std::log(xi); //MeV
      C    = (1.0 - g)* (1.0 - std::log(xi))- std::log(Urand());
  }
  while(C*C > b*Eout);
  
  return (Eout*1.0e6); //eV
}
