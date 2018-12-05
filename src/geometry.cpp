#include <cmath> 

#include "geometry.h"
#include "constant.h"
#include "error.h"


//=============================================================================
// S-Evaluate: check on which side of the surface a point is
//=============================================================================

double SurfaceSphere::evaluate(const Point& p) 
{
  const double x_t = p.x - x0;
  const double y_t = p.y - y0;
  const double z_t = p.z - z0;
  return x_t*x_t + y_t*y_t + z_t*z_t - rad_sq;
}


//=============================================================================
// S-Distance: Calculate distance for a particle to hit the surface. Return a
//             very large number if the particle is moving away or parallel to 
//             the surface.
//=============================================================================

double SurfaceSphere::distance( const Particle& P ) 
{
  Point p = P.pos;
  Point u = P.dir;

  // Put into quadratic equation form: a*s^2 + b*s + c = 0
	const double a = 1.0;
  const double b = 2.0 * ((p.x-x0)*u.x + (p.y-y0)*u.y + (p.z-z0)*u.z );
  const double c = evaluate(p);

  // Determinant
  const double D = b*b - 4.0 * a * c;
  	
  // Roots are complex, no intersection, return huge number
  // or identical roots, tangent, return huge number
  if(D <= 0.0) return MAX_float;
  else{
    const double sqrtD = std::sqrt(D);
    const double ai = 0.5 / a;

    // Roots
    double r1 = ai * (-1.0 * b - sqrtD);
    double r2 = ai * (-1.0 * b + sqrtD);

    // Negative roots return huge number (moving away from surface)
    if(r1 < 0) r1 = MAX_float;
    if(r2 < 0) r2 = MAX_float;

    return std::fmin( r1, r2 );
  }
}


//=============================================================================
// Cell-related functions
//=============================================================================

// Test if a point is inside a cell
bool test_point(const Point& p, const std::shared_ptr<Cell>& cell)
{
  // Loop over surfaces of cell, if particle not on correct side return false
  for(const auto& surface : cell->surfaces) {
    if(surface.first->evaluate(p) * surface.second < 0) return false;
  }
  return true;
}

// Search cell where a point belongs to
std::shared_ptr<Cell> 
search_cell(const Point& p, const std::vector<std::shared_ptr<Cell>>& cells)
{
  for(const auto& cell : cells) if(test_point(p,cell)) return cell;
	error("A particle is lost:\n  (x, y, z)  ("
			  +std::to_string(p.x)+", "
        +std::to_string(p.y)+", "
				+std::to_string(p.z)+")\n");
}
