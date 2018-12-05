#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <memory> 
#include <vector> 
#include <cstring>

#include "particle.h"
#include "type.h"

struct Material; // material.h

//=============================================================================
// Surface
//=============================================================================

struct Surface
{
  // Properties
  const std::string name;
  const ull         ID;
	const int         BC;
  
  Surface(const std::string n, const ull i, const int bc): 
		name(n), ID(i), BC(bc) {};
	
  // Evaluate: if particle is on the left or right side
  virtual double evaluate(const Point& p) = 0;

  // Distance for particle to hit the surface
  virtual double distance(const Particle& P) = 0;
};

// Sphere
struct SurfaceSphere : public Surface 
{
  const double x0, y0, z0, rad, rad_sq;

  SurfaceSphere(const std::string n, const ull i, const int bc, const double p1,
					      const double p2, const double p3, const double p4): 
		Surface(n,i,bc), x0(p1), y0(p2), z0(p3), rad(p4), rad_sq(p4*p4) {};

  double evaluate(const Point& p);
  double distance(const Particle& P);
};



//=============================================================================
// Cell
//=============================================================================

struct Cell
{
  // Properties
  const std::string                                          name;
  const ull                                                  ID;
  const std::vector<std::pair<std::shared_ptr<Surface>,int>> surfaces;
  const std::shared_ptr<Material>                            material;

  Cell(const std::string n, const ull i, 
       const std::vector<std::pair<std::shared_ptr<Surface>,int>> s, 
       const std::shared_ptr<Material> m): 
		name(n), ID(i), surfaces(s), material(m) {};
};

// Test if a point is inside a cell
bool test_point(const Point& p, const std::shared_ptr<Cell>& cell);
	
// Search cell where a point belongs to
std::shared_ptr<Cell>
search_cell(const Point& p, const std::vector<std::shared_ptr<Cell>>& cells);

#endif // GEOMETRY_H
