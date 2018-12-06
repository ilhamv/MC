#ifndef PARTICLE_H
#define PARTICLE_H

#include <memory>

#include "type.h"

struct Cell; // geometry.h

struct Point
{
  double x, y, z;
  Point(const double a = 0.0, const double b = 0.0, const double c = 0.0 ): 
    x(a), y(b), z(c) {};
};


struct Particle
{
  Point  pos;    // cm
  Point  dir;
  bool   alive;
  double energy; // eV
  std::shared_ptr<Cell> cell;

  Particle(const Point& p1, const Point& p2, const double E, 
           const std::shared_ptr<Cell> c): 
		alive(true), pos(p1), dir(p2), energy(E), cell(c) {};
};


struct Site
{
  Point  pos;        // cm
  Point  dir;
  double energy;     // eV
	ull    cell_index; // cell index in cells

  Site(const Point& p1, const Point& p2, const double E, const ull ci): 
		pos(p1), dir(p2), energy(E), cell_index(ci) {};
};


#endif // PARTICLE_H
