#ifndef CONSTANT_H
#define CONSTANT_H

#include <limits>
#include <cmath>

// Extreme numbers
const double MAX_float = std::numeric_limits<float>::max();

// Some constants
const double N_AVOGADRO = 6.0221409E23;
const double PI         = std::acos(-1.0);
const double PI_2       = 2.0 * PI;
const double PI_sqrt    = std::sqrt(PI);
const double PI_half    = 0.5 * PI;

// Surface boundary condition types
const int BC_VACUUM = 0;

// XS header
const int XS_ENERGY     = 0;
const int XS_TOTAL      = 1;
const int XS_CAPTURE    = 2;
const int XS_SCATTERING = 3;
const int XS_FISSION    = 4;
const int XS_NU         = 5;


#endif // CONSTANT_H
