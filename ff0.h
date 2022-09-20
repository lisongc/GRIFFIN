/* ff0.h: header file for LO form factors in classes.cc */

#include "Cplx.h"
#include "classes.h"

extern const double Qf[20];
extern const double I3f[20];

// tree-level axial-vector Z vertex factors
double az0(int type, const inval& input);
// tree-level vector Z vertex factors
double vz0(int type, const inval& input);
// tree-level Z vertex factors of form 'formt' (VEC/AXV)
double z0(int type, int formt, const inval& input);

// tree-level vector photon vertex factors
inline double vg0(int type, const inval& input) 
{ return -sqrt(4*Pi * input.get(al))*Qf[type]; }
// tree-level photon vertex factors of form 'formt' (VEC/AXV)
double g0(int type, int formt, const inval& input);

// returns real part of x if x is a regular number or zero otherwise
double realreg(Cplx x);
