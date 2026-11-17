#ifndef GRIFFIN_SM_FF0_H
#define GRIFFIN_SM_FF0_H

/* ff0.h: header file for LO form factors in classes.cc */

#include "griffin/core/Cplx.h"
#include "griffin/core/classes.h"

namespace griffin {

// electric charges of different fermion types, values are set in classes.cc
extern const double Qf[20];
// weak isospin of different fermion types, values are set in classes.cc
extern const double I3f[20];

// tree-level axial-vector Z vertex factors
double az0(Fermion type, const inval& input);
// tree-level vector Z vertex factors
double vz0(Fermion type, const inval& input);
// tree-level Z vertex factors of form 'formt' (Current::vector/Current::axial)
double z0(Fermion type, Current formt, const inval& input);

// tree-level vector photon vertex factors
inline double vg0(Fermion type, const inval& input)
{ return -sqrt(4*Pi * input.get(InputPar::alpha))*Qf[index(type)]; }
// tree-level photon vertex factors of form 'formt' (Current::vector/Current::axial)
double g0(Fermion type, Current formt, const inval& input);

// returns real part of x if x is a regular number or zero otherwise
double realreg(Cplx x);

} // namespace

#endif // GRIFFIN_SM_FF0_H
