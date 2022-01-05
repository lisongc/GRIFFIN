//#include "Cplx.h"
#include "classes.h"
#include "oneloop.h"

extern const double Qf[7];
extern const double I3f[7];

// tree-level axial-vector Z vertex factors
double az0(int type, const inval &input);
// tree-level vector Z vertex factors
double vz0(int type, const inval &input);
// tree-level Z vertex factors of form 'formt' (VEC/AXV)
double z0(int type, int formt, const inval &input);

// tree-level vector photon vertex factors
//inline double vg0(int type, double el) { return -el * Qf[type]; }
// tree-level photon vertex factors of form 'formt' (VEC/AXV)
double g0(int type, int formt, const inval &input);

//one-loop Z-boson self energy and derivatives
Cplx SigZ1(const inval &input);
Cplx SigZ1p(const inval &input); //first order derivative
Cplx SigZ1p2(const inval &input); //second order derivative

//one-loop gamma-gamma self energy
Cplx SigG1(const inval &input);
// one-loop Z-vertex factor
Cplx Z1(const int ftyp, const int off, const inval &input);

//one-loop Z-derivative factor
Cplx Z1p(const int ftyp, const int off, const inval &input);

//one-loop G-vertex factor
Cplx G1(const int ftyp, const int off, const inval &input);

//one-loop regular contribution of gauge-boson mediated boxes (S)
 Cplx B1(const int ftyp, const int iff, const int off, const int GB, const inval &input, double s);

 //one-loop principal contribution of gauge-boson mediated boxes (R)
 Cplx bgzR(const int it, const int ot, const inval &input);

 //auxiliary functions 
 Cplx deX2(const int it, const int ot, const inval &input);
 Cplx rAAI(const int it, const int ot, const inval &input);
 Cplx Ixf(const int form, const int ftyp, const inval &input);
 Cplx xIij(const int it, const int ot, const int iff, const int off, const inval &input);
 //Cplx QWf(const int form, const int ftyp, const psobs &SWf, const inval &input);
