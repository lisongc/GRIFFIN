/* ff.h: header file for ff.cc */

#include "classes.h"

// real part of fermionic one-loop Z self-energy for k^2=mz^2
double rsz1f(const inval& input);
// imag. part of fermionic one-loop Z self-energy for k^2=mz^2
double isz1f(const inval& input);
// real part of bosonic one-loop Z self-energy for k^2=mz^2
double rsz1b(const inval& input);
// imag. part of bosonic one-loop Z self-energy for k^2=mz^2
double isz1b(const inval& input);
// deriv. of real part of fermionic one-loop Z self-energy for k^2=mz^2
double rsz1fp(const inval& input);
// deriv. of imag. part of fermionic one-loop Z self-energy for k^2=mz^2
double isz1fp(const inval& input);
// deriv. of real part of bosonic one-loop Z self-energy for k^2=mz^2
double rsz1bp(const inval& input);
// deriv. of imag. part of bosonic one-loop Z self-energy for k^2=mz^2
double isz1bp(const inval& input);
// 2nd deriv. of real part of fermionic one-loop Z self-energy for k^2=mz^2
double rsz1fpp(const inval& input);
// 2nd deriv. of imag. part of fermionic one-loop Z self-energy for k^2=mz^2
double isz1fpp(const inval& input);
// 2nd deriv. of real part of bosonic one-loop Z self-energy for k^2=mz^2
double rsz1bpp(const inval& input);
// 2nd deriv. of imag. part of bosonic one-loop Z self-energy for k^2=mz^2
double isz1bpp(const inval& input);

// real part of fermionic one-loop photon self-energy for k^2=mz^2
double rsg1f(const inval& input);
// imag. part of fermionic one-loop photon self-energy for k^2=mz^2
double isg1f(const inval& input);
// real part of bosonic one-loop photon self-energy for k^2=mz^2
double rsg1b(const inval& input);
// imag. part of bosonic one-loop photon self-energy for k^2=mz^2
double isg1b(const inval& input);

/* notation for the following functions:
   type = fermion type
   form = form factor type (VEC/AXV)
*/

// real part of fermionic one-loop Z vertex factor for k^2=mz^2
double rz1f(int type, int form, const inval& input);
// imag. part of fermionic one-loop Z vertex factor for k^2=mz^2
double iz1f(int type, int form, const inval& input);
// real part of bosonic one-loop Z vertex factor for k^2=mz^2
double rz1b(int type, int form, const inval& input);
// imag. part of bosonic one-loop Z vertex factor for k^2=mz^2
double iz1b(int type, int form, const inval& input);
// deriv. of real part of fermionic one-loop Z vertex factor for k^2=mz^2
double rz1fp(int type, int form, const inval& input);
// deriv. of imag. part of fermionic one-loop Z vertex factor for k^2=mz^2
double iz1fp(int type, int form, const inval& input);
// deriv. of real part of bosonic one-loop Z vertex factor for k^2=mz^2
double rz1bp(int type, int form, const inval& input);
// deriv. of imag. part of bosonic one-loop Z vertex factor for k^2=mz^2
double iz1bp(int type, int form, const inval& input);

// real part of fermionic one-loop photon vertex factor for k^2=mz^2
double rg1f(int type, int form, const inval& input);
// imag. part of fermionic one-loop photon vertex factor for k^2=mz^2
double ig1f(int type, int form, const inval& input);
// real part of bosonic one-loop photon vertex factor for k^2=mz^2
double rg1b(int type, int form, const inval& input);
// imag. part of bosonic one-loop photon vertex factor for k^2=mz^2
double ig1b(int type, int form, const inval& input);

// one-loop box diagrams for s=mz^2
// it/ot: ini-state/fin-state fermion type
// if1/of1: ini-state/fin-state form factor type (VEC/AXV)
// s, cost: kinematic variables
// inval: model input parameters
// AA, AZ: flags to turn on/off the gam-gam and gam-Z boxes (use with caution!)
// gzflag=1: use log(1-s/s0) for singular log term in gam-Z boxes
// 0<gzflag<<<1: use log(1-s/mz^2) for singular log term in gam-Z boxes
Cplx B1(int it, int ot, int if1, int of1, double s, double cost, const inval&
input, int AA, int AZ, double gzflag = 1);
// contribution from gamma-Z box contributing to R
double bRaz1(int it, int ot, double cost, const inval& input);


// real part of fermionic one-loop Z self-energy for k^2=s
double rsz1fs(double s, const inval& input);
// imag. part of fermionic one-loop Z self-energy for k^2=s
double isz1fs(double s, const inval& input);
// real part of bosonic one-loop Z self-energy for k^2=s
double rsz1bs(double s, const inval& input);
// imag. part of bosonic one-loop Z self-energy for k^2=s
double isz1bs(double s, const inval& input);
// deriv. of real part of fermionic one-loop Z self-energy for k^2=s

// real part of fermionic one-loop photon self-energy for k^2=s
double rsg1fs(double s, const inval& input);
// imag. part of fermionic one-loop photon self-energy for k^2=s
double isg1fs(double s, const inval& input);
// real part of bosonic one-loop photon self-energy for k^2=s
double rsg1bs(double s, const inval& input);
// imag. part of bosonic one-loop photon self-energy for k^2=s
double isg1bs(double s, const inval& input);

/* notation for the following functions:
   type = fermion type
   form = form factor type (VEC/AXV)
*/

// real part of fermionic one-loop Z vertex factor for k^2=s
double rz1fs(int type, int form, double s, const inval& input);
// imag. part of fermionic one-loop Z vertex factor for k^2=s
double iz1fs(int type, int form, double s, const inval& input);
// real part of bosonic one-loop Z vertex factor for k^2=s
double rz1bs(int type, int form, double s, const inval& input);
// imag. part of bosonic one-loop Z vertex factor for k^2=s
double iz1bs(int type, int form, double s, const inval& input);

// real part of fermionic one-loop photon vertex factor for k^2=s
double rg1fs(int type, int form, double s, const inval& input);
// imag. part of fermionic one-loop photon vertex factor for k^2=s
double ig1fs(int type, int form, double s, const inval& input);
// real part of bosonic one-loop photon vertex factor for k^2=s
double rg1bs(int type, int form, double s, const inval& input);
// imag. part of bosonic one-loop photon vertex factor for k^2=s
double ig1bs(int type, int form, double s, const inval& input);

// one-loop box diagrams for arbitrary s
// it/ot: ini-state/fin-state fermion type
// if1/of1: ini-state/fin-state form factor type (VEC/AXV)
// s, cost: kinematic variables
// inval: model input parameters
// AA, AZ: flags to turn on/off the gam-gam and gam-Z boxes (use with caution!)
Cplx B1s(int it, int ot, int if1, int of1, double s, double cost, const inval&
input, int AA, int AZ);
