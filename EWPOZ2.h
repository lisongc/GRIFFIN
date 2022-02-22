/* EWPOZ2.h: header file for EWPOZ2.cc */

#ifndef __EWPOZ2__
#define __EWPOZ2__

#include "EWPOZ.h"
#include "delrho.h"

// effective weak mixing angle predicted in the SM (at NNLO+)
class SW_SMNNLO : public SW_SMNLO {
public:
  using SW_SMNLO::SW_SMNLO;
  double drho2a2(void) const   // yt^4 corrections to \Delta\rho
  { return(sqr(ival->get(MW)/ival->get(MZ)) * delrho2a2(ival)); }
  double drho2aas(void) const  // yt^2*as corrections to \Delta\rho
  { return(sqr(ival->get(MW)/ival->get(MZ)) * delrho2aas(ival)); }
  double drho3a3(void) const   // yt^6 corrections to \Delta\rho
  { return(sqr(ival->get(MW)/ival->get(MZ)) * delrho3a3(ival)); }
  double drho3a2as(void) const // yt^4*as corrections to \Delta\rho
  { return(sqr(ival->get(MW)/ival->get(MZ)) * delrho3a2as(ival)); }
  double drho3aas2(void) const // yt^2*as^2 corrections to \Delta\rho
  { return(sqr(ival->get(MW)/ival->get(MZ)) * delrho3aas2(ival)); }
  double drho4aas3(void) const // yt^2*as^3 corrections to \Delta\rho
  { return(sqr(ival->get(MW)/ival->get(MZ)) * delrho4aas3(ival)); }
  double res2ff(void) const;  // 2-loop ew. corr. with two closed fermion loops  
  double res2fb(void) const;  // 2-loop ew. corr. with one closed fermion loop  
  double res2bb(void) const;  // 2-loop ew. corr. without closed fermion loops   
  double res2aas(void) const;  // 2-loop mixed ew-QCD corrections
  double res2aasnf(void) const;  // non-factorizable mixed ew-QCD corrections
// still to add: reducible 3-loop corrections from work with Lisong Chen
  Cplx result(void) const
  {
    return(SW_SMNLO::result() + res2ff() + res2fb() + res2bb()
    	    + res2aas() + drho3aas2() + drho3a3() + drho3a2as() + drho4aas3());
  }
};

// axial-vector form factor predicted in the SM (at NNLO+)
class FA_SMNNLO : public FA_SMNLO {
  double delrhofac(void) const;
public:
  using FA_SMNLO::FA_SMNLO;
  double drho2a2(void) const   // yt^4 corrections to \Delta\rho
  { return(delrhofac() * delrho2a2(ival)); }
  double drho2aas(void) const  // yt^2*as corrections to \Delta\rho
  { return(delrhofac() * delrho2aas(ival)); }
  double drho3a3(void) const   // yt^6 corrections to \Delta\rho
  { return(delrhofac() * delrho3a3(ival)); }
  double drho3a2as(void) const // yt^4*as corrections to \Delta\rho
  { return(delrhofac() * delrho3a2as(ival)); }
  double drho3aas2(void) const // yt^2*as^2 corrections to \Delta\rho
  { return(delrhofac() * delrho3aas2(ival)); }
  double drho4aas3(void) const // yt^2*as^3 corrections to \Delta\rho
  { return(delrhofac() * delrho4aas3(ival)); }
  double res2ff(void) const;  // 2-loop ew. corr. with two closed fermion loops  
  double res2fb(void) const;  // 2-loop ew. corr. with one closed fermion loop  
  double res2bb(void) const;  // 2-loop ew. corr. without closed fermion loops   
  double res2aas(void) const;  // 2-loop mixed ew-QCD corrections
  double res2aasnf(void) const;  // non-factorizable mixed ew-QCD corrections
// still to add: reducible 3-loop corrections from work with Lisong Chen
  Cplx result(void) const
  {
    return(FA_SMNLO::result() + res2ff() + res2fb() + res2bb()
    	    + res2aas() + res2aasnf() 
	    + drho3aas2() + drho3a3() + drho3a2as() + drho4aas3());
  }
};

// vector form factor predicted in the SM (at NNLO+)
class FV_SMNNLO : public FV_SMNLO {
public:
  FV_SMNNLO(const int type, const inval& input) : FV_SMNLO(type, input)
  {
//    ftyp = type;
    fa = new FA_SMNNLO(type, input);
    sw = new SW_SMNNLO(type, input);
  }
  
  Cplx result(void) const;
};

#endif // __EWPOZ2__
