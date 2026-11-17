/*
  File: EWPOZ2.h
  Purpose: Declare cumulative NNLO+ SM pole form factors with typed inputs
           and safely owned FA/SW helpers that retain their order on copy.
  Merged-version authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-17
  Original GRIFFIN calculations/author records are retained below or in the
  including source file; the ownership repair preserves numerical expressions.
*/

/* EWPOZ2.h: header file for EWPOZ2.cc */

#ifndef __EWPOZ2__
#define __EWPOZ2__

#include "griffin/sm/EWPOZ.h"
#include "griffin/sm/delrho.h"

namespace griffin {

// The SMNNLO names denote cumulative NNLO+ results: the result() methods
// below also include the listed selected three- and four-loop terms.

// effective weak mixing angle predicted in the SM (at NNLO+)
class SW_SMNNLO : public SW_SMNLO {
protected:
  std::unique_ptr<SW_SMLO> clone() const override
  { return std::unique_ptr<SW_SMLO>(new SW_SMNNLO(*this)); }
public:
  using SW_SMNLO::SW_SMNLO;
  double drho2a2(void) const   // yt^4 corrections to \Delta\rho
  { return(sqr(ival->get(InputPar::WMassComplex)/ival->get(InputPar::ZMassComplex)) * delrho2a2(ival)); }
  double drho2aas(void) const  // yt^2*as corrections to \Delta\rho
  { return(sqr(ival->get(InputPar::WMassComplex)/ival->get(InputPar::ZMassComplex)) * delrho2aas(ival)); }
  double drho3a3(void) const   // yt^6 corrections to \Delta\rho
  { return(sqr(ival->get(InputPar::WMassComplex)/ival->get(InputPar::ZMassComplex)) * delrho3a3(ival)); }
  double drho3a2as(void) const // yt^4*as corrections to \Delta\rho
  { return(sqr(ival->get(InputPar::WMassComplex)/ival->get(InputPar::ZMassComplex)) * delrho3a2as(ival)); }
  double drho3aas2(void) const // yt^2*as^2 corrections to \Delta\rho
  { return(sqr(ival->get(InputPar::WMassComplex)/ival->get(InputPar::ZMassComplex)) * delrho3aas2(ival)); }
  double drho4aas3(void) const // yt^2*as^3 corrections to \Delta\rho
  { return(sqr(ival->get(InputPar::WMassComplex)/ival->get(InputPar::ZMassComplex)) * delrho4aas3(ival)); }
  double res2ff(void) const;  // 2-loop ew. corr. with two closed fermion loops  
  double res2fb(void) const;  // 2-loop ew. corr. with one closed fermion loop  
  double res2bb(void) const;  // 2-loop ew. corr. without closed fermion loops   
  double res2aas(void) const;  // 2-loop mixed ew-QCD corrections
  double res2aasnf(void) const;  // non-factorizable mixed ew-QCD corrections
  double res3fff(void) const;  // 3-loop ew. corr. with three closed fermion loops
  double res3ffa2as(void) const;  // 3-loop mixed ew-QCD corr. with two closed fermion loops
  Cplx result(void) const override
  {
    return(SW_SMNLO::result() + res2ff() + res2fb() + res2bb()
    	    + res2aas() + res2aasnf() + res3fff() + res3ffa2as()
	    + drho3aas2() + drho3a3() + drho3a2as() + drho4aas3());
  }
  Cplx errest(void) const override;   // error estimate from 1906.08815
};

// axial-vector form factor predicted in the SM (at NNLO+)
class FA_SMNNLO : public FA_SMNLO {
  double delrhofac(void) const;
protected:
  std::unique_ptr<FA_SMLO> clone() const override
  { return std::unique_ptr<FA_SMLO>(new FA_SMNNLO(*this)); }
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
  double res3fff(void) const;  // 3-loop ew. corr. with three closed fermion loops
  double res3ffa2as(void) const;  // 3-loop mixed ew-QCD corr. with two closed fermion loops
  Cplx result(void) const override
  {
    return(FA_SMNLO::result() + res2ff() + res2fb() + res2bb()
    	    + res2aas() + res2aasnf() + res3fff() + res3ffa2as() 
	    + drho3aas2() + drho3a3() + drho3a2as() + drho4aas3());
  }
  Cplx errest(void) const override;   // error estimate from 1906.08815
};

// Vector form factor predicted in the SM (at NNLO+); computed from F_A and
// sin^2(theta_eff), with an imaginary-ratio contribution first entering at NNLO.
class FV_SMNNLO : public FV_SMNLO {
public:
  // Allocate only the NNLO pair; inherited FV_SMLO storage owns both helpers.
  FV_SMNNLO(const Fermion type, const inval& input)
    : FV_SMNLO(type, input,
               makeFormFactor<FA_SMNNLO>(type, input),
               makeFormFactor<SW_SMNNLO>(type, input)) {}
  
  Cplx result(void) const override;
  Cplx errest(void) const override;   // error estimate from 1906.08815
};

} // namespace

#endif // __EWPOZ2__
