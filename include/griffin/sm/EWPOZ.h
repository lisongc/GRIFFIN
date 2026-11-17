/*
  File: EWPOZ.h
  Purpose: Declare NLO SM pole form factors with typed selections and
           safely owned, dynamically copied FA/SW helpers for FV.
  Merged-version authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-17
  Original GRIFFIN calculations/author records are retained below or in the
  including source file; the ownership repair preserves numerical expressions.
*/

/* EWPOZ.h: header file for EWPOZ.cc */

#ifndef __EWPOZ__
#define __EWPOZ__

#include "griffin/core/classes.h"

namespace griffin {

// effective weak mixing angle predicted in the SM (at NLO)
class SW_SMNLO : public SW_SMLO {
protected:
  std::unique_ptr<SW_SMLO> clone() const override
  { return std::unique_ptr<SW_SMLO>(new SW_SMNLO(*this)); }
public:
  using SW_SMLO::SW_SMLO;
  double res1f(void) const;  // 1-loop corrections with closed fermion loops
  double res1b(void) const;  // 1-loop corrections without closed fermion loops 
  Cplx result(void) const override
  {
    return(SW_SMLO::result()+res1f()+res1b());
  }
  Cplx errest(void) const override;   /* naive error estimate based on prefactors
  			        g^2/(4 Pi^2)*Nf and al_s/Pi*2*CF */

};

// axial-vector form factor predicted in the SM (at NLO)
class FA_SMNLO : public FA_SMLO {
protected:
  std::unique_ptr<FA_SMLO> clone() const override
  { return std::unique_ptr<FA_SMLO>(new FA_SMNLO(*this)); }
public:
  using FA_SMLO::FA_SMLO;
  double res1f(void) const;  // 1-loop corrections with closed fermion loops
  double res1b(void) const;  // 1-loop corrections without closed fermion loops 
  Cplx result(void) const override
  {
    return(FA_SMLO::result()+res1f()+res1b());
  }
  Cplx errest(void) const override;   /* naive error estimate based on prefactors
  			        g^2/(4 Pi^2)*Nf and al_s/Pi*2*CF */
};

// vector form factor predicted in the SM (at NLO); computed from F_A and sw_eff
class FV_SMNLO : public FV_SMLO {
protected:
  // FV_SMNNLO forwards its final-order pair without allocating an NLO pair.
  FV_SMNLO(const Fermion type, const inval& input,
           std::unique_ptr<FA_SMLO> axiaff, std::unique_ptr<SW_SMLO> swtheff)
    : FV_SMLO(type, input, std::move(axiaff), std::move(swtheff)) {}
public:
  // Allocate only the pair used in this prediction; FV_SMLO owns both.
  FV_SMNLO(const Fermion type, const inval& input)
    : FV_SMNLO(type, input,
               makeFormFactor<FA_SMNLO>(type, input),
               makeFormFactor<SW_SMNLO>(type, input)) {}
};

} // namespace

#endif // __EWPOZ__
