/* EWPOZSMEFT.h: header file for EWPOZ2.cc */

#ifndef __EWPOZMSEFT__
#define __EWPOZSMEFT__

#include "EWPOZ2.h"
#include "delrho.h"

namespace griffin {

// effective weak mixing angle predicted in the SMEFT (at LO)
class SW_SMEFTLO : public SW_SMNNLO {
public:
  using SW_SMNNLO::SW_SMNNLO;
  Cplx resSMEFTLO(void) const;  // tree-level SMEFT corrections
  Cplx result(void) const
  {
    return(SW_SMNNLO::result()+resSMEFTLO());
  }
};

// axial-vector form factor predicted in the SMEFT (at LO)
class FA_SMEFTLO : public FA_SMNNLO {
public:
  using FA_SMNNLO::FA_SMNNLO;
  Cplx resSMEFTLO(void) const;  // tree-level SMEFT corrections
  Cplx result(void) const
  {
    return(FA_SMNNLO::result()+resSMEFTLO());
  }
};

// vector form factor predicted in the SMEFT (at LO); computed from F_A and sw_eff
class FV_SMEFTLO : public FV_SMNNLO {
public:
  FV_SMEFTLO(const int type, const inval& input) : FV_SMNNLO(type, input)
  {
//    ftyp = type;
    fa = new FA_SMEFTLO(type, input);
    sw = new SW_SMEFTLO(type, input);
  }

  Cplx result(void) const;
};

} // namespace

#endif // __EWPOZMSEFT__
