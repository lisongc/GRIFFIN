/* xscnnlo.h: header file for xscnnlo.cc */

#include "classes.h"

// R @ NLO is same as R @ LO
#define R_SMLO coeffR
#define R_SMNLO coeffR
#define S_SMLO coeffS
#define Sp_SMLO coeffSp

class R_SMNNLO : public R_SMNLO {
protected:
  double cost;
public:
  R_SMNNLO(const int intype, const int outtype, const int inform, const int outform, 
  	   const double FAin, const double FAout, const double SWin, const double SWout,
	   const inval& input, const double costheta) : 
     R_SMNLO(intype, outtype, inform, outform, FAin, FAout, SWin, SWout, input)
  {
    cost = costheta;
  }
  
  R_SMNNLO(const int intype, const int outtype, const int inform, const int outform, 
  	   const psobs& FAin, const psobs& FAout, const psobs& SWin, const psobs& SWout,
	   const inval& input, const double costheta) : 
           R_SMNLO(intype, outtype, inform, outform, FAin, FAout, SWin, SWout, input)
  {
    cost = costheta;
  }
  
  Cplx result(void) const;
};

class S_SMNLO : public S_SMLO {
protected:
  double cost;
public:
  S_SMNLO(const int intype, const int outtype, const int inform, const int outform, 
  	  const inval& input, const double costheta) : 
	  S_SMLO(intype, outtype, inform, outform, input)
  {
    cost = costheta;
  }
  
  Cplx res1f(void) const;  // 1-loop corrections with closed fermion loops
  Cplx res1b(void) const;  // 1-loop corrections without closed fermion loops 
  Cplx result(void) const
  {
    return(S_SMLO::result()+res1f()+res1b());
  }
};


class msq_SMNNLO : public matelsq {
public:
  using matelsq::matelsq;
  Cplx result(void) const;
};
