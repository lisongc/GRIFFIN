/* deltar.h: header file for deltar.cc */

#ifndef __deltar__
#define __deltar__

#include "classes.h"
#include "delrho.h"

// Delta r predicted in the SM (at NLO)
class dr_SMNLO : public psobs {
public:
  using psobs::psobs;
  double res1f(void) const;  // 1-loop corrections with closed fermion loops
  double res1b(void) const;  // 1-loop corrections without closed fermion loops 
  Cplx result(void) const
  {
    return(res1f()+res1b());
  }
  double compmw(void) const; // compute m_W using G_mu as input
};


// Delta r predicted in the SM (at NNLO+)
class dr_SMNNLO : public dr_SMNLO {
public:
  using dr_SMNLO::dr_SMNLO;
  double drho2a2(void) const   // yt^4 corrections to \Delta\rho
  { return(1/(1-sqr(ival->get(MZ)/ival->get(MW))) * delrho2a2(ival)); }
  double drho2aas(void) const  // yt^2*as corrections to \Delta\rho
  { return(1/(1-sqr(ival->get(MZ)/ival->get(MW))) * delrho2aas(ival)); }
  double drho3a3(void) const   // yt^6 corrections to \Delta\rho
  { return(1/(1-sqr(ival->get(MZ)/ival->get(MW))) * delrho3a3(ival)); }
  double drho3a2as(void) const // yt^4*as corrections to \Delta\rho
  { return(1/(1-sqr(ival->get(MZ)/ival->get(MW))) * delrho3a2as(ival)); }
  double drho3aas2(void) const // yt^2*as^2 corrections to \Delta\rho
  { return(1/(1-sqr(ival->get(MZ)/ival->get(MW))) * delrho3aas2(ival)); }
  double drho4aas3(void) const // yt^2*as^3 corrections to \Delta\rho
  { return(1/(1-sqr(ival->get(MZ)/ival->get(MW))) * delrho4aas3(ival)); }
  double res2ff(void) const;  // 2-loop ew. corr. with two closed fermion loops  
  double res2fb(void) const;  // 2-loop ew. corr. with one closed fermion loop  
  double res2bb(void) const;  // 2-loop ew. corr. without closed fermion loops   
  double res2aas(void) const;  // 2-loop mixed ew-QCD corrections
  double res3aas2(void) const; // 3-loop O(al*als^2) corrections
  double res3fff(void) const;  // 3-loop ew. corr. with three closed fermion loops
  double res3ffa2as(void) const;  // 3-loop mixed ew-QCD corr. with two closed fermion loops
  Cplx result(void) const
  {
    return(dr_SMNLO::result() + res2ff() + res2fb() + res2bb()
    	    + res2aas() + res3aas2() + res3fff() + res3ffa2as()
	    + drho3a3() + drho3a2as() + drho4aas3());
  }
  Cplx errest(void) const
  {
    return 0.00024;  // from hep-ph/0311148
  }
};


/* input class that computes MW from Gmu */
class invalGmu : public inval {
protected:
  void compute(void);
public:
  using inval::inval;
  invalGmu(const inval& copyfrom) : inval(copyfrom) { compute(); };
};

#endif // __deltar__
