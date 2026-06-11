/* SMEFTvalG.h: define input class that includes values for Wilsonian SMEFT coeffs using G_mu as input */

#ifndef __SMEFTvalG__
#define __SMEFTvalG__

#include "deltarSMEFT.h"
#include "SMvalG.h"

namespace griffin {

class SMEFTvalGmu : public invalGmuSMEFT {
protected: 
  void compute(void)
  {
    data[MZc] = data[MZ]/sqrt(1+sqr(data[GamZ]/data[MZ]));
    data[GZc] = data[GamZ]/sqrt(1+sqr(data[GamZ]/data[MZ]));
    invalGmuSMEFT::compute();
    data[GWc] = 0.3376186*data[Gmu]*powint(data[MWc],3)*
    		 (1+ 0.2122066*data[als]);
    data[MW] = data[MWc]*sqrt(1+sqr(data[GWc]/data[MWc]));
    data[GamW] = data[GWc]*sqrt(1+sqr(data[GWc]/data[MWc]));
  }
public:
  using invalGmuSMEFT::invalGmuSMEFT;
  SMEFTvalGmu(void) : invalGmuSMEFT(SIZE1) {};
  SMEFTvalGmu(const inval& copyfrom) : invalGmuSMEFT(copyfrom) { compute(); };
};

} // namespace

#endif // __SMEFTvalG__
