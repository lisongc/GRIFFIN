/* SMvalG.h: define input class that converts PDG masses from/to complex-pole 
   masses for MW,MZ, using G_Fermi as input */

#include "deltar.h"
#include "SMval.h"

/* input class that computes MW and GamW from Gmu, 
   and in addition translates between running-with and complex pole masses */
class SMvalGmu : public invalGmu {
protected:
  void compute(void)
  {
    data[MZc] = data[MZ]/sqrt(1+sqr(data[GamZ]/data[MZ]));
    data[GZc] = data[GamZ]/sqrt(1+sqr(data[GamZ]/data[MZ]));
    invalGmu::compute();
    data[GWc] = 0.3376186*data[Gmu]*powint(data[MWc],3)*
    		 (1+ 0.2122066*data[als]);
    data[MW] = data[MWc]*sqrt(1+sqr(data[GWc]/data[MWc]));
    data[GamW] = data[GWc]*sqrt(1+sqr(data[GWc]/data[MWc]));
  }
public:
  using invalGmu::invalGmu;
  SMvalGmu(const inval& copyfrom) : invalGmu(copyfrom) {};
};

