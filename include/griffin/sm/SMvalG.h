#ifndef GRIFFIN_SM_SMVALG_H
#define GRIFFIN_SM_SMVALG_H

/* SMvalG.h: define input class that converts PDG masses from/to complex-pole 
   masses for MW,MZ, using G_Fermi as input */

#include "griffin/sm/deltar.h"
#include "griffin/sm/SMval.h"

namespace griffin {

/* input class that computes MW and GamW from Gmu, 
   and in addition translates between running-width and complex pole masses */
class SMvalGmu : public invalGmu {
protected:
  void compute(void)
  {
    data[index(InputPar::ZMassComplex)] = data[index(InputPar::ZMassExperiment)]
        /sqrt(1+sqr(data[index(InputPar::ZWidthExperiment)]
                  /data[index(InputPar::ZMassExperiment)]));
    data[index(InputPar::ZWidthComplex)] = data[index(InputPar::ZWidthExperiment)]
        /sqrt(1+sqr(data[index(InputPar::ZWidthExperiment)]
                  /data[index(InputPar::ZMassExperiment)]));
    invalGmu::compute();
    data[index(InputPar::WWidthComplex)] =
        0.3376186*data[index(InputPar::fermiConstant)]
        *powint(data[index(InputPar::WMassComplex)],3)
        *(1+ 0.2122066*data[index(InputPar::alphaS)]);
    data[index(InputPar::WMassExperiment)] = data[index(InputPar::WMassComplex)]
        *sqrt(1+sqr(data[index(InputPar::WWidthComplex)]
                  /data[index(InputPar::WMassComplex)]));
    data[index(InputPar::WWidthExperiment)] = data[index(InputPar::WWidthComplex)]
        *sqrt(1+sqr(data[index(InputPar::WWidthComplex)]
                  /data[index(InputPar::WMassComplex)]));
  }
public:
  using invalGmu::invalGmu;
  SMvalGmu(void) : invalGmu() {};
  SMvalGmu(const inval& copyfrom) : invalGmu(copyfrom) { compute(); };
};

} // namespace

#endif // GRIFFIN_SM_SMVALG_H
