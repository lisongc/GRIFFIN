/*
  File: SMval.h
  Purpose: Convert typed experimental (running-width) W/Z masses and widths
           into the complex-pole input values consumed by SM calculations.
           InputPar labels never change meaning with header inclusion order.
  Merged-version authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#ifndef __SMval__
#define __SMval__

#include "griffin/core/classes.h"

namespace griffin {

/* input parameter class for SM computations in the complex pole scheme;
   it computes the complex poles masses from user-provided masses in the
   running-width scheme */

class SMval : public inval {
protected:
  void compute(void)
  {
    data[index(InputPar::WMassComplex)] = data[index(InputPar::WMassExperiment)]/
        sqrt(1+sqr(data[index(InputPar::WWidthExperiment)]/data[index(InputPar::WMassExperiment)]));
    data[index(InputPar::WWidthComplex)] = data[index(InputPar::WWidthExperiment)]/
        sqrt(1+sqr(data[index(InputPar::WWidthExperiment)]/data[index(InputPar::WMassExperiment)]));
    data[index(InputPar::ZMassComplex)] = data[index(InputPar::ZMassExperiment)]/
        sqrt(1+sqr(data[index(InputPar::ZWidthExperiment)]/data[index(InputPar::ZMassExperiment)]));
    data[index(InputPar::ZWidthComplex)] = data[index(InputPar::ZWidthExperiment)]/
        sqrt(1+sqr(data[index(InputPar::ZWidthExperiment)]/data[index(InputPar::ZMassExperiment)]));
  }
public:
  using inval::inval;
  SMval(void) : inval() {};
  SMval(const inval& copyfrom) : inval(copyfrom) {};
};

} // namespace

#endif // __SMval__
