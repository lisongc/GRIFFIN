// Generate the frozen reference with the unmodified upstream v1.2 library.
// See README.md in this directory for provenance and reproduction commands.
#include "EWPOZ2.h"
#include "xscnnlo2.h"
#include "SMval.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace griffin;

int main()
{
  SMval input;
  input.set(MZ,91.1876);
  input.set(MW,80.377);
  input.set(GamZ,2.4952);
  input.set(GamW,2.085);
  input.set(al,1.0/137.03599976);
  input.set(als,0.1179);
  input.set(MH,125.1);
  input.set(MT,172.5);
  input.set(MB,2.87);
  input.set(MC,1.27);
  input.set(ML,1.777);
  input.set(MS,0.342);
  input.set(MU,0.246);
  input.set(MD,0.246);
  input.set(MM,0.1057);
  input.set(ME,0.000511);
  input.set(Delal,0.059);

  const int ci[] = {VEC,AXV,VEC,AXV};
  const int cf[] = {VEC,VEC,AXV,AXV};
  const double energies[] = {50.0,90.0,160.0};
  const double cosines[] = {-0.45,0.3,0.6};
  std::cout << std::scientific
            << std::setprecision(std::numeric_limits<double>::max_digits10);
  for(int final : {MUO,UQU,DQU})
  {
    const FA_SMNNLO fai(ELE,input), faf(final,input);
    const SW_SMNNLO swi(ELE,input), swf(final,input);
    for(unsigned p=0;p<3;++p)
      for(unsigned k=0;k<4;++k)
      {
        const double s=energies[p]*energies[p];
        const mat_SMNNLO2 upgraded(ELE,final,ci[k],cf[k],
            fai,faf,swi,swf,s,cosines[p],input);
        const mat_SMaas previous(ELE,final,ci[k],cf[k],
            fai,faf,swi,swf,s,cosines[p],input);
        const Cplx values[] = {upgraded.result(),previous.result(),
            upgraded.resoffZ2f(),upgraded.resoffZfb(),upgraded.resvpol()};
        std::cout << "  {" << final << "," << ci[k] << "," << cf[k]
                  << "," << energies[p] << "," << cosines[p] << ",{\n    ";
        for(unsigned j=0;j<5;++j)
        {
          if(!std::isfinite(real(values[j])) || !std::isfinite(imag(values[j])))
            throw std::runtime_error("nonfinite upstream reference value");
          if(j) std::cout << ",\n    ";
          std::cout << real(values[j]) << "," << imag(values[j]);
        }
        std::cout << "}},\n";
      }
  }
}
