/* test file for computing differential cross
section with respect to cost*/

#include <iostream>
using namespace std;

#include "EWPOZ.h"
#include "EWPOZ2.h"
#include "xscnnlo.h"
#include "SMval.h"

int main()
{
  SMval myinput; // convert masses from PDG values to complex pole scheme
  myinput.set(al, 1 / 137.03599976);
  myinput.set(als, 0.118);
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.35321);
  myinput.set(GamZ, 2.494684);
  myinput.set(GamW, 2.089792);
  myinput.set(MH, 125.0);
  myinput.set(MT, 173.1);
  myinput.set(MB, 0);
  myinput.set(Delal, 0.059);

  // process ee -> ll
  int ini = ELE, fin = MUE;


  double cme = 86.1876, // center-of-mass energy
      cost = 0.0,       // scattering angle
      pe = 0.,          // polarization degree of the incoming electron beam
      Nc = 1;           // the color factor of the outgoing fermion
  double mz = myinput.get(MZc), gz = myinput.get(GZc);

  // compute SM form factors and Laurent coeffs
  FA_SMNNLO FAe(ini, myinput), FAb(fin, myinput);
  SW_SMNNLO SWe(ini, myinput), SWb(fin, myinput);

  Cplx amp[2][2];
  string formlist[2] = {"VEC", "AXV"};

  mat_SMNNLO Mij(ini, fin, VEC, VEC, FAe, FAb, SWe, SWb, cme*cme, cost, myinput);
  
  for (int i = 1; i < 21; i++)
  {
    cme += 0.5;
    Cplx sminuss0 = cme * cme - Cplx(mz * mz, -mz * gz);
    for (int m = VEC; m <= AXV; m++)
    {
      for (int n = VEC; n <= AXV; n++)
      {
        Mij.setkinvar(cme*cme, cost);
        Mij.setform(m, n);
        amp[m][n] = Mij.result();
      }
    }
   
    /* compute the differential cross-section by R, S, and S'*/ 
     cout << cme << " " << real(((Nc * cme * cme) / (32 * Pi)) * ((1 + cost * cost) * (sqr(abs(amp[0][0])) + sqr(abs(amp[1][1])) + sqr(abs(amp[0][1])) + sqr(abs(amp[1][0]))) + 4 * cost * (amp[0][0] * conj(amp[1][1]) + amp[0][1] * conj(amp[1][0])) - 2 * pe * (1 + cost * cost) * (amp[0][0] * conj(amp[1][0]) + amp[0][1] * conj(amp[1][1])) - 4 * pe * cost * (amp[0][0] * conj(amp[0][1]) + amp[1][0] * conj(amp[1][1])))) << endl;
  }
  return 0;
}

