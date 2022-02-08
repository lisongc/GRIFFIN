#include <iostream>
using namespace std;

#include "EWPOZ.h"
#include "xscnnlo.h"
#include "SMval.h"

int main()
{
  SMval myinput;	// properly convert masses for complex pole scheme
  myinput.set(al, 1/137.03599976);
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.4044);
  myinput.set(GamZ, 2.4952);
  myinput.set(GamW, 2.115);
  myinput.set(MH, 125.1);
  myinput.set(MT, 172.5);
  myinput.set(MB, 0);
  myinput.set(Delal, 0.059);
  
  int ini = LEP, fin = BQU, iff = VEC, off = VEC;
  string typenam[5] = {"l", "n", "u", "d", "b"};
  string formnam[2] = {"V", "A"};
  
  cout << endl << "Complex-pole masses: MW=" << myinput.get(MWc) << ", MZ=" 
    << myinput.get(MZc) << endl << endl;
 
 for(fin = LEP; fin <= BQU; fin++)
 {
//  if(fin == NEU) continue;
  cout << "i = " << typenam[ini] << ", f = " << typenam[fin] <<
    " ========================================" << endl << endl;
  
  // compute vertex form factors:
  FA_SMNLO FAe(ini, myinput), FAb(fin, myinput);
  SW_SMNLO SWe(ini, myinput), SWb(fin, myinput);
  cout << "F_A^i (NLO) = " << FAe.result() << endl;
  cout << "F_A^f (NLO) = " << FAb.result() << endl;
  cout << "sineff^i (NLO) = " << SWe.result() << endl;
  cout << "sineff^f (NLO) = " << SWb.result() << endl;
  cout << endl;
  
  double cme = 88.,  // center-of-mass energy
         cost = 0.5; // scattering angle

 for(iff = VEC; iff <= AXV; iff++)
 for(off = VEC; off <= AXV; off++)
 {
  cout << formnam[iff] << formnam[off] << ":" << endl;

  // compute matrix element for ee->bb using SM form factors:
  R_SMNNLO R(ini, fin, iff, off, FAe, FAb, SWe, SWb, myinput, cost);
  S_SMNLO S(ini, fin, iff, off, myinput, cost);
  Sp_SMLO Sp(ini, fin, iff, off, myinput);
  cout << "R_" << formnam[iff] << formnam[off] << " = " << R.result() << endl;
  cout << "S_" << formnam[iff] << formnam[off] << " = " << S.result() << endl;
  cout << "S'_" << formnam[iff] << formnam[off] << "= " << Sp.result() << endl;
  
  double mz = myinput.get(MZ), gz = myinput.get(GamZ);
  Cplx sminuss0 = cme*cme - Cplx(mz*mz,-mz*gz);
  Cplx matel1 = R.result()/sminuss0 + S.result() + Sp.result()*sminuss0;
  cout << "SM matrix element = " << matel1 << endl;
  cout << " ( \" ) squared = " << sqr(abs(matel1)) << endl;
  cout << " ( \" ) squared, term by term and truncated = " 
    << sqr(abs(R.result()/sminuss0)) +
      2*real(R.result()/sminuss0*conj(S.result() + Sp.result()*sminuss0)) 
      + sqr(abs(S.result())) << endl;
  
  // directly compute squared matrix element for ee->bb using SM form factors:
  msq_SMNNLO Msq(ini, fin, iff, off, iff, off, FAe, FAb, SWe, SWb, 
  		 cme*cme, cost, myinput);
  cout << "Direct computation of SM matrix element squared = " 
  	<< Msq.result() << endl;
  cout << endl;
 }
 }
  
  return 0;
}
