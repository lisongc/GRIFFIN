#include <iostream>
using namespace std;

#include "EWPOZ.h"
#include "xscnnlo.h"
#include "SMval.h"

int main()
{
  SMval myinput;	// convert masses from PDG values to complex pole scheme
  myinput.set(al, 1/137.03599976);
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.4044);
  myinput.set(GamZ, 2.4952);
  myinput.set(GamW, 2.115);
  myinput.set(MH, 125.1);
  myinput.set(MT, 172.5);
  myinput.set(MB, 0);
  myinput.set(Delal, 0.059);
  
  cout << endl << "Complex-pole masses: MW=" << myinput.get(MWc) << ", MZ=" 
    << myinput.get(MZc) << endl << endl;
 
  // compute matrix element for ee->bb with axial-vector coupling in initial
  // state and vector coupling in final state
  int ini = LEP, fin = BQU, iff = AXV, off = VEC;
  string typenam[5] = {"l", "n", "u", "d", "b"};
  string formnam[2] = {"V", "A"};
  
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

  // compute matrix element for ee->ff using SM form factors:
  mat_SMNNLO M(ini, fin, iff, off, FAe, FAb, SWe, SWb, cme*cme, cost, myinput);
  cout << "R_" << formnam[iff] << formnam[off] << " = " << M.coeffR() << endl;
  cout << "S_" << formnam[iff] << formnam[off] << " = " << M.coeffS() << endl;
  cout << "S'_" << formnam[iff] << formnam[off] << "= " << M.coeffSp() << endl;
  
  Cplx matel1 = M.result();
  cout << "SM matrix element M_" << formnam[iff] << formnam[off] << " = " << matel1 << endl;
  cout << " ( \" ) squared |M_" << formnam[iff] << formnam[off] << "|^2 = " << sqr(abs(matel1)) << endl;
  
  // as above, but now with user-supplied SM form factors, 
  // as would be used, e.g., for fitting:
  double FAue = 0.0345, FAub = 0.0340, SWue = 0.231, SWub = 0.233;
  cout << "F_A^i (user) = " << FAue << endl;
  cout << "F_A^f (user) = " << FAub << endl;
  cout << "sineff^i (user) = " << SWue << endl;
  cout << "sineff^f (user) = " << SWub << endl;
  cout << endl;
  
  // compute matrix element for ee->ff using user form factors:
  mat_SMNNLO Mu(ini, fin, iff, off, FAue, FAub, SWue, SWub, cme*cme, cost, myinput);
  cout << "R_" << formnam[iff] << formnam[off] << " = " << Mu.coeffR() << endl;
  cout << "S_" << formnam[iff] << formnam[off] << " = " << Mu.coeffS() << endl;
  cout << "S'_" << formnam[iff] << formnam[off] << "= " << Mu.coeffSp() << endl;
  
  matel1 = Mu.result();
  cout << "Matrix element with user form factors M_" << formnam[iff] << formnam[off] << " = " << matel1 << endl;
  cout << " ( \" ) squared |M_" << formnam[iff] << formnam[off] << "|^2 = " << sqr(abs(matel1)) << endl;
  
  return 0;
}
