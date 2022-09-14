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
  int ini = LEP, fin = DQU, iff = VEC, off = VEC;
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
  
  double cme = 10.,  // center-of-mass energy
         cost = 0.5; // scattering angle
  Cplx res1, res2, res20;

  cout << "SM matrix element M_" << formnam[iff] << formnam[off] << ": " << endl;
  // compute matrix element for ee->ff using SM form factors:
  mat_SMNNLO M(ini, fin, iff, off, FAe, FAb, SWe, SWb, cme*cme, cost, myinput);
  matel M0(ini, fin, iff, off, FAe, FAb, SWe, SWb, cme*cme, cost, myinput);

  for(cme = 10.; cme <= 200.; cme += 5.)
  {
    M.setkinvar(cme*cme, cost);
    M0.setkinvar(cme*cme, cost);
    res1 = M.result();
    res20 = M0.resoffZ();
    res2 = M.resoffZ();
    cout << "sqrt(s) = " << cme << " " << res1-res2 << " " << res20 << " " << res2 <<
    	endl;
  }

  return 0;
}
