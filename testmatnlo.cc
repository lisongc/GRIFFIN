#include <iostream>
using namespace std;

#include "EWPOZ2.h"
#include "xscenlo.h"
#include "SMval.h"

int main()
{
  SMval myinput;	// convert masses from PDG values to complex pole scheme
  myinput.set(al, 1/137.03599976);
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.385);
  myinput.set(GamZ, 2.4952);
  myinput.set(GamW, 2.085);
  myinput.set(MH, 125.0);
  myinput.set(MT, 173.0);
  myinput.set(MB, 0);
  myinput.set(Delal, 0.059);
  myinput.set(als, 0.118);
  
//  cout << endl << "Complex-pole masses: MW=" << myinput.get(MWc) << ", MZ=" 
//    << myinput.get(MZc) << endl << endl;
 
  // compute matrix element for ee->dd with vector coupling in initial
  // state and vector coupling in final state
  int ini = DQU, fin = ELE, iff = VEC, off = VEC;
  
  cout << "=== Matrix element for dd->ee (i=d, f=e) ===" << endl << endl;
  cout.precision(8);
  
  // compute vertex form factors:
  FA_SMLO FAi0(ini, myinput), FAf0(fin, myinput);
  SW_SMLO SWi0(ini, myinput), SWf0(fin, myinput);
  FA_SMNLO FAi(ini, myinput), FAf(fin, myinput);
  SW_SMNLO SWi(ini, myinput), SWf(fin, myinput);
  cout << "F_A^i (NLO) = " << FAi.result() << endl;
  cout << "F_A^f (NLO) = " << FAf.result() << endl;
  cout << "sineff^i (NLO) = " << SWi.result() << endl;
  cout << "sineff^f (NLO) = " << SWf.result() << endl;
  cout << endl;
  
  double cme,        // center-of-mass energy
         cost = 0.5; // scattering angle
  Cplx res1, res2;

  matel M0(ini, fin, iff, off, FAi0, FAf0, SWi0, SWf0, cme*cme, cost, myinput);
  mat_SMeNLO M1(ini, fin, iff, off, FAi, FAf, SWi, SWf, cme*cme, cost, myinput);

  cout << "diff. cross-section for cos(theta)=" << cost << ": " << endl;
  // compute diff. cross-section for unpolarized beams from matrix element:  
  Cplx res0vv, res0va, res0av, res0aa,
       resvv, resva, resav, resaa;
  double xsec,
         GeVtoNB = 0.38937966e6;  // unit conversion from GeV^-2 to nb
  cout << "sqrt(s)\tdsig/dcos [nb]" << endl;
  double cmel[] = {15, 40, 91, 200, 1000};
  int i;
  for(i=0; i<5; i++)
  {
    cme = cmel[i];
    M0.setkinvar(cme*cme, cost);
    M0.setform(VEC, VEC);
    res0vv = M0.result();
    M0.setform(AXV, VEC);
    res0av = M0.result();
    M0.setform(VEC, AXV);
    res0va = M0.result();
    M0.setform(AXV, AXV);
    res0aa = M0.result();
    M1.setkinvar(cme*cme, cost);
    M1.setform(VEC, VEC);
    resvv = M1.result();
    M1.setform(AXV, VEC);
    resav = M1.result();
    M1.setform(VEC, AXV);
    resva = M1.result();
    M1.setform(AXV, AXV);
    resaa = M1.result();
    xsec = real((1+cost*cost)*
                ((res0vv+2*(resvv-res0vv))*conj(res0vv) 
		 + (res0av+2*(resav-res0av))*conj(res0av)
    		 + (res0va+2*(resva-res0va))*conj(res0va) 
		 + (res0aa+2*(resaa-res0aa))*conj(res0aa)) +
		+ 4*cost*(res0vv*conj(res0aa) + (resvv-res0vv)*conj(res0aa)
		          + res0vv*conj(resaa-res0aa) 
		        + res0va*conj(res0av) + (resva-res0va)*conj(res0av)
			  + res0va*conj(resav-res0av)));
    xsec *= cme*cme*cme*cme/3;
    cout << cme << " \t" << xsec << endl;
  }
  
  return 0;
}
