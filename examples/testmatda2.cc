// Demonstrate use of NNLO cross-section calculation with (\Delta\alpha)^2 terms

#include <iostream>
using namespace std;

#include "EWPOZ2.h"
#include "xscenlo.h"
#include "xscda2.h"
#include "SMval.h"
using namespace griffin;

int main()
{
  SMval myinput;	// convert masses from PDG values to complex pole scheme
  myinput.set(al, 1/137.03599976);  // for scheme 0
  myinput.set(Delal, 0.059);
//  myinput.set(al, 1/137.03599976/(1-0.059));  // for scheme 1
//  myinput.set(Delal, 0);
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.385);
  myinput.set(GamZ, 2.4952);
  myinput.set(GamW, 2.085);
  myinput.set(MH, 125.0);
  myinput.set(MT, 173.0);
  myinput.set(MB, 0);
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
  FA_SMNLO FAi1(ini, myinput), FAf1(fin, myinput);
  SW_SMNLO SWi1(ini, myinput), SWf1(fin, myinput);
  FA_SMNNLO FAi2(ini, myinput), FAf2(fin, myinput);
  SW_SMNNLO SWi2(ini, myinput), SWf2(fin, myinput);
  cout << "F_A^i (NLO) = " << FAi1.result() << endl;
  cout << "F_A^f (NLO) = " << FAf1.result() << endl;
  cout << "sineff^i (NLO) = " << SWi1.result() << endl;
  cout << "sineff^f (NLO) = " << SWf1.result() << endl;
  cout << endl;
  cout << "F_A^i (NNLO) = " << FAi2.result() << endl;
  cout << "F_A^f (NNLO) = " << FAf2.result() << endl;
  cout << "sineff^i (NNLO) = " << SWi2.result() << endl;
  cout << "sineff^f (NNLO) = " << SWf2.result() << endl;
  cout << endl;
  
  double cme,        // center-of-mass energy
         cost = 0.5; // scattering angle
  Cplx res1, res2;

  matel M0(ini, fin, iff, off, FAi0, FAf0, SWi0, SWf0, cme*cme, cost, myinput);
  mat_SMeNLO M1(ini, fin, iff, off, FAi1, FAf1, SWi1, SWf1, cme*cme, cost, myinput);
  mat_SMNNLO M2(ini, fin, iff, off, FAi2, FAf2, SWi2, SWf2, cme*cme, cost, myinput);
  mat_SMda2 Mda2(ini, fin, iff, off, FAi2, FAf2, SWi2, SWf2, cme*cme, cost, myinput);

  cout << "diff. cross-section for cos(theta)=" << cost << ": " << endl;
  // compute diff. cross-section for unpolarized beams from matrix element:  
  Cplx res0vv, res0va, res0av, res0aa,
       resvv, resva, resav, resaa;
  double xsec,
         GeVtoNB = 0.38937966e6;  // unit conversion from GeV^-2 to nb
  double cmel[] = {15, 40, 91, 200, 1000};
  int i;
  cout << "sqrt(s)\t|M^2| (NLO)\t|M^2| (NNLO)\t|M^2| (NNLO, lin.exp.)" << endl;
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
    cout << cme << " \t" << xsec;
    M2.setkinvar(cme*cme, cost);
    M2.setform(VEC, VEC);
    resvv = M2.result();
    M2.setform(AXV, VEC);
    resav = M2.result();
    M2.setform(VEC, AXV);
    resva = M2.result();
    M2.setform(AXV, AXV);
    resaa = M2.result();
    xsec = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
                             + resva*conj(resva) + resaa*conj(resaa)) +
                + 4*cost*(resvv*conj(resaa) + resva*conj(resav)));
    xsec *= cme*cme*cme*cme/3;
    cout << " \t" << xsec;
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
    cout << " \t" << xsec << endl;
  }
  
  cout << endl << "sqrt(s)\t|M^2| (NNLO w.da^2)\t|M^2| (NNLO w.da^2, lin.exp.)" << endl;
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
    Mda2.setkinvar(cme*cme, cost);
    Mda2.setform(VEC, VEC);
    resvv = Mda2.result();
    Mda2.setform(AXV, VEC);
    resav = Mda2.result();
    Mda2.setform(VEC, AXV);
    resva = Mda2.result();
    Mda2.setform(AXV, AXV);
    resaa = Mda2.result();
    xsec = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
                             + resva*conj(resva) + resaa*conj(resaa)) +
                + 4*cost*(resvv*conj(resaa) + resva*conj(resav)));
    xsec *= cme*cme*cme*cme/3;
    cout << cme << " \t" << xsec;
// add one \Delta\alpha^2 term by hand, which would normally come from |M_1|^2:
    double gi = g0(ini,VEC,myinput), gf = g0(fin,VEC,myinput);
    double resvvda2 = gi*gf*sqr(myinput.get(Delal))/sqr(cme);
    xsec = real((1+cost*cost)*
                ((res0vv+2*(resvv-res0vv) + resvvda2)*conj(res0vv) 
		 + (res0av+2*(resav-res0av))*conj(res0av)
    		 + (res0va+2*(resva-res0va))*conj(res0va) 
		 + (res0aa+2*(resaa-res0aa))*conj(res0aa)) +
		+ 4*cost*(res0vv*conj(res0aa) + (resvv-res0vv)*conj(res0aa)
		          + res0vv*conj(resaa-res0aa) 
		        + res0va*conj(res0av) + (resva-res0va)*conj(res0av)
			  + res0va*conj(resav-res0av)));
    xsec *= cme*cme*cme*cme/3;
    cout << " \t\t" << xsec << endl;
  }
  
  return 0;
}
