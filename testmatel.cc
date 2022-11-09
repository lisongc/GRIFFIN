#include <iostream>
using namespace std;

#include "EWPOZ2.h"
#include "xscnnlo.h"
#include "SMval.h"

int main()
{
  cout<<"======================================================"<<endl;
  cout<<"======================================================"<<endl;

  cout<<"     ______ ____   ____ ______ ______ ____ _   __"<<endl;
  cout<<"    / ____// __ \\ /  _// ____// ____//  _// | / /"<<endl;
  cout<<"   / / __ / /_/ / / / / /_   / /_    / / /  |/ / "<<endl;
  cout<<"  / /_/ // _, _/_/ / / __/  / __/  _/ / / /|  /  "<<endl;
  cout<<"  \\____//_/ |_|/___//_/    /_/    /___//_/ |_/   "<<endl;

  cout<<"======================================================"<<endl;
  cout<<"======================================================"<<endl;

  SMval myinput;	// convert masses from PDG values to complex pole scheme
  myinput.set(al, 1/137.03599976);
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.377);
  myinput.set(GamZ, 2.4952);
  myinput.set(GamW, 2.085);
  myinput.set(MH, 125.1);
  myinput.set(MT, 172.5);
  myinput.set(MB, 2.87);
  myinput.set(Delal, 0.059);
  myinput.set(als, 0.1179);
  
  cout << endl << "Complex-pole masses: MW=" << myinput.get(MWc) << ", MZ=" 
    << myinput.get(MZc) << endl << endl;
 
  // compute matrix element for ee->dd with vector coupling in initial
  // state and vector coupling in final state
  int ini = ELE, fin = DQU, iff = VEC, off = VEC;
  
  cout << "=== Matrix element for ee->dd (i=e, f=d) ===" << endl << endl;
  
  // compute vertex form factors:
  FA_SMNNLO FAi(ini, myinput), FAf(fin, myinput);
  SW_SMNNLO SWi(ini, myinput), SWf(fin, myinput);
  cout << "F_A^i (NNLO+) = " << FAi.result() << endl;
  cout << "F_A^f (NNLO+) = " << FAf.result() << endl;
  cout << "sineff^i (NNLO+) = " << SWi.result() << endl;
  cout << "sineff^f (NNLO+) = " << SWf.result() << endl;
  cout << endl;
  
  double cme,        // center-of-mass energy
         cost = 0.5; // scattering angle
  Cplx res1, res2;

  cout << "SM matrix element M_VV for cos(theta)=" << cost << ": " << endl;
  // compute matrix element for ee->dd using SM form factors:
  mat_SMNNLO M(ini, fin, iff, off, FAi, FAf, SWi, SWf, cme*cme, cost, myinput);
  cout << "sqrt(s)\t\ttot. result\t\toff-resonance contrib." << endl;
  for(cme = 10.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    res1 = M.result();
    res2 = M.resoffZ();
    cout << cme << " \t" << res1 << " \t" << res2 << endl;
  }
  cout << endl;

  cout << "diff. cross-section for cos(theta)=" << cost << ": " << endl;
  // compute diff. cross-section for unpolarized beams from matrix element:  
  Cplx resvv, resva, resav, resaa;
  double xsec,
         GeVtoNB = 0.38937966e6;  // unit conversion from GeV^-2 to nb
  cout << "sqrt(s)\tdsig/dcos [nb]" << endl;
  for(cme = 10.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    M.setform(VEC, VEC);
    resvv = M.result();
    M.setform(AXV, VEC);
    resav = M.result();
    M.setform(VEC, AXV);
    resva = M.result();
    M.setform(AXV, AXV);
    resaa = M.result();
    xsec = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
    			     + resva*conj(resva) + resaa*conj(resaa)) +
		+ 4*cost*(resvv*conj(resaa) + resva*conj(resav)));
    xsec *= 3*cme*cme/(32*Pi) * GeVtoNB;
    cout << cme << " \t" << xsec << endl;
  }
  
  return 0;
}
