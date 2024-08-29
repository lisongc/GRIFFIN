#include <iostream>
using namespace std;

#include "EWPOZ2.h"
#include "xscaas.h"
#include "SMval.h"
using namespace griffin;

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

  SMval myinput;        // convert masses from PDG values to complex pole scheme
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

  // compute matrix element for ee->dd with vector coupling in initial
  // state and vector coupling in final state
  int ini = ELE, fin = DQU, iff = VEC, off = VEC;
  
  cout << "=== Matrix element for ee->dd (i=e, f=d) ===" << endl << endl;

  // compute vertex form factors:
  FA_SMNNLO FAi(ini, myinput), FAf(fin, myinput);
  SW_SMNNLO SWi(ini, myinput), SWf(fin, myinput);

  double cme,        // center-of-mass energy
         cost = 0.5; // scattering angle
  Cplx res1, res2;

  cout << "SM matrix element M_VV for cos(theta)=" << cost << ": " << endl;
  // compute matrix element for ee->dd using SM form factors:
  mat_SMaas M(ini, fin, iff, off, FAi, FAf, SWi, SWf, cme*cme, cost, myinput);
  cout << "sqrt(s)\ttot. result without O(aas)\ttot. result with O(aas)" << endl;
  for(cme = 10.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    res1 = M.mat_SMNNLO::result();
    res2 = M.result();
    cout << cme << " \t" << res1 << " \t" << res2 << endl;
  }
  cout << endl;
  cout << "Without O(aas) in matrix el.: R=" << M.mat_SMNNLO::coeffR() << endl;
  cout << "With O(aas) in matrix el.: R=" << M.coeffR() << endl << endl;
  cout << "sqrt(s)\toff-resonance without O(aas)\toff-resonance with O(aas)" << endl;
  for(cme = 10.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    res1 = M.resoffZ();
    res2 = M.resoffZaas();
    cout << cme << " \t" << res1-res2 << " \t" << res1 << endl;
  }
  cout << endl;

  cout << "diff. cross-section for cos(theta)=" << cost << ": " << endl;
  // compute diff. cross-section for unpolarized beams from matrix element:  
  Cplx res1vv, res1va, res1av, res1aa;
  Cplx res2vv, res2va, res2av, res2aa;
  double xsec,
         GeVtoNB = 0.38937966e6;  // unit conversion from GeV^-2 to nb
  cout << "sqrt(s)\tdsig/dcos [nb] without O(aas)\tdsig/dcos [nb] with O(aas)" << endl;
  for(cme = 10.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    M.setform(VEC, VEC);
    res1vv = M.mat_SMNNLO::result();
    res2vv = M.result();
    M.setform(AXV, VEC);
    res1av = M.mat_SMNNLO::result();
    res2av = M.result();
    M.setform(VEC, AXV);
    res1va = M.mat_SMNNLO::result();
    res2va = M.result();
    M.setform(AXV, AXV);
    res1aa = M.mat_SMNNLO::result();
    res2aa = M.result();
    xsec = real((1+cost*cost)*(res1vv*conj(res1vv) + res1av*conj(res1av)
    			     + res1va*conj(res1va) + res1aa*conj(res1aa)) +
		+ 4*cost*(res1vv*conj(res1aa) + res1va*conj(res1av)));
    xsec *= 3*cme*cme/(32*Pi) * GeVtoNB;
    cout << cme << " \t" << xsec << " \t\t\t";
    xsec = real((1+cost*cost)*(res2vv*conj(res2vv) + res2av*conj(res2av)
    			     + res2va*conj(res2va) + res2aa*conj(res2aa)) +
		+ 4*cost*(res2vv*conj(res2aa) + res2va*conj(res2av)));
    xsec *= 3*cme*cme/(32*Pi) * GeVtoNB;
    cout << xsec << endl;
  }

  return 0;
}
