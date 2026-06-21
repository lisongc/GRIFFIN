#include <iostream>
using namespace std;

#include "EWPOZSMEFT.h"
#include "xscnnlo.h"
#include "xscSMEFT.h"
#include "SMEFTval.h"
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

  SMval myinput1;   // convert masses from PDG values to complex pole scheme
  myinput1.set(al, 1/137.03599976);
  myinput1.set(MZ, 91.1876);
  myinput1.set(MW, 80.377);
  myinput1.set(GamZ, 2.4952);
  myinput1.set(GamW, 2.085);
  myinput1.set(MH, 125.1);
  myinput1.set(MT, 172.5);
  myinput1.set(MB, 2.87);
  myinput1.set(Delal, 0.059);
  myinput1.set(als, 0.1179);   

  SMEFTval myinput2;        
  myinput2.set(al, 1/137.03599976);
  myinput2.set(MZ, 91.1876);
  myinput2.set(MW, 80.377);
  myinput2.set(GamZ, 2.4952);
  myinput2.set(GamW, 2.085);
  myinput2.set(MH, 125.1);
  myinput2.set(MT, 172.5);
  myinput2.set(MB, 2.87);
  myinput2.set(Delal, 0.059);
  myinput2.set(als, 0.1179);
  myinput2.set(Cll,0.060516,1,1,2,2);
  myinput2.set(Cll,0.060516,1,2,2,1);
  myinput2.set(Cll,0.060516,2,2,1,1);
  myinput2.set(Cll,0.060516,2,1,1,2);
  myinput2.set(Cee,0.060516/3.,1,1,2,2);
  myinput2.set(Cee,0.060516/3.,1,2,2,1);
  myinput2.set(Cee,0.060516/3.,2,2,1,1);
  myinput2.set(Cee,0.060516/3.,2,1,1,2);
  myinput2.set(Cle,0.060516/2.,1,1,2,2);
  myinput2.set(Cle,0.060516/4.,2,2,1,1);
  myinput2.set(CphiD,0);
  myinput2.set(CphiWB,0);
  myinput2.set(Cphie,0,1,1);
  myinput2.set(Cphil1,0,1,1);
  myinput2.set(Cphil3,0,1,1);
  myinput2.set(Cphie,0,2,2);
  myinput2.set(Cphil1,0,2,2);
  myinput2.set(Cphil3,0,2,2);


  cout << endl << "Complex-pole masses: MW=" << myinput1.get(MWc) << ", MZ=" 
    << myinput1.get(MZc) << endl << endl;
  cout << endl << "Complex-pole masses: MW=" << myinput2.get(MWc) << ", MZ=" 
    << myinput2.get(MZc) << endl << endl;
 
  // compute matrix element for ee->mumu with vector coupling in initial
  // state and vector coupling in final state
  int ini = ELE, fin = MUO, iff = VEC, off = VEC;
  
  cout << "=== Matrix element for ee->ee (i=e, f=ee) ===" << endl << endl;
  
  // compute vertex form factors:
  FA_SMNNLO FAi(ini, myinput1), FAf(fin, myinput1);
  SW_SMNNLO SWi(ini, myinput1), SWf(fin, myinput1);
  FA_SMEFTLO FAi2(ini, myinput2), FAf2(fin, myinput2);
  SW_SMEFTLO SWi2(ini, myinput2), SWf2(fin, myinput2);
  cout << "F_A^i (NNLO+) = " << FAi.result() << endl;
  cout << "F_A^i (SMEFTLO+) = " << FAi2.result() << endl;
  cout << "F_A^f (NNLO+) = " << FAf.result() << endl;
  cout << "F_A^f (SMEFTLO+) = " << FAf2.result() << endl;
  cout << "sineff^i (NNLO+) = " << SWi.result() << endl;
  cout << "sineff^i (SMEFTLO+) = " << SWi2.result() << endl;
  cout << "sineff^f (NNLO+) = " << SWf.result() << endl;
  cout << "sineff^f (SMEFTLO+) = " << SWf2.result() << endl;
  cout << endl;
  
  double cme,        // center-of-mass energy
         cost = 0.5; // scattering angle
  Cplx res1, res2;

  cout << "SM matrix element M_VV for cos(theta)=" << cost << ": " << endl;
  // compute matrix element for ee->mumu using SM form factors:
  mat_SMNNLO M1(ini, fin, iff, off, FAi, FAf, SWi, SWf, cme*cme, cost, myinput1);
  cout << "sqrt(s)\t\ttot. result\t\toff-resonance contrib." << endl;
  for(cme = 10.; cme <= 190.; cme += 20.)
  {
    M1.setkinvar(cme*cme, cost);
    res1 = M1.result();
    res2 = M1.resoffZ();
    cout << cme << " \t" << res1 << " \t" << res2 << endl;
  }
  cout << endl;

  cout << "SMEFT matrix element M_VV for cos(theta)=" << cost << ": " << endl;
  // compute matrix element for ee->mumu using SMEFT form factors:
  mat_SMEFTLO M2(ini, fin, iff, off, FAi2, FAf2, SWi2, SWf2, cme*cme, cost, myinput2);
  cout << "sqrt(s)\t\ttot. result\t\toff-resonance contrib." << endl;
  for(cme = 10.; cme <= 190.; cme += 20.)
  {
    M2.setkinvar(cme*cme, cost);
    res1 = M2.result();
    res2 = M2.resoffZ();
    cout << cme << " \t" << res1 << " \t" << res2 << endl;
  }
  cout << endl;

  return 0;
}