/*
  File: testdeltar.cc
  Purpose: Demonstrate SM Delta-r and the Gmu W-mass calculation.
  Original SM example retained without numerical changes.
  Merged-version authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16 (typed-only input labels)
*/

#include <iostream>
using namespace std;

#include "griffin/sm/deltar.h"
#include "griffin/sm/SMvalG.h"
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

  int i;
  SMval myinput;	// properly convert masses for complex pole scheme
  myinput.set(InputPar::ZMassExperiment, 91.1876);
  myinput.set(InputPar::WMassExperiment, 80.358);
  myinput.set(InputPar::alpha, 1/137.03599976);
  myinput.set(InputPar::alphaS, 0.1179);
  myinput.set(InputPar::ZWidthExperiment, 2.4952);
  myinput.set(InputPar::WWidthExperiment, 2.089);
  myinput.set(InputPar::HiggsMass, 125.1);
  myinput.set(InputPar::topMass, 173.0);
  myinput.set(InputPar::bottomMass, 2.87);	// MSbar mass at scale mu=MZ
  myinput.set(InputPar::DeltaAlpha, 0.05900);
  myinput.set(InputPar::fermiConstant, 1.166379e-5);
  
  dr_SMNNLO dr2(myinput);
  
  // show dependence on input parameters:
  for(i=-2; i<=2; i++)
  {
    myinput.set(InputPar::topMass, 173+i);
    cout << "mt = " << 173+i << ":  ";
    cout << "Delta r = " << dr2.res1f()+dr2.res1b() << " (NLO)     ";
    cout <<                 dr2.result() << " (NNLO+)" << endl;
  }
  myinput.set(InputPar::topMass, 173.0);
  
  // illustration of Gmu input scheme, where MW is an output:
  SMvalGmu myinput2(myinput);
  cout << endl;
  cout << "complex-pole mass: mw = " << myinput2.get(InputPar::WMassComplex) << endl;
  cout << "PDG mass:          mw = " << myinput2.get(InputPar::WMassExperiment) << endl;
  cout << endl;

  return 0;
}
