#include <iostream>
using namespace std;

#include "EWPOZSMEFT.h"
#include "xscaas.h"
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
  myinput2.set(CphiD,0.060516);
  myinput2.set(CphiWB,0.060516);
  myinput2.set(Cphil1,0.060516,1,1);
  myinput2.set(Cphil3,0.060516,1,1);
  myinput2.set(Cphiq3,0.060516,1,1);
  myinput2.set(Cphie,0.060516,1,1);


  // compute effective weak mixing angle and vertex form factors for ee with vector coupling in initial
  // state
  int ini = ELE;
  
  // compute vertex form factors:
  SW_SMNNLO SW1(ini, myinput1);
  SW_SMEFTLO SW2(ini, myinput2);

  FA_SMNNLO FA1(ini, myinput1);
  FA_SMEFTLO FA2(ini, myinput2);

  FV_SMNNLO FV1(ini, myinput1);
  FV_SMEFTLO FV2(ini, myinput2);

  cout << myinput1.get(MWc) << endl;
  cout << myinput1.get(MZc) << endl;

  cout << myinput2.get(MWc) << endl;
  cout << myinput2.get(MZc) << endl;

  cout << "NNLO SM SW for initial e+e-: " << SW1.result() << endl;
  cout << "SMEFT SW for initial e+e-: " << SW2.result() << endl;

  cout << "NNLO SM FA for initial e+e-: " << FA1.result() << endl;
  cout << "SMEFT FA for initial e+e-: " << FA2.result() << endl;

  cout << "NNLO SM FV for initial e+e-: " << FV1.result() << endl;
  cout << "SMEFT FV for initial e+e-: " << FV2.result() << endl;

  return 0;
}