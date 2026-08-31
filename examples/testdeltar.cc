#include <iostream>
#include <iomanip>
using namespace std;

#include "deltar.h"
#include "EWPOZSMEFT.h"
#include "SMEFTvalG.h"
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

  SMval myinput1;	// properly convert masses for complex pole scheme
  myinput1.set(MZ, 91.1876);
  myinput1.set(MW, 80.377);
  myinput1.set(al, 1/137.03599976);
  myinput1.set(als, 0.1179);
  myinput1.set(GamZ, 2.4952);
  myinput1.set(GamW, 2.085);
  myinput1.set(MH, 125.1);
  myinput1.set(MT, 172.5);
  myinput1.set(MB, 2.87);	// MSbar mass at scale mu=MZ
  myinput1.set(Delal, 0.05900);
  myinput1.set(Gmu, 1.166379e-5);
  
  // illustration of Gmu input scheme, where MW is an output:
  SMvalGmu myinput2(myinput1);
  cout << endl;
  cout << "complex-pole mass: mw = " << myinput2.get(MWc) << endl;
  cout << "PDG mass:          mw = " << myinput2.get(MW) << endl;
  cout << endl;

  SMEFTvalGmu myinput3(myinput1);
  myinput3.set(Cll,0,1,2,2,1);
  myinput3.set(Cll,0,2,1,1,2);
  myinput3.set(Cphil3,0.060516*0.001,1,1);
  myinput3.set(Cphil3,0,2,2);
  myinput3.set(CphiD,0);
  myinput3.set(CphiWB,0);
  myinput3.set(Cphie,0,1,1);
  myinput3.set(Cphil1,0,1,1);
  cout << endl;
  cout << "complex-pole mass: mw = " << myinput3.get(MWc) << endl;
  cout << "PDG mass:          mw = " << myinput3.get(MW) << endl;
  cout << endl;

  SMEFTval myinput4(myinput1);
  myinput4.set(CphiD,0);
  myinput4.set(CphiWB,0.001);
  myinput4.set(Cphie,0,1,1);
  myinput4.set(Cphil1,0,1,1);
  myinput4.set(Cphil3,0,1,1);

  cout << endl;
  cout << "complex-pole mass: mw = " << myinput4.get(MWc) << endl;
  cout << "PDG mass:          mw = " << myinput4.get(MW) << endl;
  cout << endl;

  int ini = ELE;

  SW_SMNNLO SW1(ini, myinput1);
  SW_SMNNLO SW2(ini, myinput2);
  SW_SMEFTLO SW3(ini, myinput3);
  SW_SMEFTLO SW4(ini, myinput4);

  FA_SMNNLO FA1(ini, myinput1);
  FA_SMNNLO FA2(ini, myinput2);
  FA_SMEFTLO FA3(ini, myinput3);
  FA_SMEFTLO FA4(ini, myinput4);

  FV_SMNNLO FV1(ini, myinput1);
  FV_SMNNLO FV2(ini, myinput2);
  FV_SMEFTLO FV3(ini, myinput3);
  FV_SMEFTLO FV4(ini, myinput4);

  Cplx resultSW1 = SW3.result();
  Cplx resultFA1 = FA3.result();
  Cplx resultFV1 = FV3.result();
  myinput3.set(Cphil3,-0.001*0.060516);

  cout << std::setprecision(9) << "Jackson: " << (resultSW1-SW3.result())/(0.002*0.060516) << endl;
  cout << "Jackson: " << (resultFA1-FA3.result())/(0.002*0.060516) << endl;
  cout << "Jackson: " << (resultFV1-FV3.result())/(0.002*0.060516) << endl;


  /*cout << std::setprecision (15) <<  "NNLO SM SW (alpha) for initial e+e-: " << SW1.result() << endl;
  cout << "NNLO SM SW (Gmu) for initial e+e-: " << SW2.result() << endl;
  cout << "SMEFT SW (alpha) for initial e+e-: " << SW4.result() << endl;
  cout << "SMEFT SW (Gmu) for initial e+e-: " << SW3.result() << endl;

  cout << "NNLO SM FA (alpha) for initial e+e-: " << FA1.result() << endl;
  cout << "NNLO SM FA (Gmu) for initial e+e-: " << FA2.result() << endl;
  cout << "SMEFT FA (alpha) for initial e+e-: " << FA4.result() << endl;
  cout << "SMEFT FA (Gmu) for initial e+e-: " << FA3.result() << endl;

  cout << "NNLO SM FV (alpha) for initial e+e-: " << FV1.result() << endl;
  cout << "NNLO SM FV (Gmu) for initial e+e-: " << FV2.result() << endl;
  cout << "SMEFT FV (alpha) for initial e+e-: " << FV4.result() << endl;
  cout << "SMEFT FV (Gmu) for initial e+e-: " << FV3.result() << endl;*/

  return 0;
}
