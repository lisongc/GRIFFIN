#include <iostream>
using namespace std;

#include "deltar.h"
#include "SMvalG.h"

int main()
{
  cout<<"======================================================"<<endl;
  cout<<"                                                      "<<endl;
  cout<<"======================================================"<<endl;

  cout<<"     ______ ____   ____ ______ ______ ____ _   __"<<endl;
  cout<<"    / ____// __ \\ /  _// ____// ____//  _// | / /"<<endl;
  cout<<"   / / __ / /_/ / / / / /_   / /_    / / /  |/ / "<<endl;
  cout<<"  / /_/ // _, _/_/ / / __/  / __/  _/ / / /|  /  "<<endl;
  cout<<"  \\____//_/ |_|/___//_/    /_/    /___//_/ |_/   "<<endl;

  cout<<"======================================================"<<endl;
  cout<<"                 version 1.0                          "<<endl;
  cout<<"======================================================"<<endl;
  int i;
  SMval myinput;	// properly convert masses for complex pole scheme
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.358);
  myinput.set(al, 1/137.03599976);
  myinput.set(als, 0.1179);
  myinput.set(GamZ, 2.4952);
  myinput.set(GamW, 2.089);
  myinput.set(MH, 125.1);
  myinput.set(MT, 173.0);
  myinput.set(MB, 2.87);	// MSbar mass at scale mu=MZ
  myinput.set(Delal, 0.05900);
  myinput.set(Gmu, 1.166379e-5);
  
  dr_SMNNLO dr2(myinput);
  
  // show dependence on input parameters:
  for(i=-2; i<=2; i++)
  {
    myinput.set(MT, 173+i);
    cout << "mt = " << 173+i << ":  ";
    cout << "Delta r = " << dr2.res1f()+dr2.res1b() << " (NLO)     ";
    cout <<                 dr2.result() << " (NNLO+)" << endl;
  }
  myinput.set(MT, 173.0);
  
  // illustration of Gmu input scheme, where MW is an output:
  SMvalGmu myinput2(myinput);
  cout << endl;
  cout << "complex-pole mass: mw = " << myinput2.get(MWc) << endl;
  cout << "PDG mass:          mw = " << myinput2.get(MW) << endl;
  cout << endl;

  return 0;
}
