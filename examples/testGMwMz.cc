#include <iostream>
using namespace std;

#include "SMvalGMwMz.h"
using namespace griffin;

int main()
{
  int i;
  SMval myinput;	// properly convert masses for complex pole scheme
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.358);
//  myinput.set(al, 1/137.03599976);
  myinput.set(als, 0.1179);
  myinput.set(GamZ, 2.4952);
  myinput.set(GamW, 2.089);
  myinput.set(MH, 125.1);
  myinput.set(MT, 173.0);
  myinput.set(MB, 2.87);	// MSbar mass at scale mu=MZ
  myinput.set(Delal, 0.05900);
  myinput.set(Gmu, 1.166379e-5);
    
  // illustration of Gmu-MW-MZ input scheme, where alpha is an output:
  SMvalGMwMz myinput2(myinput);
  cout << endl;
  cout << "complex-pole mass: mw = " << myinput2.get(MWc) << endl;
  cout << "PDG mass:          mw = " << myinput2.get(MW) << endl;
  cout << "alpha(0):          al = " << myinput2.get(al) << endl;
  cout << "                   1/al = " << 1/myinput2.get(al) << endl;
  cout << endl;

  return 0;
}
