#include <iostream>
using namespace std;

#include "EWPOZ2.h"
#include "SMval.h"
#include "tools.h"

int main()
{
  static double mhlist[4] = {100, 200, 600, 1000};
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
  myinput.set(MB, 2.877); // MSbar mass at scale mu=MZ for mb(mb)=4.20
  myinput.set(MC, 0.652); // MSbar mass at scale mu=MZ
  myinput.set(MS, 0);
  myinput.set(MU, 0);
  myinput.set(MD, 0);
  myinput.set(ML, 1.777);
  myinput.set(MM, 0);
  myinput.set(ME, 0);
  myinput.set(Delal, 0.05900);
  myinput.set(Gmu, 1.166379e-5);
  
  FV_SMNNLO fv2(ELE,myinput);
  FA_SMNNLO fa2(ELE,myinput);  
  // flavor selection here will be overwritten by partzwidth/zwidth below
  
  string typenam[5] = {"l", "n", "u", "d", "b"};
  int typenum[5] = {ELE, NUE, UQU, DQU, BQU};
  
  for(i = 0; i < 5; i++)
    cout << "Gamma[Z->" << typenam[i] << typenam[i] << "] = "
         << partzwidth(fa2, fv2, typenum[i], myinput, RUNWIDTHSCHEME) << endl;
  cout << "Gamma[Z] = " << zwidth(fa2, fv2, myinput, RUNWIDTHSCHEME) << endl;
  cout << endl;
  
  return 0;
}
