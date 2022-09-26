#include <iostream>
using namespace std;

#include "EWPOZ2.h"
#include "SMval.h"

int main()
{
  static double mhlist[4] = {100, 200, 600, 1000};
  int i, j;
  SMval myinput;	// properly convert masses for complex pole scheme
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.4044);
  myinput.set(al, 1/137.03599976);
  myinput.set(als, 0.119);
  myinput.set(GamZ, 2.4952);
  myinput.set(GamW, 2.115);
  myinput.set(MT, 172.5);
//  myinput.set(MB, 4.85);
  myinput.set(MB, 2.91); // MSbar mass at scale mu=MZ for MB_OS=4.85
//  myinput.set(MB, 2.87); // MSbar mass at scale mu=MZ for mb(mb)=4.20
  myinput.set(Delal, 0.05907);
  myinput.set(Gmu, 1.16637e-5);
//  myinput.set(Gmu, 1.166379e-5);
  
  double sw0 = 1-sqr(myinput.get(MWc)/myinput.get(MZc));
/*
  SW_SMNLO sw1(ELE,myinput);
  for(i=0; i<4; i++)
  {
    myinput.set(MH, mhlist[i]);
    cout << "MH = " << mhlist[i] << ":  ";
    cout << "1f = " << sw1.res1f()/sw0 << ",  ";
    cout << "1b = " << sw1.res1b()/sw0 << ",  ";
    cout << "1  = " << sw1.result()/sw0-1 << endl;
  }
  cout << endl;
*/
  SW_SMNNLO sw2(ELE,myinput);
  cout << "aas = " << sw2.drho2aas()/sw0 << " " << sw2.res2aas()/sw0 << endl;
  cout << "aas2= " << sw2.drho3aas2()/sw0 /*<< " " << sw2.res3aas2()/sw0*/ << endl;
  cout << "aas3= " << sw2.drho4aas3()/sw0 << endl << endl;
  for(i=0; i<4; i++)
  {
    myinput.set(MH, mhlist[i]);
    cout << "MH = " << mhlist[i] << ":  ";
    cout << "a2as = " << sw2.drho3a2as()/sw0 << "  " <<
            "a3 = " << sw2.drho3a3()/sw0 << endl;
  }
  
  string typenam[5] = {"LEP", "NEU", "UQU", "DQU", "BQU"};
  int typenum[5] = {MUO, NUM, CQU, SQU, BQU};
  for(j=0; j<5; j++)
  {
    if(typenum[j]==NUM) continue;
    sw2.setftype(typenum[j]);
    cout << endl << typenam[j] << ":" << endl;
    for(i=0; i<4; i++)
    {
      myinput.set(MH, mhlist[i]);
      cout << "MH = " << mhlist[i] << ":  ";
      cout << "a = " << (sw2.res1f()+sw2.res1b())/sw0 << "  " <<
    	      "a2 = " << sw2.res2ff()/sw0+sw2.res2fb()/sw0 << " " << 
		 sw2.res2bb()/sw0 << "  " <<
	      "sw_eff = " << real(sw2.result()) << endl;
    }
  }
// input values for sineff_b:
  myinput.set(MW, 80.385);
  myinput.set(GamW, 2.085);
  myinput.set(MT, 173.2);
  myinput.set(MH, 125.1);
  sw0 = 1-sqr(myinput.get(MWc)/myinput.get(MZc));
  cout << "(2016 paper) MH = 125.1:  a2_bos = " << sw2.res2bb()/sw0 << endl;
  cout << endl;
  
  return 0;
}
