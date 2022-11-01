#include <iostream>
using namespace std;

#include "deltar.h"
#include "SMvalG.h"

int main()
{
  static double mhlist[5] = {100, 200, 300, 600, 1000};
  int i;
  SMval myinput;	// properly convert masses for complex pole scheme
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.358);
  myinput.set(al, 1/137.03599976);
  myinput.set(als, 0.1179);
  myinput.set(GamZ, 2.4952);
  myinput.set(GamW, 2.089);
  myinput.set(MT, 173.0);
//  myinput.set(MB, 4.7);
  myinput.set(MB, 2.87);	// MSbar mass at scale mu=MZ
  myinput.set(Delal, 0.05900);
  myinput.set(Gmu, 1.166379e-5);
  
  dr_SMNLO dr1(myinput);
/*  for(i=0; i<5; i++)
  {
    myinput.set(MH, mhlist[i]);
    cout << "MH = " << mhlist[i] << ":  ";
    cout << "1f = " << dr1.res1f() << ",  ";
    cout << "1b = " << dr1.res1b() << ",  ";
    cout << "1  = " << dr1.result() << endl;
  }
  cout << endl;
*/
  dr_SMNNLO dr2(myinput);
  for(i=0; i<5; i++)
  {
    myinput.set(MH, mhlist[i]);
    cout << "MH = " << mhlist[i] << ":  ";
    cout << "Delta r = " << dr2.result() << endl;
  }
  cout << endl;
  cout << "      from DelRho    all" << endl;
  cout << "aas = " << dr2.drho2aas() << "  " << dr2.res2aas() << endl;
  cout << "aas2= " << dr2.drho3aas2() << "  " << dr2.res3aas2() << endl;
  cout << "aas3= " << dr2.drho4aas3() << endl;
  cout << "a3(Nf=3)=          " << dr2.res3fff() << endl;
  cout << "a2as(Nf=2)=        " << dr2.res3ffa2as() << endl << endl;
  for(i=0; i<5; i++)
  {
    myinput.set(MH, mhlist[i]);
    cout << "MH = " << mhlist[i] << ":  ";
    cout << "a = " << dr2.res1f()+dr2.res1b() << "  " <<
    	    "a2as = " << dr2.drho3a2as() << "  " <<
            "a3 = " << dr2.drho3a3() << endl;
  }
  cout << endl;
  cout << "                  Nf=2         Nf=1        Nf=1+2      Nf=0" << endl;
  for(i=0; i<5; i++)
  {
    myinput.set(MH, mhlist[i]);
    cout << "MH = " << mhlist[i] << ":  ";
    cout << "a2 = " << dr2.res2ff() << "  " << dr2.res2fb() << "  " << 
        dr2.res2ff()+dr2.res2fb() << "  " << 
	dr2.res2bb() << endl;
  }
  cout << endl;
  
  SMvalGmu myinput2(myinput);
  cout << "                comp.-pole mass  PDG mass" << endl;
  for(i=0; i<5; i++)
  {
    myinput2.set(MH, mhlist[i]);
    cout << "MH = " << mhlist[i] << ":  ";
    cout << "mw = " << myinput2.get(MWc) << "          " << myinput2.get(MW) << endl;
  }
  cout << endl;
/*  
  invalGmu inp2;
  inp2.set(al, 1/137.03599976);
  inp2.set(MZc, 91.1534);
//  inp2.set(MWc, sqrt(6431.85706042165));
  inp2.set(MH, 100);
  inp2.set(MT, 174.3);
  inp2.set(als, 0.119);
  inp2.set(MB, 2.87);
  inp2.set(Delal, 0.05907);
  inp2.set(Gmu, 1.166379e-5);
  dr_SMNNLO dr3(inp2);
  cout << inp2.get(MWc) << " " << dr3.result() << endl;
*/  
  return 0;
}
