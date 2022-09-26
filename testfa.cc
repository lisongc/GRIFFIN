#include <iostream>
using namespace std;

#include "EWPOZ2.h"
#include "SMval.h"

int main()
{
  SMval myinput;	// properly convert masses for complex pole scheme
  myinput.set(MZ, 91.1876);
  myinput.set(MW, 80.385);
  myinput.set(al, 1/137.03599976);
  myinput.set(als, 0.1184);
  myinput.set(GamZ, 2.4952);
  myinput.set(GamW, 2.085);
  myinput.set(MT, 173.2);
  myinput.set(MH, 125.1);
  myinput.set(MB, 2.87); // MSbar mass at scale mu=MZ for mb(mb)=4.20
  myinput.set(Delal, 0.059);
  myinput.set(Gmu, 1.16638e-5);
  
  cout << myinput.get(MWc) << endl;
  
  FA_SMNNLO fa2(ELE,myinput);
  cout << "aas = " << fa2.drho2aas() << " " << fa2.res2aas() << endl;
  cout << "aas2= " << fa2.drho3aas2() /*<< " " << fa2.res3aas2()*/ << endl;
  cout << "aas3= " << fa2.drho4aas3() << endl << endl;
  cout << "a2as = " << fa2.drho3a2as() << "  " <<
            "a3 = " << fa2.drho3a3() << endl;
  cout << endl;
  
  string typenam[5] = {"LEP", "NEU", "UQU", "DQU", "BQU"};
  int typenum[5] = {MUO, NUM, CQU, SQU, BQU};
  int i;
  for(i=0; i<5; i++)
  {
    fa2.setftype(typenum[i]);
    cout << typenam[i] << ":  ";
    cout << "a = " << (fa2.res1f()+fa2.res1b()) << "  " <<
    	      "a2_ff = " << fa2.res2ff() << "  a2_fb = " << fa2.res2fb() << endl;
    cout << "       a2_bb = " << fa2.res2bb() << "  aas_nf = " << fa2.res2aasnf() << "   " <<
	      "F_A = " << real(fa2.result()) << endl;
    cout << endl;
  }
  
  return 0;
}
