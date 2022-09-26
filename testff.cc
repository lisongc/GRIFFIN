#include <iostream>
#include <iomanip>
using namespace std;

#include "ff.h"

//#include "clooptools.h"

int main()
{
  inval myinput;
/*  myinput.set(al, 1/137.037);
  myinput.set(MZ, 91.2);
  myinput.set(MW, 80.4);
  myinput.set(MH, 125.1);
  myinput.set(MT, 173.2);
  myinput.set(Delal, 0.059);
*/  
// Lisong's numbers:
  myinput.set(al, 0.303*0.303/(4*Pi));
  myinput.set(MZ, 91.0);
  myinput.set(GamZ, 2.5);
  myinput.set(MW, 80.0);
  myinput.set(MH, 125.0);
  myinput.set(MT, 173.0);
  myinput.set(MB, 0);
  myinput.set(Delal, 0.059);
  
  cout << "Re SigZ_f = " << rsz1f(myinput) << endl;
  cout << "Im SigZ_f = " << isz1f(myinput) << endl;
  cout << "Re SigZ_b = " << rsz1b(myinput) << endl;
  cout << "Im SigZ_b = " << isz1b(myinput) << endl;
  cout << "Re SigZ'_f = " << rsz1fp(myinput) << endl;
  cout << "Im SigZ'_f = " << isz1fp(myinput) << endl;
  cout << "Re SigZ'_b = " << rsz1bp(myinput) << endl;
  cout << "Im SigZ'_b = " << isz1bp(myinput) << endl;
  cout << "Re SigZ\"_f = " << rsz1fpp(myinput) << endl;
  cout << "Im SigZ\"_f = " << isz1fpp(myinput) << endl;
  cout << "Re SigZ\"_b = " << rsz1bpp(myinput) << endl;
  cout << "Im SigZ\"_b = " << isz1bpp(myinput) << endl;
  cout << endl;
  
  string typenam[5] = {"l", "n", "u", "d", "b"};
  int typenum[5] = {MUO, NUM, CQU, SQU, BQU};
  string formnam[2] = {"V", "A"};
  int i,j,k;
  cout << "      ";
  for(i=0; i<5; i++)
    cout << setw(13) << typenam[i];
  cout << endl;
  for(j=0; j<2; j++)
  {
    cout << "Re Z" << formnam[j] << "_f = ";
    for(i=0; i<5; i++)
      cout << setw(13) << rz1f(typenum[i],j,myinput);
    cout << endl;
    cout << "Im Z" << formnam[j] << "_f = ";
    for(i=0; i<5; i++)
      cout << setw(13) << iz1f(typenum[i],j,myinput);
    cout << endl;
    cout << "Re Z" << formnam[j] << "_b = ";
    for(i=0; i<5; i++)
      cout << setw(13) << rz1b(typenum[i],j,myinput);
    cout << endl;
    cout << "Im Z" << formnam[j] << "_b = ";
    for(i=0; i<5; i++)
      cout << setw(13) << iz1b(typenum[i],j,myinput);
    cout << endl;
  }
  cout << endl;
  for(j=0; j<2; j++)
  {
    cout << "Re Z'" << formnam[j] << "_f= ";
    for(i=0; i<5; i++)
      cout << setw(13) << rz1fp(typenum[i],j,myinput);
    cout << endl;
    cout << "Im Z'" << formnam[j] << "_f= ";
    for(i=0; i<5; i++)
      cout << setw(13) << iz1fp(typenum[i],j,myinput);
    cout << endl;
    cout << "Re Z'" << formnam[j] << "_b= ";
    for(i=0; i<5; i++)
      cout << setw(13) << rz1bp(typenum[i],j,myinput);
    cout << endl;
    cout << "Im Z'" << formnam[j] << "_b= ";
    for(i=0; i<5; i++)
      cout << setw(13) << iz1bp(typenum[i],j,myinput);
    cout << endl;
  }
  cout << endl;
  for(j=0; j<2; j++)
  {
    cout << "Re G" << formnam[j] << "_f = ";
    for(i=0; i<5; i++)
      cout << setw(13) << rg1f(typenum[i],j,myinput);
    cout << endl;
    cout << "Im G" << formnam[j] << "_f = ";
    for(i=0; i<5; i++)
      cout << setw(13) << ig1f(typenum[i],j,myinput);
    cout << endl;
    cout << "Re G" << formnam[j] << "_b = ";
    for(i=0; i<5; i++)
      cout << setw(13) << rg1b(typenum[i],j,myinput);
    cout << endl;
    cout << "Im G" << formnam[j] << "_b = ";
    for(i=0; i<5; i++)
      cout << setw(13) << ig1b(typenum[i],j,myinput);
    cout << endl;
  }
  cout << endl;
/*
  cout << "Re SigZ = " << rsz1f(myinput)+rsz1b(myinput) << endl;
  cout << "Im SigZ = " << isz1f(myinput)+isz1b(myinput) << endl;
  cout << "Re SigZ' = " << rsz1fp(myinput)+rsz1bp(myinput) << endl;
  cout << "Im SigZ' = " << isz1fp(myinput)+isz1bp(myinput) << endl;
  cout << "Re SigZ\" = " << rsz1fpp(myinput)+rsz1bpp(myinput) << endl;
  cout << "Im SigZ\" = " << isz1fpp(myinput)+isz1bpp(myinput) << endl;
  cout << endl;
  
  string typenam[5] = {"l", "n", "u", "d", "b"};
  string formnam[2] = {"V", "A"};
  int typenum[i],j;
  cout << "      ";
  for(i=0; i<5; i++)
    cout << setw(13) << typenam[i];
  cout << endl;
  for(j=0; j<2; j++)
  {
    cout << "Re Z" << formnam[j] << " = ";
    for(i=0; i<5; i++)
      cout << setw(13) << rz1f(typenum[i],j,myinput)+rz1b(typenum[i],j,myinput);
    cout << endl;
    cout << "Im Z" << formnam[j] << " = ";
    for(i=0; i<5; i++)
      cout << setw(13) << iz1f(typenum[i],j,myinput)+iz1b(typenum[i],j,myinput);
    cout << endl;
  }
  cout << endl;
  for(j=0; j<2; j++)
  {
    cout << "Re Z'" << formnam[j] << "= ";
    for(i=0; i<5; i++)
      cout << setw(13) << rz1fp(typenum[i],j,myinput)+rz1bp(typenum[i],j,myinput);
    cout << endl;
    cout << "Im Z'" << formnam[j] << "= ";
    for(i=0; i<5; i++)
      cout << setw(13) << iz1fp(typenum[i],j,myinput)+iz1bp(typenum[i],j,myinput);
    cout << endl;
  }
  cout << endl;
  for(j=0; j<2; j++)
  {
    cout << "Re G" << formnam[j] << " = ";
    for(i=0; i<5; i++)
      cout << setw(13) << rg1f(typenum[i],j,myinput)+rg1b(typenum[i],j,myinput);
    cout << endl;
    cout << "Im G" << formnam[j] << " = ";
    for(i=0; i<5; i++)
      cout << setw(13) << ig1f(typenum[i],j,myinput)+ig1b(typenum[i],j,myinput);
    cout << endl;
  }
  cout << endl;
*/     

//  ltini();
  cout << "      ";
  for(i=0; i<5; i++)
    cout << setw(13) << typenam[i];
  cout << endl;
  
  for(j=0; j<2; j++)
    for(k=0; k<2; k++)
    {
      cout << "Re B" << formnam[j] << formnam[k] << " =  ";
      for(i=0; i<5; i++)
        cout << setw(13) << real(B1(ELE,typenum[i],j,k,88*88,0.5,myinput,0,1)-
			B1(ELE,typenum[i],j,k,88*88,0.5,myinput,0,0));
      cout << endl;
    }
  cout << endl;
//  ltexi();
  
  return 0;
}
