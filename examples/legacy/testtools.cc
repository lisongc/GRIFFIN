#include <iostream>
using namespace std;

#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/SMval.h"
#include "griffin/sm/tools.h"
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

  static double mhlist[4] = {100, 200, 600, 1000};
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
  myinput.set(InputPar::bottomMass, 2.877); // MSbar mass at scale mu=MZ for mb(mb)=4.20
  myinput.set(InputPar::charmMass, 0.652); // MSbar mass at scale mu=MZ
  myinput.set(InputPar::strangeMass, 0);
  myinput.set(InputPar::upMass, 0);
  myinput.set(InputPar::downMass, 0);
  myinput.set(InputPar::tauMass, 1.777);
  myinput.set(InputPar::muonMass, 0);
  myinput.set(InputPar::electronMass, 0);
  myinput.set(InputPar::DeltaAlpha, 0.05900);
  myinput.set(InputPar::fermiConstant, 1.166379e-5);

  FV_SMNNLO fv2(Fermion::electron,myinput);
  FA_SMNNLO fa2(Fermion::electron,myinput);
  // flavor selection here will be overwritten by partzwidth/zwidth below

  string typenam[5] = {"l", "n", "u", "d", "b"};
  Fermion typenum[5] = {Fermion::electron, Fermion::nuElectron, Fermion::u, Fermion::d, Fermion::b};

  for(i = 0; i < 5; i++)
    cout << "Gamma[Z->" << typenam[i] << typenam[i] << "] = "
         << partzwidth(fa2, fv2, typenum[i], myinput, RUNWIDTHSCHEME) << endl;
  cout << "Gamma[Z] = " << zwidth(fa2, fv2, myinput, RUNWIDTHSCHEME) << endl;
  cout << endl;

  return 0;
}
