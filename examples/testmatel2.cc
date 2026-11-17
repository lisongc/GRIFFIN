// Demonstrate the active v1.2 fermionic two-loop additions for ee -> dd.
// Compare mat_SMaas and mat_SMNNLO2 amplitudes and d sigma / d cos(theta).
// This prints a scan; automated reference assertions are in tests/test_sm_nnlo2.cc.
// mat_SMNNLO2 retains upstream's disabled compoffZfb and resaaa sectors.
#include <iostream>
using namespace std;

#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/xscnnlo2.h"
#include "griffin/sm/SMval.h"
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

  SMval myinput;        // convert masses from PDG values to complex pole scheme
  myinput.set(InputPar::alpha, 1/137.03599976);
  myinput.set(InputPar::ZMassExperiment, 91.1876);
  myinput.set(InputPar::WMassExperiment, 80.377);
  myinput.set(InputPar::ZWidthExperiment, 2.4952);
  myinput.set(InputPar::WWidthExperiment, 2.085);
  myinput.set(InputPar::HiggsMass, 125.1);
  myinput.set(InputPar::topMass, 172.5);
  myinput.set(InputPar::bottomMass, 2.87);
  myinput.set(InputPar::tauMass, 1.777);
  myinput.set(InputPar::charmMass, 1.27);
  myinput.set(InputPar::strangeMass, 0.342);  // for light quarks use thrsh. masses from 1712.09146
  myinput.set(InputPar::muonMass, 0.1057);
  myinput.set(InputPar::upMass, 0.246);
  myinput.set(InputPar::downMass, 0.246);
  myinput.set(InputPar::electronMass, 0.000511);
  myinput.set(InputPar::DeltaAlpha, 0.059);
  myinput.set(InputPar::alphaS, 0.1179);

  // compute matrix element for ee->dd with vector coupling in initial
  // state and vector coupling in final state
  Fermion ini = Fermion::electron, fin = Fermion::d;
  Current iff = Current::vector, off = Current::vector;

  cout << "=== Matrix element for ee->dd (i=e, f=d) ===" << endl << endl;

  // compute vertex form factors:
  FA_SMNNLO FAi(ini, myinput), FAf(fin, myinput);
  SW_SMNNLO SWi(ini, myinput), SWf(fin, myinput);

  double cme = 30.,        // center-of-mass energy
         cost = 0.5; // scattering angle
  Cplx res1, res2;

  cout << "SM matrix element M_VV for cos(theta)=" << cost << ": " << endl;
  // compute matrix element for ee->dd using SM form factors:
  mat_SMNNLO2 M(ini, fin, iff, off, FAi, FAf, SWi, SWf, cme*cme, cost, myinput);
  cout << "sqrt(s)\tmat_SMaas total\tmat_SMNNLO2 total (active sectors)" << endl;
  for(cme = 30.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    res1 = M.mat_SMaas::result();
    res2 = M.result();
    cout << cme << " \t" << res1 << " \t" << res2 << endl;
  }
  cout << endl;
  cout << "sqrt(s)\tmat_SMaas off-resonance\tmat_SMNNLO2 off-resonance (active sectors)" << endl;
  for(cme = 30.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    res1 = M.mat_SMaas::resoffZ();
    res2 = M.resoffZ();
    cout << cme << " \t" << res1 << " \t" << res2 << endl;
  }
  cout << endl;

  cout << "sqrt(s)\tresvpol (active resaaz term)" << endl;
  for(cme = 30.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    res2 = M.resvpol();
    cout << cme << " \t" << res2 << endl;
  }
  cout << endl;

  cout << "diff. cross-section [nb] for cos(theta)=" << cost << ": " << endl;
  // compute diff. cross-section for unpolarized beams from matrix element:
  Cplx res1vv, res1va, res1av, res1aa;
  Cplx res2vv, res2va, res2av, res2aa;
  double xsec,
         GeVtoNB = 0.38937966e6;  // unit conversion from GeV^-2 to nb
  cout << "sqrt(s)\tdsig/dcos mat_SMaas\tdsig/dcos mat_SMNNLO2 (active sectors)" << endl;
  for(cme = 30.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    M.setform(Current::vector, Current::vector);
    res1vv = M.mat_SMaas::result();
    res2vv = M.result();
    M.setform(Current::axial, Current::vector);
    res1av = M.mat_SMaas::result();
    res2av = M.result();
    M.setform(Current::vector, Current::axial);
    res1va = M.mat_SMaas::result();
    res2va = M.result();
    M.setform(Current::axial, Current::axial);
    res1aa = M.mat_SMaas::result();
    res2aa = M.result();
    xsec = real((1+cost*cost)*(res1vv*conj(res1vv) + res1av*conj(res1av)
                             + res1va*conj(res1va) + res1aa*conj(res1aa)) +
		+ 4*cost*(res1vv*conj(res1aa) + res1va*conj(res1av)));
    xsec *= 3*cme*cme/(32*Pi) * GeVtoNB;
    cout << cme << " \t" << xsec << " \t\t\t";
    xsec = real((1+cost*cost)*(res2vv*conj(res2vv) + res2av*conj(res2av)
                             + res2va*conj(res2va) + res2aa*conj(res2aa)) +
		+ 4*cost*(res2vv*conj(res2aa) + res2va*conj(res2av)));
    xsec *= 3*cme*cme/(32*Pi) * GeVtoNB;
    cout << xsec << endl;
  }

  return 0;
}
