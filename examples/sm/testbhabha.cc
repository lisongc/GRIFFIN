/*
  File: testbhabha.cc
  Purpose: Demonstrate the legacy SM Bhabha-scattering calculation.
  Original SM example retained without numerical changes.
  Merged-version authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16 (typed-only input labels)
*/

#include <iostream>
using namespace std;

#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/xscnnlo.h"
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

  SMval myinput;	// convert masses from PDG values to complex pole scheme
  myinput.set(InputPar::alpha, 1/137.03599976);
  myinput.set(InputPar::ZMassExperiment, 91.1876);
  myinput.set(InputPar::WMassExperiment, 80.377);
  myinput.set(InputPar::ZWidthExperiment, 2.4952);
  myinput.set(InputPar::WWidthExperiment, 2.085);
  myinput.set(InputPar::HiggsMass, 125.1);
  myinput.set(InputPar::topMass, 172.5);
  myinput.set(InputPar::bottomMass, 2.87);
  myinput.set(InputPar::DeltaAlpha, 0.059);
  myinput.set(InputPar::alphaS, 0.1179);
  
  cout << endl << "Complex-pole masses: MW=" << myinput.get(InputPar::WMassComplex) << ", MZ="
    << myinput.get(InputPar::ZMassComplex) << endl << endl;
 
  // compute matrix element for ee->ee with vector coupling in initial
  // state and vector coupling in final state
  const Fermion ini = Fermion::electron, fin = Fermion::electron;
  const Current iff = Current::vector, off = Current::vector;
  
  cout << "=== Matrix element for ee->ee (i=e, f=e) ===" << endl << endl;
  
  // compute vertex form factors:
  FA_SMNNLO FAi(ini, myinput), FAf(fin, myinput);
  SW_SMNNLO SWi(ini, myinput), SWf(fin, myinput);
  cout << "F_A^i (NNLO+) = " << FAi.result() << endl;
  cout << "F_A^f (NNLO+) = " << FAf.result() << endl;
  cout << "sineff^i (NNLO+) = " << SWi.result() << endl;
  cout << "sineff^f (NNLO+) = " << SWf.result() << endl;
  cout << endl;
  
  double cme = 10.,  // center-of-mass energy
         cost = 0.5; // scattering angle
  Cplx res1, res2;

  cout << "SM matrix element M_VV for cos(theta)=" << cost << ": " << endl;
  // compute matrix element for ee->ee using SM form factors:
  mat_SMNNLO M(ini, fin, iff, off, FAi, FAf, SWi, SWf, cme*cme, cost, myinput);
  cout << "sqrt(s)\t\ttot. result\t\toff-resonance contrib." << endl;
  for(cme = 10.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    res1 = M.result();
    res2 = M.resoffZ();
    cout << cme << " \t" << res1 << " \t" << res2 << endl;
  }
  cout << endl;

  cout << "diff. cross-section for cos(theta)=" << cost << ": " << endl;
  // compute diff. cross-section for unpolarized beams from matrix element:  
  Cplx resvv, resva, resav, resaa, resss, respp;
  double xsec,
         GeVtoNB = 0.38937966e6;  // unit conversion from GeV^-2 to nb
  cout << "sqrt(s)\tdsig/dcos [nb]" << endl;
  for(cme = 10.; cme <= 190.; cme += 20.)
  {
    M.setkinvar(cme*cme, cost);
    M.setform(Current::vector, Current::vector);
    resvv = M.result();
    M.setform(Current::axial, Current::vector);
    resav = M.result();
    M.setform(Current::vector, Current::axial);
    resva = M.result();
    M.setform(Current::axial, Current::axial);
    resaa = M.result();
    M.setform(Current::scalar, Current::scalar);
    resss = M.result();
    M.setform(Current::pseudoscalar, Current::pseudoscalar);
    respp = M.result();
    xsec = real((1+cost*cost)*(resvv*conj(resvv) + resav*conj(resav)
    			     + resva*conj(resva) + resaa*conj(resaa)) +
		+ 4*cost*(resvv*conj(resaa) + resva*conj(resav))
		+ resss*conj(resss) + respp*conj(respp));
    xsec *= cme*cme/(32*Pi) * GeVtoNB;
    cout << cme << " \t" << xsec << endl;
  }
  return 0;
}
