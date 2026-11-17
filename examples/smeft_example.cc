/*
  File: smeft_example.cc
  Purpose: Demonstrate SM input values, Wilson coefficients, the finite SMEFT
           calculation, and full/SM observable results through the public API.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

// Complete example of the finite SMEFT interface: NNLO+ SM form factors,
// mat_SMaas residue/background and prescription-generated nonlinear terms.
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <griffin/SMEFT.h>

int main()
{
  using namespace griffin;
  using namespace griffin::prediction;
  try {
    SMvalGmu sm;
    sm.set(InputPar::alpha,1.0/137.03599976);
    sm.set(InputPar::ZMassExperiment,91.1876);
    sm.set(InputPar::ZWidthExperiment,2.4952);
    sm.set(InputPar::HiggsMass,125.1);
    sm.set(InputPar::topMass,172.5);
    sm.set(InputPar::bottomMass,2.87);
    sm.set(InputPar::charmMass,0.652);
    sm.set(InputPar::tauMass,1.777);
    sm.set(InputPar::electronMass,0.0);
    sm.set(InputPar::muonMass,0.0);
    sm.set(InputPar::upMass,0.0);
    sm.set(InputPar::downMass,0.0);
    sm.set(InputPar::strangeMass,0.0);
    sm.set(InputPar::alphaS,0.1179);
    sm.set(InputPar::DeltaAlpha,0.059);
    sm.set(InputPar::fermiConstant,1.166379e-5);

    // Dimensionless Wilson coefficients, an explicit scale in GeV, and zero
    // for every coefficient not set below.
    smeft::Input coefficients;
    coefficients.setLambdaGeV(1000.0);
    coefficients.setCphiWB(0.01); // C/Lambda^2 = 0.01 TeV^-2
    const ElectroweakSetup setup(EWInputScheme::GmuAlphaMZ,sm);
    const SMEFTCalculator calculator(setup,coefficients);
    const Fermion initialFermion = Fermion::electron;
    const Fermion finalFermion = Fermion::muon;
    const double toPicobarns = 0.38937966e9; // convert GeV^-2 to pb
    const ScatteringPoint point(initialFermion,finalFermion,
        90.0*90.0,0.5,toPicobarns);
    const auto sigma = calculator.differentialCrossSection(point);
    const auto fv = calculator.pseudoObservable(PseudoObservable::FV,Fermion::electron);
    std::cout << std::setprecision(15)
      << "Finite Gmu pole MW [GeV]: "
      << calculator.resolvedInput().get(InputPar::WMassComplex) << '\n'
      << "SM d sigma/d cos(theta) [pb]: " << sigma.smResult() << '\n'
      << "Full SMEFT d sigma/d cos(theta) [pb]: " << sigma.total() << '\n'
      << "Full minus zero-coefficient result [pb]: " << sigma.modelCorrection() << '\n'
      << "FV from the same FA and sweff (including SM imaginary part): "
      << real(fv.total()) << '\n';
    return 0;
  } catch(const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
