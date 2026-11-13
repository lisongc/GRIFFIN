#ifndef GRIFFIN_SMEFT_CALCULATOR_H
#define GRIFFIN_SMEFT_CALCULATOR_H

#include "griffin/ElectroweakSetup.h"
#include "griffin/Observables.h"
#include "griffin/smeft/NativeInputPoint.h"
#include "griffin/core/FermionPairCrossSection.h"

namespace griffin {
namespace prediction {

// Two complete evaluations: reference() has all Wilson coefficients zero;
// total() uses the supplied coefficients. modelCorrection() is their finite
// difference, not just the term linear in Wilson coefficients. The total is
// stored directly, not reconstructed by adding the difference to the reference.
template<typename Value>
class FinitePrediction {
  Value reference_, total_;
public:
  FinitePrediction(const Value& reference, const Value& total)
    : reference_(reference), total_(total) {}
  const Value& reference() const { return reference_; }
  const Value& total() const { return total_; }
  Value modelCorrection() const { return total_-reference_; }
};

/*
  SMEFT evaluator using typed coefficient inputs, Jackson's finite Gmu mass
  solution, shifted pole observables and the square of the complete amplitude.

  Fixed SM content: FA_SMNNLO and SW_SMNNLO (NNLO plus selected higher-order
  terms), dr_SMNNLO for Gmu mass resolution, and mat_SMaas for scattering.
  mat_SMaas adds selected fermionic mixed electroweak-QCD terms to both the
  pole residue and the off-resonance background; it is not a complete NNLO
  off-resonance calculation. Full and zero-coefficient predictions use the
  same prescription, with independently resolved inputs.

  The explicit dimension-six FA/SW input shifts are linear at fixed masses;
  FA contains 2*a0*deltaA but no explicit deltaA^2. Finite mass evaluation,
  pole construction and full squaring retain their generated nonlinear
  coefficient dependence. This is not a complete expansion to every higher
  power of the coefficients. FV uses the same shifted FA/SW and the existing
  SM imaginary-ratio term. See docs/SMEFT_IMPLEMENTATION.md for the formulas.

  Both inputs are copied at construction: later caller mutations do not alter
  this prediction. Construct a new object for a new fit point. Copies/moves
  are safe because this object owns values, not pointers into another object.
  The direct vertex helpers contain no separate linear W-mass response; all
  mass dependence is evaluated at the resolved finite input point.
*/
class SMEFTCalculator {
public:
  typedef FermionPairCrossSection::Amplitudes Amplitudes;

  SMEFTCalculator(
      const ElectroweakSetup& setup, const smeft::Input& coefficients,
      const smeft::NativeSolverOptions& options = smeft::NativeSolverOptions());

  EWInputScheme inputScheme() const { return scheme_; }
  const smeft::Input& wilsonInput() const { return coefficients_; }
  const inval& resolvedInput() const { return point_.smInput(); }
  const inval& referenceInput() const { return referencePoint_.smInput(); }
  const smeft::NativeInputPoint& massSolution() const { return point_; }

  // WMass returns the pole mass. SW is undefined for neutrinos.
  // DeltaR is deliberately not exposed: Jackson's additive mass-equation
  // displacement must not be mislabeled as a finite Delta-r observable.
  FinitePrediction<Cplx> pseudoObservable(
      PseudoObservable observable, Fermion fermion) const;

  // Order of channels: VV, AV, VA, AA. amplitudes() includes contacts once;
  // referenceAmplitudes() uses zero coefficients; contactAmplitudes() returns
  // only the contact term and must not be added to amplitudes() a second time.
  Amplitudes amplitudes(const ScatteringPoint& point) const;
  Amplitudes referenceAmplitudes(const ScatteringPoint& point) const;
  Amplitudes contactAmplitudes(const ScatteringPoint& point) const;
  FinitePrediction<Cplx> amplitude(
      const ScatteringPoint& point, const CurrentChannel& channel) const;

  // Unpolarized d sigma / d cos(theta), with massless external fermions;
  // GeV^-2 unless point.unitConversion is supplied. This is the angular square
  // of the amplitude, not a radiation-inclusive detector-level observable.
  // The exact real Z pole is rejected, but nearby mat_SMaas background values
  // also have unresolved limitations; see docs/KNOWN_ISSUES.md before scanning.
  FinitePrediction<double> differentialCrossSection(
      const ScatteringPoint& point) const;

private:
  EWInputScheme scheme_;
  smeft::Input coefficients_, zeroCoefficients_;
  smeft::NativeInputPoint referencePoint_, point_;
};

} // namespace prediction
} // namespace griffin
#endif
