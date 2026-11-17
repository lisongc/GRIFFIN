/*
  File: EWPO_smeft.cc
  Purpose: Evaluate SMEFT FA/SW/FV at already-resolved electroweak inputs.
           The calculator uses the same FA/SW for its observable results and
           scattering pole inputs. These private helpers neither solve for MW
           nor introduce a separate public evaluator.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

#include "griffin/smeft/SMEFTCalculator.h"

#include <cmath>
#include <stdexcept>

#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/ff.h"
#include "griffin/sm/ff0.h"
#include "griffin/smeft/ff0_smeft.h"

namespace griffin {
namespace prediction {

SMEFTCalculator::PoleInputs SMEFTCalculator::poleInputs(
    Fermion fermion, const inval& input, const smeft::Input& coefficients)
{
  const Fermion type = fermion;
  // The input already contains the resolved MW, including its finite Gmu
  // shift.  These helpers add only the explicit vertex changes at that point.
  // They validate the particle before the charge arrays are indexed below.
  const double da = smeft::directZAxialShift(type,input,coefficients);
  const double dv = smeft::directZVectorShift(type,input,coefficients);
  const double a0 = az0(type,input);
  const double v0 = vz0(type,input);
  PoleInputs result;
  // These explicit FA/SW additions are linear in Wilson coefficients at a
  // fixed mass point. In particular FA has no explicit da*da term. The SM
  // functions are NNLO+ and use the resolved mass. Downstream pole/FV
  // construction and full amplitude squaring retain their generated nonlinear
  // terms; this is not a prescription for all higher powers of the coefficients.
  result.fa = FA_SMNNLO(type,input).result()+2.0*a0*da;
  result.sw = 0.0; // unused for neutrinos; their public SW request is rejected
  if(!isNeutrino(fermion))
    result.sw = SW_SMNNLO(type,input).result()
      -(dv/a0-v0*da/(a0*a0))/(4.0*std::fabs(Qf[index(type)]));
  requireFinite(result.fa,"FA");
  requireFinite(result.sw,"SW");
  if(real(result.fa) <= 0.0)
    throw std::domain_error("SMEFT FA must be positive for the pole construction");
  return result;
}

Cplx SMEFTCalculator::observableAt(PseudoObservable observable, Fermion fermion,
                                  const inval& input,
                                  const smeft::Input& coefficients)
{
  if(observable == PseudoObservable::WMass)
    return input.get(InputPar::WMassComplex);
  if(observable == PseudoObservable::DeltaR)
    throw std::domain_error("finite DeltaR is not defined by SMEFTCalculator; "
                            "use massSolution() for the Gmu equation diagnostics");
  if(observable == PseudoObservable::SW && isNeutrino(fermion))
    throw std::domain_error("SW is undefined for a neutral fermion");
  const PoleInputs pole = poleInputs(fermion,input,coefficients);
  if(observable == PseudoObservable::FA) return pole.fa;
  if(observable == PseudoObservable::SW) return pole.sw;
  if(observable != PseudoObservable::FV)
    throw std::invalid_argument("unsupported SMEFT pseudo-observable");

  // Use the FV_SMNNLO construction with the same shifted FA and SW
  // as the pole amplitude. Retain the SM imaginary-ratio term at resolved MW;
  // otherwise even the zero-coefficient limit would differ from GRIFFIN SM.
  const Fermion type = fermion;
  const double a0 = az0(type,input), v0 = vz0(type,input);
  const double imaginaryRatio =
      (a0*(iz1f(type,Current::vector,input)+iz1b(type,Current::vector,input))
       -v0*(iz1f(type,Current::axial,input)+iz1b(type,Current::axial,input)))/(a0*a0);
  const double q = 1.0-4.0*std::fabs(Qf[index(type)])*real(pole.sw);
  const Cplx fv = pole.fa*(q*q+imaginaryRatio*imaginaryRatio);
  requireFinite(fv,"FV");
  return fv;
}

} // namespace prediction
} // namespace griffin
