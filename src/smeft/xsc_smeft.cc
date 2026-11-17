/*
  File: xsc_smeft.cc
  Purpose: Match four-fermion Wilson coefficients to local vector/axial contact
           amplitudes without repeating electroweak mass resolution.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

/*
  Local dimension-six four-fermion contact amplitudes.

  These fixed-strength contributions are independent of electroweak mass
  resolution. SMEFTCalculator adds them once to its finite pole/background
  amplitude before taking the full square. No linearized vertex residue or
  W-mass response is implemented here.
*/

#include "griffin/smeft/xsc_smeft.h"

#include <stdexcept>

namespace griffin {
namespace smeft {
namespace detail {

struct ChiralContact {
  double ll;
  double lr;
  double rl;
  double rr;
  static ChiralContact forProcess(Fermion intype, Fermion outtype, const Input& input);
};

ChiralContact ChiralContact::forProcess(
    Fermion intype, Fermion outtype, const Input& input)
{
  if(intype != Fermion::electron)
    throw std::invalid_argument(
        "SMEFT contact amplitudes require an electron initial state");
  if(outtype == Fermion::electron)
    throw std::domain_error(
        "Bhabha SMEFT contact terms require separate s/t-channel treatment");

  const Input::FermionInfo info = Input::lookupFermion(outtype);
  const Generation one = Generation::first;
  const Generation i = info.generation;

  switch(info.category)
  {
    case Input::FermionCategory::upQuark:
      return {
        input.Clq1(one,one,i,i)-input.Clq3(one,one,i,i),
        input.Clu(one,one,i,i),
        input.Cqe(i,i,one,one),
        input.Ceu(one,one,i,i)
      };

    case Input::FermionCategory::downQuark:
      return {
        input.Clq1(one,one,i,i)+input.Clq3(one,one,i,i),
        input.Cld(one,one,i,i),
        input.Cqe(i,i,one,one),
        input.Ced(one,one,i,i)
      };

    case Input::FermionCategory::chargedLepton:
      // The electron final state was rejected above.  For muons and taus the
      // two Cll contractions are independent gauge-invariant inputs.
      // Each getter is already the total coefficient of its unique operator;
      // do not also add the equivalent reversed-current/Fierz spellings.
      return {
        input.Cll(one,one,i,i)+input.Cll(one,i,i,one),
        input.Cle(one,one,i,i),
        input.Cle(i,i,one,one),
        input.Cee(one,one,i,i)
      };

    case Input::FermionCategory::neutrino:
      if(i == Generation::first)
      {
        // Qll[1111] contains two identical ee-nu_e-nu_e components in the
        // operator normalization used by the contact matching below.
        return {
          2.0*input.Cll(one,one,one,one),
          0.0,
          input.Cle(one,one,one,one),
          0.0
        };
      }
      return {
        input.Cll(one,one,i,i),
        0.0,
        input.Cle(i,i,one,one),
        0.0
      };
  }

  throw std::logic_error("unreachable SMEFT fermion category");
}

} // namespace detail

namespace {

void validateVectorAxialCurrents(Current inform, Current outform)
{
  if((inform != Current::vector && inform != Current::axial)
      || (outform != Current::vector && outform != Current::axial))
    throw std::invalid_argument(
        "SMEFT contact amplitudes support only vector and axial-vector currents");
}

double projectContact(const detail::ChiralContact& coefficient,
                      Current inform, Current outform)
{
  if(inform == Current::vector && outform == Current::vector)
    return coefficient.ll+coefficient.lr+coefficient.rl+coefficient.rr;
  if(inform == Current::vector && outform == Current::axial)
    return coefficient.ll-coefficient.lr+coefficient.rl-coefficient.rr;
  if(inform == Current::axial && outform == Current::vector)
    return coefficient.ll+coefficient.lr-coefficient.rl-coefficient.rr;
  return coefficient.ll-coefficient.lr-coefficient.rl+coefficient.rr;
}

} // namespace

Cplx fourFermionContactAmplitude(Fermion intype, Fermion outtype,
                                Current inform, Current outform,
                                const Input& eftInput)
{
  validateVectorAxialCurrents(inform,outform);
  const detail::ChiralContact coefficient =
      detail::ChiralContact::forProcess(intype,outtype,eftInput);
  return 0.25*eftInput.inverseLambdaSquaredGeV()
      *projectContact(coefficient,inform,outform);
}

} // namespace smeft
} // namespace griffin
