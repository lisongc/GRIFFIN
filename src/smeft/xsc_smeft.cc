/*
  Local dimension-six four-fermion contact amplitudes.

  These fixed-strength contributions are independent of electroweak mass
  resolution. SMEFTCalculator adds them once to its finite pole/background
  amplitude before taking the full square. No linearized vertex residue or
  W-mass response is implemented here.
*/

#include "griffin/smeft/xsc_smeft.h"

#include <stdexcept>

#include "smeft/detail/SMEFTFlavor.h"

namespace griffin {
namespace smeft {
namespace {

struct ChiralContact {
  double ll;
  double lr;
  double rl;
  double rr;
};

void validateVectorAxialCurrents(int inform, int outform)
{
  if((inform != VEC && inform != AXV)
      || (outform != VEC && outform != AXV))
    throw std::invalid_argument(
        "SMEFT contact S supports only vector and axial-vector currents");
}

ChiralContact chiralContact(int intype, int outtype, const Input& input)
{
  if(intype != ELE)
    throw std::invalid_argument(
        "SMEFT contact S formulas require an electron initial state");
  if(outtype == ELE)
    throw std::domain_error(
        "Bhabha SMEFT contact terms require separate s/t-channel treatment");

  const detail::FermionInfo info = detail::fermionInfo(outtype);
  const Generation one = Generation::first;
  const Generation i = info.generation;

  switch(info.family)
  {
    case detail::FermionFamily::upQuark:
      return {
        input.Clq1(one,one,i,i)-input.Clq3(one,one,i,i),
        input.Clu(one,one,i,i),
        input.Cqe(i,i,one,one),
        input.Ceu(one,one,i,i)
      };

    case detail::FermionFamily::downQuark:
      return {
        input.Clq1(one,one,i,i)+input.Clq3(one,one,i,i),
        input.Cld(one,one,i,i),
        input.Cqe(i,i,one,one),
        input.Ced(one,one,i,i)
      };

    case detail::FermionFamily::chargedLepton:
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

    case detail::FermionFamily::neutrino:
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

  throw std::logic_error("unreachable SMEFT fermion family");
}

double projectContact(const ChiralContact& coefficient,
                      int inform, int outform)
{
  if(inform == VEC && outform == VEC)
    return coefficient.ll+coefficient.lr+coefficient.rl+coefficient.rr;
  if(inform == VEC && outform == AXV)
    return coefficient.ll-coefficient.lr+coefficient.rl-coefficient.rr;
  if(inform == AXV && outform == VEC)
    return coefficient.ll+coefficient.lr-coefficient.rl-coefficient.rr;
  return coefficient.ll-coefficient.lr-coefficient.rl+coefficient.rr;
}

} // namespace

Cplx coeffS4fD6Tree(int intype, int outtype, int inform, int outform,
                    const Input& eftInput)
{
  validateVectorAxialCurrents(inform,outform);
  const ChiralContact coefficient = chiralContact(intype,outtype,eftInput);
  return 0.25*eftInput.inverseLambdaSquaredGeV()
      *projectContact(coefficient,inform,outform);
}

Cplx coeffSp4fD6Tree(int intype, int outtype, int inform, int outform,
                     const Input& eftInput)
{
  validateVectorAxialCurrents(inform,outform);
  // Validate the flavor/process domain even though the local derivative is 0.
  chiralContact(intype,outtype,eftInput);
  return 0.0;
}

} // namespace smeft
} // namespace griffin
