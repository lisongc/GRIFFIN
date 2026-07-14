/*
  xsc_smeft.cc: SM LO plus linear dimension-six matrix elements

  The implementation separates three physically different pieces:

    * coeffRD6Tree: shifts of the two resonant Z vertices;
    * coeffS4fD6Tree: local vector four-fermion contact interactions;
    * coeffSp4fD6Tree: zero because those local contacts are s independent.

  The numerical result is assembled only as poleD6Tree()+resoffZD6Tree().
  The separately exposed S and S' coefficients are never added to that path,
  preventing a local contact term from being counted both as S and as an exact
  background.  The background interface also leaves room for process-specific
  nonlocal terms such as t-channel W exchange in e+e- -> nu_e nubar_e.
*/

#include "smeft/xsc_smeft.h"

#include <stdexcept>

#include "ff0.h"
#include "smeft/SMEFTFlavor.h"
#include "smeft/ff0_smeft.h"
#include "xscmassless.h"

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
        // operator normalization used by the revised matching result.
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

void requireCompleteRegularBackground(int intype, int outtype)
{
  if(intype == ELE && outtype == NUE)
    throw std::domain_error(
        "e+e- -> nu_e nubar_e requires the process-specific t-channel W "
        "background; use the pole and local-contact coefficients separately");
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

Cplx mat_SMEFTLO::coeffRD6Tree() const
{
  if(iff > AXV || off > AXV)
    return 0.0;

  // The D6 residue is a tree-level object.  Its reference factors must remain
  // tree level even when the supplied FA/SW objects carry a higher-order SM
  // prediction for the separate SM amplitude.
  const Cplx ziSM = z0(it,iff,*ival);
  const Cplx zfSM = z0(ot,off,*ival);
  const double dzi = zD6Tree(it,iff,*ival,*eftInput,scheme);
  const double dzf = zD6Tree(ot,off,*ival,*eftInput,scheme);

  // Keep exactly one D6 insertion.  The omitted dzi*dzf term is O(Lambda^-4).
  return dzi*zfSM+ziSM*dzf;
}

Cplx mat_SMEFTLO::coeffR() const
{
  return matel::coeffR()+coeffRD6Tree();
}

Cplx mat_SMEFTLO::coeffS4fD6Tree() const
{
  return smeft::coeffS4fD6Tree(it,ot,iff,off,*eftInput);
}

Cplx mat_SMEFTLO::coeffSp4fD6Tree() const
{
  return smeft::coeffSp4fD6Tree(it,ot,iff,off,*eftInput);
}

Cplx mat_SMEFTLO::coeffS() const
{
  return matel::coeffS()+coeffS4fD6Tree();
}

Cplx mat_SMEFTLO::coeffSp() const
{
  return matel::coeffSp()+coeffSp4fD6Tree();
}

Cplx mat_SMEFTLO::poleD6Tree() const
{
  const double mz = ival->get(MZ);
  const double gz = ival->get(GamZ);
  const Cplx sMinusPole(s-mz*mz,mz*gz);
  return coeffRD6Tree()/sMinusPole;
}

Cplx mat_SMEFTLO::resoffZ4fD6Tree() const
{
  // A local tree-level contact term is constant, so its exact regular
  // amplitude equals its Laurent S coefficient at every s and angle.
  return coeffS4fD6Tree();
}

Cplx mat_SMEFTLO::resoffZD6Tree() const
{
  requireCompleteRegularBackground(it,ot);
  return resoffZ4fD6Tree();
}

Cplx mat_SMEFTLO::resultD6Tree() const
{
  requireCompleteRegularBackground(it,ot);
  return poleD6Tree()+resoffZD6Tree();
}

Cplx mat_SMEFTLO::resoffZ() const
{
  return matel::resoffZ()+resoffZD6Tree();
}

Cplx mat_SMEFTLO::result() const
{
  // matel::result() contains the SM pole and exact SM background.  The D6
  // correction is then added once through its own pole/background split.
  return matel::result()+resultD6Tree();
}

DifferentialCrossSectionLinear differentialCrossSectionLinear(
    int intype, int outtype,
    const psobs& FAin, const psobs& FAout,
    const psobs& SWin, const psobs& SWout,
    double sval, double costheta, const inval& smInput,
    const Input& eftInput, EWInputScheme inputScheme,
    double unitConversion)
{
  const MasslessFermionCrossSectionBuilder builder(
      static_cast<Fermion>(intype),static_cast<Fermion>(outtype),
      sval,costheta,unitConversion);

  matel smSelected(intype,outtype,VEC,VEC,FAin,FAout,SWin,SWout,
                   sval,costheta,smInput);

  FA_SMLO faTreeIn(intype,smInput);
  FA_SMLO faTreeOut(outtype,smInput);
  SW_SMLO swTreeIn(intype,smInput);
  SW_SMLO swTreeOut(outtype,smInput);
  matel smTree(intype,outtype,VEC,VEC,
               faTreeIn,faTreeOut,swTreeIn,swTreeOut,
               sval,costheta,smInput);
  mat_SMEFTLO eftTree(intype,outtype,VEC,VEC,
                      faTreeIn,faTreeOut,swTreeIn,swTreeOut,
                      sval,costheta,smInput,eftInput,inputScheme);

  MasslessFermionCrossSectionBuilder::Amplitudes selected = {};
  MasslessFermionCrossSectionBuilder::Amplitudes tree = {};
  MasslessFermionCrossSectionBuilder::Amplitudes d6 = {};
  for(std::size_t current = 0; current < builder.channelCount(); ++current)
  {
    const Current initialCurrent = builder.initialCurrent(current);
    const Current finalCurrent = builder.finalCurrent(current);
    smSelected.setform(initialCurrent,finalCurrent);
    smTree.setform(initialCurrent,finalCurrent);
    eftTree.setform(initialCurrent,finalCurrent);
    selected[current] = smSelected.result();
    tree[current] = smTree.result();
    d6[current] = eftTree.resultD6Tree();
  }

  DifferentialCrossSectionLinear result;
  result.sm = builder.squaredPrediction(selected);
  result.d6 = builder.linearInterference(tree,d6);
  return result;
}

} // namespace smeft
} // namespace griffin
