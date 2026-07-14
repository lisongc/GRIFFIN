/*
  ff0_smeft.cc: tree-level SMEFT shifts of vector and axial-vector Z factors

  Conventions follow L_SMEFT = L_SM + Sum_i C_i Q_i/Lambda^2 with
  dimensionless Wilson coefficients.  Both supported electroweak input
  schemes are truncated consistently at O(1/Lambda^2).

  The direct vertex correction is written in the {alpha, MW, MZ} scheme.  For
  GmuAlphaMZ it is supplemented by the first-order displacement of the derived
  MW.  This follows Jackson's input bookkeeping while keeping the truncation
  explicit: the SM input is not mutated and no product of Wilson coefficients
  is ever evaluated.
*/

#include "smeft/ff0_smeft.h"

#include <cmath>
#include <stdexcept>

#include "ff0.h"
#include "smeft/SMEFTFlavor.h"

namespace griffin {
namespace smeft {
namespace {

struct ChiralWilsonCoefficients {
  double right;
  double leftSinglet;
  double leftTriplet;
};

struct VertexShift {
  double vector;
  double axial;
};

ChiralWilsonCoefficients chiralWilsonCoefficients(
    const detail::FermionInfo& info, const Input& input)
{
  const Generation i = info.generation;
  ChiralWilsonCoefficients coefficients;
  switch(info.family)
  {
    case detail::FermionFamily::neutrino:
      coefficients.right = 0.0;
      coefficients.leftSinglet = input.Cphil1(i,i);
      coefficients.leftTriplet = input.Cphil3(i,i);
      return coefficients;
    case detail::FermionFamily::chargedLepton:
      coefficients.right = input.Cphie(i,i);
      coefficients.leftSinglet = input.Cphil1(i,i);
      coefficients.leftTriplet = input.Cphil3(i,i);
      return coefficients;
    case detail::FermionFamily::upQuark:
      coefficients.right = input.Cphiu(i,i);
      coefficients.leftSinglet = input.Cphiq1(i,i);
      coefficients.leftTriplet = input.Cphiq3(i,i);
      return coefficients;
    case detail::FermionFamily::downQuark:
      coefficients.right = input.Cphid(i,i);
      coefficients.leftSinglet = input.Cphiq1(i,i);
      coefficients.leftTriplet = input.Cphiq3(i,i);
      return coefficients;
  }
  throw std::logic_error("unreachable fermion family");
}

/*
  Direct SMEFT correction with {alpha, MW, MZ} held fixed.  This is the signed,
  arbitrary-Lambda version of Jackson's z0SMEFT() chiral formula.  The Gmu
  conversion is added separately through inducedDerivedMWShift() below.
*/
VertexShift directAlphaMWMZShift(int type, const inval& smInput,
                                 const Input& eftInput)
{
  const detail::FermionInfo info = detail::fermionInfo(type);
  const ChiralWilsonCoefficients coefficients =
      chiralWilsonCoefficients(info,eftInput);

  const double alpha = smInput.get(InputPar::alpha);
  const double mw = smInput.get(InputPar::WMassComplex);
  const double mz = smInput.get(InputPar::ZMassComplex);
  if(alpha <= 0.0 || mw <= 0.0 || mz <= 0.0)
    throw std::domain_error("electroweak vertex inputs must be positive");

  const double mw2 = mw*mw;
  const double mz2 = mz*mz;
  if(mw2 >= mz2)
    throw std::domain_error("MW/MZ does not define a physical weak angle");

  const double x = mw2/mz2;
  const double sw2 = 1.0-x;
  const double sw = std::sqrt(sw2);
  const double cw = std::sqrt(x);
  const double e = std::sqrt(4.0*Pi*alpha);
  const double invLambda2 = eftInput.inverseLambdaSquaredGeV();
  const double t3 = I3f[type];
  const double charge = Qf[type];

  const double deltaLeft =
      eftInput.CphiD()*invLambda2*mw2/(8.0*Pi*alpha)
          *(2.0*(t3-charge)-2.0*(2.0*t3-charge)*x)
      + coefficients.leftSinglet*invLambda2*sw2*mw2/(2.0*Pi*alpha)
      - t3*coefficients.leftTriplet*invLambda2*sw2*mw2/(Pi*alpha)
      - t3*eftInput.CphiWB()*invLambda2*mw2*mw
          *std::sqrt(mz2-mw2)/(Pi*alpha*mz2);
  const double deltaRight =
      coefficients.right*invLambda2*sw2*mw2/(2.0*Pi*alpha)
      - charge*eftInput.CphiD()*invLambda2*sw2*mw2/(4.0*Pi*alpha);

  const double normalization = -e/(2.0*sw*cw);
  VertexShift shift;
  shift.vector = normalization*(deltaLeft+deltaRight);
  shift.axial = normalization*(deltaLeft-deltaRight);
  return shift;
}

/*
  Shift of the ordinary SM tree vertex induced by the Gmu-derived
  delta x, x=MW^2/MZ^2.  The derivative is analytic, so the returned term is
  exactly linear even though the unexpanded SM vertex is nonlinear in MW.
*/

VertexShift inducedDerivedMWShift(int type, const inval& smInput,
                                  const Input& eftInput,
                                  EWInputScheme scheme)
{
  VertexShift shift = {0.0,0.0};
  if(scheme == EWInputScheme::AlphaMWMZ)
    return shift;
  if(scheme != EWInputScheme::GmuAlphaMZ)
    throw std::invalid_argument("unsupported SMEFT electroweak input scheme");

  const double alpha = smInput.get(InputPar::alpha);
  const double mw = smInput.get(InputPar::WMassComplex);
  const double mz = smInput.get(InputPar::ZMassComplex);
  const double x = mw*mw/(mz*mz);
  if(x <= 0.0 || x >= 1.0)
    throw std::domain_error("derived MW/MZ does not define a physical weak angle");

  const double e = std::sqrt(4.0*Pi*alpha);
  const double sw = std::sqrt(1.0-x);
  const double cw = std::sqrt(x);
  const double axial0 = I3f[type]*e/(2.0*sw*cw);
  const double vectorRatio =
      1.0-4.0*std::fabs(Qf[type])*(1.0-x);
  const double deltaX = derivedMWShiftD6(smInput,eftInput).deltaMW2OverMZ2;
  const double deltaAxial = axial0*(2.0*x-1.0)
      /(2.0*x*(1.0-x))*deltaX;

  shift.axial = deltaAxial;
  shift.vector = vectorRatio*deltaAxial
      + 4.0*std::fabs(Qf[type])*axial0*deltaX;
  return shift;
}

} // namespace

double azD6Tree(int type, const inval& smInput, const Input& eftInput)
{
  return azD6Tree(type,smInput,eftInput,EWInputScheme::GmuAlphaMZ);
}

double azD6Tree(int type, const inval& smInput, const Input& eftInput,
                EWInputScheme scheme)
{
  const VertexShift direct = directAlphaMWMZShift(type,smInput,eftInput);
  const VertexShift inputShift =
      inducedDerivedMWShift(type,smInput,eftInput,scheme);
  return direct.axial+inputShift.axial;
}

double vzD6Tree(int type, const inval& smInput, const Input& eftInput)
{
  return vzD6Tree(type,smInput,eftInput,EWInputScheme::GmuAlphaMZ);
}

double vzD6Tree(int type, const inval& smInput, const Input& eftInput,
                EWInputScheme scheme)
{
  const VertexShift direct = directAlphaMWMZShift(type,smInput,eftInput);
  const VertexShift inputShift =
      inducedDerivedMWShift(type,smInput,eftInput,scheme);
  return direct.vector+inputShift.vector;
}

double zD6Tree(int type, int form, const inval& smInput,
               const Input& eftInput)
{
  return zD6Tree(type,form,smInput,eftInput,
                 EWInputScheme::GmuAlphaMZ);
}

double zD6Tree(int type, int form, const inval& smInput,
               const Input& eftInput, EWInputScheme scheme)
{
  if(form == VEC)
    return vzD6Tree(type,smInput,eftInput,scheme);
  if(form == AXV)
    return azD6Tree(type,smInput,eftInput,scheme);
  return 0.0;
}

} // namespace smeft
} // namespace griffin
