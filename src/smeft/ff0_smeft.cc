/*
  ff0_smeft.cc: direct SMEFT changes to vector and axial-vector Z couplings

  Conventions follow L_SMEFT = L_SM + Sum_i C_i Q_i/Lambda^2 with
  dimensionless Wilson coefficients.  These are the explicit dimension-six
  tree vertex formulas at the supplied alpha, MW and MZ.  The finite input
  resolution belongs to WMassSolution; it is not repeated or linearized
  here.  SMEFTCalculator uses these shifts in its full finite prediction.
*/

#include "griffin/smeft/ff0_smeft.h"

#include <cmath>
#include <stdexcept>

#include "griffin/sm/ff0.h"
#include "smeft/detail/SMEFTFlavor.h"

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
  arbitrary-Lambda version of Jackson's z0SMEFT() chiral formula.  A Gmu
  calculation supplies the W mass already resolved by WMassSolution.
*/
VertexShift directZShift(int type, const inval& smInput,
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

} // namespace

double directZAxialShift(int type, const inval& smInput, const Input& eftInput)
{
  return directZShift(type,smInput,eftInput).axial;
}

double directZVectorShift(int type, const inval& smInput, const Input& eftInput)
{
  return directZShift(type,smInput,eftInput).vector;
}

} // namespace smeft
} // namespace griffin
