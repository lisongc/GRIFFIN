/* EWInputScheme.cc: common electroweak reference parameters for SMEFT */

#include "smeft/EWInputScheme.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "deltar.h"

namespace griffin {
namespace smeft {

const char* inputSchemeName(EWInputScheme scheme)
{
  switch(scheme)
  {
    case EWInputScheme::GmuAlphaMZ:
      return "{Gmu, alpha, MZ}";
    case EWInputScheme::AlphaMWMZ:
      return "{alpha, MW, MZ}";
  }
  throw std::invalid_argument("unsupported SMEFT electroweak input scheme");
}

namespace detail {

EWTreeParameters ewTreeParameters(EWInputScheme scheme,
                                  const inval& smInput,
                                  const Input& eftInput)
{
  const double alpha = smInput.get(InputPar::alpha);
  const double mw = smInput.get(InputPar::WMassComplex);
  const double mz = smInput.get(InputPar::ZMassComplex);

  if(alpha <= 0.0 || mw <= 0.0 || mz <= 0.0)
    throw std::domain_error("electroweak inputs must be positive");

  const double cw0 = mw/mz;
  const double cw02 = cw0*cw0;
  if(cw02 <= 0.0 || cw02 >= 1.0)
    throw std::domain_error("MW/MZ does not define a physical weak angle");

  const double sw0 = std::sqrt(1.0-cw02);
  const double e0 = std::sqrt(4.0*Pi*alpha);
  double vevSquared = 0.0;

  switch(scheme)
  {
    case EWInputScheme::GmuAlphaMZ:
    {
      const double gmu = smInput.get(InputPar::fermiConstant);
      if(gmu <= 0.0)
        throw std::domain_error("Gmu must be positive in GmuAlphaMZ");
      vevSquared = 1.0/(std::sqrt(2.0)*gmu);
      break;
    }
    case EWInputScheme::AlphaMWMZ:
    {
      const double vAlpha = 2.0*mw*sw0/e0;
      vevSquared = vAlpha*vAlpha;
      break;
    }
    default:
      throw std::invalid_argument("unsupported SMEFT electroweak input scheme");
  }

  EWTreeParameters parameters;
  parameters.e0 = e0;
  parameters.sw0 = sw0;
  parameters.cw0 = cw0;
  parameters.gZ0 = e0/(sw0*cw0);
  parameters.vevSquared = vevSquared;
  parameters.epsilonD6 =
      vevSquared*eftInput.inverseLambdaSquaredGeV();
  return parameters;
}

double deltaMu(const Input& input)
{
  const Generation first = Generation::first;
  const Generation second = Generation::second;
  return input.Cphil3(first,first) + input.Cphil3(second,second)
      - 0.5*(input.Cll(first,second,second,first)
             + input.Cll(second,first,first,second));
}

double deltaRSMEFTTree(EWInputScheme scheme,
                       const inval& smInput,
                       const Input& eftInput)
{
  double sw02 = 0.0;
  double cw02 = 0.0;
  double epsilonD6 = 0.0;

  switch(scheme)
  {
    case EWInputScheme::GmuAlphaMZ:
    {
      const double alpha = smInput.get(InputPar::alpha);
      const double mz = smInput.get(InputPar::ZMassComplex);
      const double gmu = smInput.get(InputPar::fermiConstant);
      if(alpha <= 0.0 || mz <= 0.0 || gmu <= 0.0)
        throw std::domain_error("GmuAlphaMZ Delta r inputs must be positive");

      const double sw2cw2 = Pi*alpha/(std::sqrt(2.0)*gmu*mz*mz);
      if(sw2cw2 <= 0.0 || sw2cw2 >= 0.25)
        throw std::domain_error(
            "GmuAlphaMZ inputs do not define a physical tree-level weak angle");

      const double root = std::sqrt(0.25-sw2cw2);
      cw02 = 0.5+root;
      sw02 = 0.5-root;
      epsilonD6 = eftInput.inverseLambdaSquaredGeV()
          /(std::sqrt(2.0)*gmu);
      break;
    }
    case EWInputScheme::AlphaMWMZ:
    {
      const EWTreeParameters ew =
          ewTreeParameters(scheme,smInput,eftInput);
      sw02 = ew.sw0*ew.sw0;
      cw02 = ew.cw0*ew.cw0;
      epsilonD6 = ew.epsilonD6;
      break;
    }
    default:
      throw std::invalid_argument("unsupported SMEFT Delta r input scheme");
  }

  const double bracket = 2.0*std::sqrt(cw02/sw02)*eftInput.CphiWB()
      + cw02/(2.0*sw02)*eftInput.CphiD()
      + deltaMu(eftInput);
  return epsilonD6*bracket;
}

} // namespace detail

DerivedMWShiftD6 derivedMWShiftD6(const inval& smInput,
                                  const Input& eftInput)
{
  const double alpha = smInput.get(InputPar::alpha);
  const double mz = smInput.get(InputPar::ZMassComplex);
  const double mw = smInput.get(InputPar::WMassComplex);
  const double gmu = smInput.get(InputPar::fermiConstant);
  if(alpha <= 0.0 || mz <= 0.0 || mw <= 0.0 || gmu <= 0.0)
    throw std::domain_error("GmuAlphaMZ derived-MW inputs must be positive");

  const double k = Pi*alpha/(std::sqrt(2.0)*gmu*mz*mz);
  if(k <= 0.0 || k >= 0.25)
    throw std::domain_error(
        "GmuAlphaMZ inputs do not define a physical tree-level MW shift");

  const double deltaR = detail::deltaRSMEFTTree(
      EWInputScheme::GmuAlphaMZ,smInput,eftInput);
  const double inputDeltaX =
      -k/(2.0*std::sqrt(0.25-k))*deltaR;

  /*
    Jackson solves the loop-improved fixed-point equation

      y = MZ^2 [ 1/2 + sqrt(1/4-k(1+Delta r_SM(y))) + delta x_6 ],
      y = MW^2.

    Linearizing the equation, rather than inserting Wilson coefficients into
    the iterative input object, gives

      delta y = MZ^2 delta x_6 / (1-d F_SM/dy).

    Evaluate the purely SM derivative numerically at the stored SM solution.
    This retains Jackson's derived-MW convention while ensuring that the
    returned displacement is exactly first order in the Wilson coefficients.
  */
  const double y0 = mw*mw;
  const double stepY = std::max(1.0e-4,1.0e-6*y0);
  if(y0 <= stepY)
    throw std::domain_error("derived MW is too small to linearize");

  const auto smFixedPointMap = [&](double y) {
    inval shifted(smInput);
    shifted.set(InputPar::WMassComplex,std::sqrt(y));
    const double deltaRSM = real(dr_SMNNLO(shifted).result());
    const double radicand = 0.25-k*(1.0+deltaRSM);
    if(radicand <= 0.0)
      throw std::domain_error("SM Delta r gives a non-physical derived MW");
    return mz*mz*(0.5+std::sqrt(radicand));
  };

  const double mapDerivative =
      (smFixedPointMap(y0+stepY)-smFixedPointMap(y0-stepY))
      /(2.0*stepY);
  const double responseDenominator = 1.0-mapDerivative;
  if(std::fabs(responseDenominator) < 1.0e-8)
    throw std::domain_error("derived-MW fixed point is singular");
  const double responseFactor = 1.0/responseDenominator;
  const double finalDeltaX = responseFactor*inputDeltaX;

  DerivedMWShiftD6 shift;
  shift.referenceMW = mw;
  shift.inputDeltaMW2OverMZ2 = inputDeltaX;
  shift.responseFactor = responseFactor;
  shift.deltaMW2OverMZ2 = finalDeltaX;
  shift.deltaMW2 = mz*mz*finalDeltaX;
  shift.deltaMW = shift.deltaMW2/(2.0*mw);
  return shift;
}

} // namespace smeft
} // namespace griffin
