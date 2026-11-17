/*
  File: WMassSolution.cc
  Purpose: Solve and store the W mass and associated typed SM input values.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16

  GmuAlphaMZ uses Jackson's finite equation; AlphaMWMZ retains the supplied
  mass. See WMassSolution.h for ownership, user access and solver diagnostics.
  Finite prescription adapted from E. Jackson Wallace's deltarSMEFT.cc and
  SMEFTvalG.h, jackson_smeft revision a25076d; typed-input implementation 2026.
*/

#include "griffin/smeft/WMassSolution.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "griffin/sm/deltar.h"

namespace griffin {
namespace smeft {
namespace {

// Legacy inval::get() terminates the process on missing values.  Validate a
// plain copy through protected storage first, so this public resolver reports
// invalid inputs by exception instead.  No derived compute() is inherited.
class CheckedInput : public inval {
public:
  explicit CheckedInput(const inval& input) : inval(input) {}

  double finite(InputPar parameter, const char* name) const
  {
    const std::size_t slot = static_cast<std::size_t>(index(parameter));
    if(slot >= data.size() || !std::isfinite(data[slot]))
      throw std::invalid_argument(std::string("SMEFT input is missing or nonfinite: ")
                                  + name);
    return data[slot];
  }

  double positive(InputPar parameter, const char* name) const
  {
    const double value = finite(parameter,name);
    if(value <= 0.0)
      throw std::domain_error(std::string("SMEFT input must be positive: ")
                               + name);
    return value;
  }

  double nonnegative(InputPar parameter, const char* name) const
  {
    const double value = finite(parameter,name);
    if(value < 0.0)
      throw std::domain_error(std::string("SMEFT input must be nonnegative: ")
                               + name);
    return value;
  }
};

// Dimension-six muon-decay/input matching that drives Jackson's finite
// fixed-point equation.  This is not a derived-MW displacement or a response
// derivative.  Its additive mass-squared term is applied outside the square
// root below, and the resulting mass is solved without linearization.
double gmuInputMatching(const inval& smInput, const Input& eftInput)
{
  const double alpha = smInput.get(InputPar::alpha);
  const double mz = smInput.get(InputPar::ZMassComplex);
  const double gmu = smInput.get(InputPar::fermiConstant);
  if(alpha <= 0.0 || mz <= 0.0 || gmu <= 0.0)
    throw std::domain_error("Gmu input matching requires positive alpha, MZ and Gmu");

  // Preserve the validated expression's arithmetic, including the order of
  // the two MZ factors rather than reusing a separately rounded squared MZ.
  const double sw2cw2 = Pi*alpha/(std::sqrt(2.0)*gmu*mz*mz);
  if(sw2cw2 <= 0.0 || sw2cw2 >= 0.25)
    throw std::domain_error("Gmu inputs do not define a physical tree weak angle");
  const double root = std::sqrt(0.25-sw2cw2);
  const double cw02 = 0.5+root;
  const double sw02 = 0.5-root;
  const double epsilon = eftInput.inverseLambdaSquaredGeV()
      /(std::sqrt(2.0)*gmu);
  const Generation first = Generation::first;
  const Generation second = Generation::second;
  // Input stores the total coefficient of the unique Cll operator. The
  // 1221 and 2112 getters are aliases, not independent raw entries.
  // Jackson's raw sum maps to this one coefficient; summing both getters
  // here would incorrectly double its muon-decay effect.
  const double muonDecayShift =
      eftInput.Cphil3(first,first) + eftInput.Cphil3(second,second)
      - 0.5*eftInput.Cll(first,second,second,first);
  const double bracket = 2.0*std::sqrt(cw02/sw02)*eftInput.CphiWB()
      + cw02/(2.0*sw02)*eftInput.CphiD()
      + muonDecayShift;
  return epsilon*bracket;
}

void updateWWidthAndRunningParameters(inval& point, double gmu, double alphaS)
{
  const double mass = point.get(InputPar::WMassComplex);
  // Preserve Jackson's arithmetic and width prescription, including powint.
  const double width = 0.3376186*gmu*powint(mass,3)*(1.0+0.2122066*alphaS);
  const double conversion = std::sqrt(1.0+sqr(width/mass));
  const double runningMass = mass*conversion;
  const double runningWidth = width*conversion;
  if(!std::isfinite(width) || width < 0.0
      || !std::isfinite(runningMass) || !std::isfinite(runningWidth))
    throw std::domain_error("SMEFT W-width conversion is nonphysical");

  point.set(InputPar::WWidthComplex,width);
  point.set(InputPar::WMassExperiment,runningMass);
  point.set(InputPar::WWidthExperiment,runningWidth);
}

} // namespace

WMassSolution::WMassSolution(const inval& smInput, const Input& eftInput,
                            EWInputScheme scheme,
                            const WMassSolverOptions& options)
  : smInput_(smInput), additiveShift_(0.0), iterations_(0),
    finalResidualGeV2_(0.0)
{
  if(scheme != EWInputScheme::AlphaMWMZ
      && scheme != EWInputScheme::GmuAlphaMZ)
    throw std::invalid_argument("Unsupported SMEFT electroweak input scheme");
  if(!std::isfinite(options.massSquaredToleranceGeV2)
      || options.massSquaredToleranceGeV2 <= 0.0)
    throw std::invalid_argument("SMEFT mass-squared tolerance must be finite and positive");
  if(options.maxIterations == 0)
    throw std::invalid_argument("SMEFT maximum iteration count must be positive");

  // Lambda must be explicitly supplied even at the zero-Wilson point.  The
  // input container already rejects nonpositive and nonfinite Lambda values.
  if(!std::isfinite(eftInput.inverseLambdaSquaredGeV()))
    throw std::domain_error("SMEFT inverse Lambda squared is nonfinite");

  const CheckedInput checked(smInput_);
  const double alpha = checked.positive(InputPar::alpha,"alpha");
  const double mz = checked.positive(InputPar::ZMassComplex,"pole MZ");
  const double mzSquared = mz*mz;
  if(!std::isfinite(mzSquared))
    throw std::domain_error("SMEFT pole MZ squared is nonfinite");

  if(scheme == EWInputScheme::AlphaMWMZ)
  {
    const double mw = checked.positive(InputPar::WMassComplex,"pole MW");
    if(mw >= mz)
      throw std::domain_error("SMEFT MW/MZ does not define a physical weak angle");
    return;
  }

  // Validate all fields read by the SM Delta-r evaluator before calling it.
  const double gmu = checked.positive(InputPar::fermiConstant,"Gmu");
  const double alphaS = checked.nonnegative(InputPar::alphaS,"alphaS");
  checked.positive(InputPar::HiggsMass,"MH");
  checked.positive(InputPar::topMass,"MT");
  checked.nonnegative(InputPar::bottomMass,"MB");
  checked.finite(InputPar::DeltaAlpha,"DeltaAlpha");
  if(smInput_.getsize() <= static_cast<unsigned>(index(InputPar::WWidthExperiment)))
    throw std::invalid_argument("SMEFT input has no storage for the derived W parameters");

  const double k = Pi*alpha/(std::sqrt(2.0)*gmu*mzSquared);
  if(!std::isfinite(k) || k <= 0.0 || k >= 0.25)
    throw std::domain_error("SMEFT Gmu inputs do not define a physical tree weak angle");
  const double matching = gmuInputMatching(smInput_,eftInput);
  additiveShift_ = -k/(2.0*std::sqrt(0.25-k))*matching;
  if(!std::isfinite(additiveShift_))
    throw std::domain_error("SMEFT additive mass shift is nonfinite");

  // Jackson restarts every finite solve at MW=80 GeV.  The SM part of the
  // fixed-point map is evaluated at the current iterate; its SMEFT shift is
  // additive AFTER the square root and is independent of that iterate.
  double massSquared = 80.0*80.0;
  const dr_SMNNLO smDeltaR(smInput_);
  for(unsigned iteration = 0; iteration < options.maxIterations; ++iteration)
  {
    if(!std::isfinite(massSquared) || massSquared <= 0.0 || massSquared >= mzSquared)
      throw std::domain_error("SMEFT W-mass iterate is outside 0 < MW < MZ");
    smInput_.set(InputPar::WMassComplex,std::sqrt(massSquared));

    const Cplx deltaR = smDeltaR.result();
    if(!std::isfinite(real(deltaR)) || !std::isfinite(imag(deltaR)))
      throw std::domain_error("SMEFT SM Delta r is nonfinite");
    // Use the same decimal constant and operation order as Jackson's solver.
    const double radicand = 0.25-2.221441469079183*alpha/(gmu*mzSquared)
                                  *(1.0+real(deltaR));
    if(!std::isfinite(radicand) || radicand <= 0.0)
      throw std::domain_error("SMEFT W-mass fixed point has a nonphysical radicand");
    const double candidate = mzSquared*(0.5+std::sqrt(radicand)+additiveShift_);
    if(!std::isfinite(candidate) || candidate <= 0.0 || candidate >= mzSquared)
      throw std::domain_error("SMEFT derived W mass is outside 0 < MW < MZ");

    iterations_ = iteration+1;
    finalResidualGeV2_ = std::fabs(candidate-massSquared);
    if(finalResidualGeV2_ <= options.massSquaredToleranceGeV2)
    {
      // Do not store candidate: Jackson's solver exits with the last
      // evaluated W pole mass, and derives its width from that mass.
      updateWWidthAndRunningParameters(smInput_,gmu,alphaS);
      return;
    }
    massSquared = candidate;
  }
  throw std::runtime_error("SMEFT W-mass fixed point did not converge within maxIterations");
}

} // namespace smeft
} // namespace griffin
