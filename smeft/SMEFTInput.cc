/* SMEFTInput.cc: scheme-independent SMEFT input parameters */

#include "smeft/SMEFTInput.h"

#include <cmath>
#include <stdexcept>

namespace griffin {
namespace smeft {

// Wilson coefficients start at the SM point. Lambda remains deliberately unset.
Input::Input()
  : lambdaGeV_(0.0), lambdaIsSet_(false), CphiWB_(0.0), CphiD_(0.0)
{
  clearWilsonCoefficients();
}

void Input::setLambdaGeV(double value)
{
  if(!std::isfinite(value) || value <= 0.0)
    throw std::invalid_argument("SMEFT Lambda must be finite and positive");

  lambdaGeV_ = value;
  lambdaIsSet_ = true;
}

double Input::lambdaGeV() const
{
  if(!lambdaIsSet_)
    throw std::logic_error("SMEFT Lambda has not been set");

  return lambdaGeV_;
}

double Input::inverseLambdaSquaredGeV() const
{
  const double lambda = lambdaGeV();
  return 1.0/(lambda*lambda);
}

std::size_t Input::generationIndex(Generation generation)
{
  // Physics notation is one-based; std::array storage is zero-based.
  const int value = static_cast<int>(generation);
  if(value < 1 || value > 3)
    throw std::out_of_range("SMEFT generation must be 1, 2, or 3");

  return static_cast<std::size_t>(value - 1);
}

std::size_t Input::matrixIndex(Generation i, Generation j)
{
  // Row-major layout for the 3x3 matrix: [i][j] -> 3*i+j.
  return 3*generationIndex(i) + generationIndex(j);
}

std::size_t Input::tensorIndex(Generation p, Generation r,
                               Generation s, Generation t)
{
  // Row-major layout for a 3x3x3x3 tensor.
  return 27*generationIndex(p) + 9*generationIndex(r)
       + 3*generationIndex(s) + generationIndex(t);
}

void Input::validateWilsonValue(double value)
{
  if(!std::isfinite(value))
    throw std::invalid_argument("SMEFT Wilson coefficient must be finite");
}

void Input::setMatrix(Matrix3& coefficient, Generation i, Generation j,
                      double value)
{
  validateWilsonValue(value);
  coefficient[matrixIndex(i,j)] = value;
}

double Input::getMatrix(const Matrix3& coefficient, Generation i,
                        Generation j)
{
  return coefficient[matrixIndex(i,j)];
}

void Input::setTensor(Tensor4& coefficient, Generation p, Generation r,
                      Generation s, Generation t, double value)
{
  validateWilsonValue(value);
  coefficient[tensorIndex(p,r,s,t)] = value;
}

double Input::getTensor(const Tensor4& coefficient, Generation p,
                        Generation r, Generation s, Generation t)
{
  return coefficient[tensorIndex(p,r,s,t)];
}

void Input::setCphiWB(double value)
{
  validateWilsonValue(value);
  CphiWB_ = value;
}

double Input::CphiWB() const
{
  return CphiWB_;
}

void Input::setCphiD(double value)
{
  validateWilsonValue(value);
  CphiD_ = value;
}

double Input::CphiD() const
{
  return CphiD_;
}

void Input::setCphil1(Generation i, Generation j, double value)
{
  setMatrix(Cphil1_,i,j,value);
}

double Input::Cphil1(Generation i, Generation j) const
{
  return getMatrix(Cphil1_,i,j);
}

void Input::setCphil3(Generation i, Generation j, double value)
{
  setMatrix(Cphil3_,i,j,value);
}

double Input::Cphil3(Generation i, Generation j) const
{
  return getMatrix(Cphil3_,i,j);
}

void Input::setCphie(Generation i, Generation j, double value)
{
  setMatrix(Cphie_,i,j,value);
}

double Input::Cphie(Generation i, Generation j) const
{
  return getMatrix(Cphie_,i,j);
}

void Input::setCphiq1(Generation i, Generation j, double value)
{
  setMatrix(Cphiq1_,i,j,value);
}

double Input::Cphiq1(Generation i, Generation j) const
{
  return getMatrix(Cphiq1_,i,j);
}

void Input::setCphiq3(Generation i, Generation j, double value)
{
  setMatrix(Cphiq3_,i,j,value);
}

double Input::Cphiq3(Generation i, Generation j) const
{
  return getMatrix(Cphiq3_,i,j);
}

void Input::setCphiu(Generation i, Generation j, double value)
{
  setMatrix(Cphiu_,i,j,value);
}

double Input::Cphiu(Generation i, Generation j) const
{
  return getMatrix(Cphiu_,i,j);
}

void Input::setCphid(Generation i, Generation j, double value)
{
  setMatrix(Cphid_,i,j,value);
}

double Input::Cphid(Generation i, Generation j) const
{
  return getMatrix(Cphid_,i,j);
}

void Input::setCll(Generation p, Generation r, Generation s, Generation t,
                   double value)
{
  setTensor(Cll_,p,r,s,t,value);
}

double Input::Cll(Generation p, Generation r, Generation s,
                  Generation t) const
{
  return getTensor(Cll_,p,r,s,t);
}

void Input::setClq1(Generation p, Generation r, Generation s, Generation t,
                    double value)
{
  setTensor(Clq1_,p,r,s,t,value);
}

double Input::Clq1(Generation p, Generation r, Generation s,
                   Generation t) const
{
  return getTensor(Clq1_,p,r,s,t);
}

void Input::setClq3(Generation p, Generation r, Generation s, Generation t,
                    double value)
{
  setTensor(Clq3_,p,r,s,t,value);
}

double Input::Clq3(Generation p, Generation r, Generation s,
                   Generation t) const
{
  return getTensor(Clq3_,p,r,s,t);
}

void Input::setClu(Generation p, Generation r, Generation s, Generation t,
                   double value)
{
  setTensor(Clu_,p,r,s,t,value);
}

double Input::Clu(Generation p, Generation r, Generation s,
                  Generation t) const
{
  return getTensor(Clu_,p,r,s,t);
}

void Input::setCld(Generation p, Generation r, Generation s, Generation t,
                   double value)
{
  setTensor(Cld_,p,r,s,t,value);
}

double Input::Cld(Generation p, Generation r, Generation s,
                  Generation t) const
{
  return getTensor(Cld_,p,r,s,t);
}

void Input::setCqe(Generation p, Generation r, Generation s, Generation t,
                   double value)
{
  setTensor(Cqe_,p,r,s,t,value);
}

double Input::Cqe(Generation p, Generation r, Generation s,
                  Generation t) const
{
  return getTensor(Cqe_,p,r,s,t);
}

void Input::setCeq(Generation p, Generation r, Generation s, Generation t,
                   double value)
{
  // Ceq[p,r,s,t] is the legacy current-reordered spelling of Cqe[s,t,p,r].
  setCqe(s,t,p,r,value);
}

double Input::Ceq(Generation p, Generation r, Generation s,
                  Generation t) const
{
  return Cqe(s,t,p,r);
}

void Input::setCeu(Generation p, Generation r, Generation s, Generation t,
                   double value)
{
  setTensor(Ceu_,p,r,s,t,value);
}

double Input::Ceu(Generation p, Generation r, Generation s,
                  Generation t) const
{
  return getTensor(Ceu_,p,r,s,t);
}

void Input::setCed(Generation p, Generation r, Generation s, Generation t,
                   double value)
{
  setTensor(Ced_,p,r,s,t,value);
}

double Input::Ced(Generation p, Generation r, Generation s,
                  Generation t) const
{
  return getTensor(Ced_,p,r,s,t);
}

void Input::setCle(Generation p, Generation r, Generation s, Generation t,
                   double value)
{
  setTensor(Cle_,p,r,s,t,value);
}

double Input::Cle(Generation p, Generation r, Generation s,
                  Generation t) const
{
  return getTensor(Cle_,p,r,s,t);
}

void Input::setCee(Generation p, Generation r, Generation s, Generation t,
                   double value)
{
  // In the four-dimensional tree convention used for the massless S
  // coefficient, Cee[1,i,i,1] is the same operator as Cee[1,1,i,i].
  // Map both public spellings to one slot so they can never be double counted.
  if(p == Generation::first && t == Generation::first && r == s)
  {
    setTensor(Cee_,Generation::first,Generation::first,r,r,value);
    return;
  }
  setTensor(Cee_,p,r,s,t,value);
}

double Input::Cee(Generation p, Generation r, Generation s,
                  Generation t) const
{
  if(p == Generation::first && t == Generation::first && r == s)
    return getTensor(Cee_,Generation::first,Generation::first,r,r);
  return getTensor(Cee_,p,r,s,t);
}

void Input::clearWilsonCoefficients()
{
  // Lambda is a separate model scale and is intentionally preserved.
  CphiWB_ = 0.0;
  CphiD_ = 0.0;

  Cphil1_.fill(0.0);
  Cphil3_.fill(0.0);
  Cphie_.fill(0.0);
  Cphiq1_.fill(0.0);
  Cphiq3_.fill(0.0);
  Cphiu_.fill(0.0);
  Cphid_.fill(0.0);

  Cll_.fill(0.0);
  Clq1_.fill(0.0);
  Clq3_.fill(0.0);
  Clu_.fill(0.0);
  Cld_.fill(0.0);
  Cqe_.fill(0.0);
  Ceu_.fill(0.0);
  Ced_.fill(0.0);
  Cle_.fill(0.0);
  Cee_.fill(0.0);
}

} // namespace smeft
} // namespace griffin
