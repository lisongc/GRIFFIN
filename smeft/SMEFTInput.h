/*
  SMEFTInput.h: scheme-independent SMEFT input parameters

  This module stores only the parameters of

    L_SMEFT = L_SM + Sum_i C_i Q_i / Lambda^2,

  with dimensionless Wilson coefficients C_i and Lambda expressed in GeV.
  It deliberately does not inherit from inval/SMvalGmu and does not compute
  MW, electroweak input-scheme shifts, form factors, or observables. Those
  calculations consume an Input object together with a separate SM input.

  The public API follows the coefficient names and flavor-index conventions
  documented in SMEFT_IMPLEMENTATION.md. All coefficients default to zero;
  Lambda has no default and must be set explicitly before it is used.
*/

#ifndef GRIFFIN_SMEFT_INPUT_H
#define GRIFFIN_SMEFT_INPUT_H

#include <array>
#include <cstddef>

namespace griffin {
namespace smeft {

// Public generation labels follow the one-based Warsaw-basis convention.
enum class Generation {
  first = 1,
  second = 2,
  third = 3
};

/*
  Scheme-independent SMEFT input container.

  Named accessors are intentional: calls such as setCphil3(...) and setCll(...)
  make the physics coefficient visible at the call site and avoid exposing the
  internal flat-array layout. Two-index coefficients are stored for all 3x3
  flavor entries and four-fermion coefficients for all 3^4 entries. No flavor
  symmetry, Hermiticity relation, or coefficient equality is imposed here.
*/
class Input {
public:
  Input();

  // Lambda is mandatory, uses GeV, and must be finite and positive.
  void setLambdaGeV(double value);
  double lambdaGeV() const;
  // Returns 1/Lambda^2 in GeV^-2 and also checks that Lambda was set.
  double inverseLambdaSquaredGeV() const;

  // Bosonic coefficients with no generation indices.
  void setCphiWB(double value);
  double CphiWB() const;

  void setCphiD(double value);
  double CphiD() const;

  // Two-index current coefficients C[i,j], with i,j in {1,2,3}.
  void setCphil1(Generation i, Generation j, double value);
  double Cphil1(Generation i, Generation j) const;

  void setCphil3(Generation i, Generation j, double value);
  double Cphil3(Generation i, Generation j) const;

  void setCphie(Generation i, Generation j, double value);
  double Cphie(Generation i, Generation j) const;

  void setCphiq1(Generation i, Generation j, double value);
  double Cphiq1(Generation i, Generation j) const;

  void setCphiq3(Generation i, Generation j, double value);
  double Cphiq3(Generation i, Generation j) const;

  void setCphiu(Generation i, Generation j, double value);
  double Cphiu(Generation i, Generation j) const;

  void setCphid(Generation i, Generation j, double value);
  double Cphid(Generation i, Generation j) const;

  /*
    Four-fermion coefficients C[p,r,s,t]. For an operator written as
    (bar psi_p ... psi_r)(bar chi_s ... chi_t), p,r label the first
    bilinear and s,t the second.
  */
  void setCll(Generation p, Generation r, Generation s, Generation t,
              double value);
  double Cll(Generation p, Generation r, Generation s, Generation t) const;

  void setClq1(Generation p, Generation r, Generation s, Generation t,
               double value);
  double Clq1(Generation p, Generation r, Generation s, Generation t) const;

  void setClq3(Generation p, Generation r, Generation s, Generation t,
               double value);
  double Clq3(Generation p, Generation r, Generation s, Generation t) const;

  void setClu(Generation p, Generation r, Generation s, Generation t,
              double value);
  double Clu(Generation p, Generation r, Generation s, Generation t) const;

  void setCld(Generation p, Generation r, Generation s, Generation t,
              double value);
  double Cld(Generation p, Generation r, Generation s, Generation t) const;

  /*
    Canonical Warsaw/WCxf coefficient

      Qqe[p,r,s,t] = (bar q_p gamma q_r)(bar e_s gamma e_t).

    The legacy Ceq API below denotes the same operator with the two currents
    reversed: Ceq[p,r,s,t] := Cqe[s,t,p,r].  It is an alias, not an
    independent Wilson coefficient.
  */
  void setCqe(Generation p, Generation r, Generation s, Generation t,
              double value);
  double Cqe(Generation p, Generation r, Generation s, Generation t) const;

  // Backward-compatible current-reordered alias for Cqe.
  void setCeq(Generation p, Generation r, Generation s, Generation t,
              double value);
  double Ceq(Generation p, Generation r, Generation s, Generation t) const;

  void setCeu(Generation p, Generation r, Generation s, Generation t,
              double value);
  double Ceu(Generation p, Generation r, Generation s, Generation t) const;

  void setCed(Generation p, Generation r, Generation s, Generation t,
              double value);
  double Ced(Generation p, Generation r, Generation s, Generation t) const;

  void setCle(Generation p, Generation r, Generation s, Generation t,
              double value);
  double Cle(Generation p, Generation r, Generation s, Generation t) const;

  /*
    Cee uses one canonical slot for the electron-pair contractions relevant
    here: Cee[1,i,i,1] is mapped to the Fierz-equivalent Cee[1,1,i,i].
  */
  void setCee(Generation p, Generation r, Generation s, Generation t,
              double value);
  double Cee(Generation p, Generation r, Generation s, Generation t) const;

  // Reset every Wilson coefficient to zero without changing Lambda.
  void clearWilsonCoefficients();

private:
  /*
    Fixed-size flat storage keeps the object allocation-free and compact:
      Matrix3: 3^2 entries for C[i,j]
      Tensor4: 3^4 entries for C[p,r,s,t]
    Users never depend on this representation; access is through named methods.
  */
  typedef std::array<double, 9> Matrix3;
  typedef std::array<double, 81> Tensor4;

  // Convert public one-based generations to validated zero-based array indices.
  static std::size_t generationIndex(Generation generation);
  // Row-major flattening: (i,j) -> 3*i+j after zero-based conversion.
  static std::size_t matrixIndex(Generation i, Generation j);
  // Row-major flattening: (p,r,s,t) -> 27*p+9*r+3*s+t.
  static std::size_t tensorIndex(Generation p, Generation r,
                                 Generation s, Generation t);
  // NaN and infinity are rejected at the input boundary.
  static void validateWilsonValue(double value);

  // Shared storage helpers keep every named coefficient accessor consistent.
  static void setMatrix(Matrix3& coefficient, Generation i, Generation j,
                        double value);
  static double getMatrix(const Matrix3& coefficient, Generation i,
                          Generation j);
  static void setTensor(Tensor4& coefficient, Generation p, Generation r,
                        Generation s, Generation t, double value);
  static double getTensor(const Tensor4& coefficient, Generation p,
                          Generation r, Generation s, Generation t);

  // lambdaIsSet_ distinguishes an omitted scale from any numerical value.
  double lambdaGeV_;
  bool lambdaIsSet_;

  // Generation-independent bosonic coefficients.
  double CphiWB_;
  double CphiD_;

  // Complete two-index flavor matrices.
  Matrix3 Cphil1_;
  Matrix3 Cphil3_;
  Matrix3 Cphie_;
  Matrix3 Cphiq1_;
  Matrix3 Cphiq3_;
  Matrix3 Cphiu_;
  Matrix3 Cphid_;

  // Complete four-index flavor tensors.
  Tensor4 Cll_;
  Tensor4 Clq1_;
  Tensor4 Clq3_;
  Tensor4 Clu_;
  Tensor4 Cld_;
  Tensor4 Cqe_;
  Tensor4 Ceu_;
  Tensor4 Ced_;
  Tensor4 Cle_;
  Tensor4 Cee_;
};

} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_INPUT_H
