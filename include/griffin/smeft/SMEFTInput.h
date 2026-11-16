/*
  File: SMEFTInput.h
  Purpose: Declare typed Wilson-coefficient inputs, scale conventions and their
           private storage and particle-to-generation lookup.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

/*
  SMEFTInput.h: scheme-independent SMEFT input parameters

  This module stores only the parameters of

    L_SMEFT = L_SM + Sum_i C_i Q_i / Lambda^2,

  with dimensionless Wilson coefficients C_i and Lambda expressed in GeV.
  It deliberately does not inherit from inval/SMvalGmu and does not compute
  MW, electroweak input-scheme shifts, form factors, or observables. Those
  calculations consume an Input object together with a separate SM input.

  The public API follows the coefficient names and flavor-index conventions
  documented in docs/SMEFT_IMPLEMENTATION.md. All coefficients default to zero;
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

// Defined only in the vertex/contact implementation files, not user interfaces.
namespace detail {
struct ChiralWilsonCoefficients;
struct ChiralContact;
}

/*
  Scheme-independent SMEFT input container.

  Named accessors are intentional: calls such as setCphil3(...) and setCll(...)
  make the physics coefficient visible at the call site and avoid exposing the
  internal flat-array layout. All 3x3 and 3^4 flavor labels are accepted.
  Operator-identical Cll/Cee permutations share one coefficient, as specified
  below; setters replace that value rather than adding another contribution.
  No extra flavor-symmetry or Hermiticity assumptions are imposed. Other
  tensors retain their independent entries (apart from the Ceq/Cqe alias).
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

    Cll[p,r,s,t] = Cll[s,t,p,r] by exchange of the identical currents.
    One stored value is the total coefficient of that unique operator, not
    a raw entry to sum over equivalent labels. Setters overwrite; the latest
    assignment wins, including zero. Cll[1,1,2,2] and Cll[1,2,2,1] remain
    independent: for weak-doublet fields the Fierz rearrangement also changes
    the weak-index contraction, so it does not identify these two Cll terms
    as it does for the weak-singlet Cee terms.
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
    Cee[p,r,s,t], Cee[s,t,p,r], Cee[p,t,s,r], Cee[s,r,p,t] share one slot
    by current exchange and the four-dimensional chiral Fierz identity.
    This holds for all generations, not only electron-pair contractions.
    As for Cll, this is one total coefficient: aliases overwrite, never add.
    Importing independently populated redundant raw entries requires summing
    them explicitly before one setter call; see docs/SMEFT_IMPLEMENTATION.md.
  */
  void setCee(Generation p, Generation r, Generation s, Generation t,
              double value);
  double Cee(Generation p, Generation r, Generation s, Generation t) const;

  // Reset every Wilson coefficient to zero without changing Lambda.
  void clearWilsonCoefficients();

private:
  // Internal particle lookup used to select coefficient generation indices.
  // It neither changes coefficients nor calculates observables. Only the
  // vertex/contact matching helpers may use it; it is not a public Input API.
  friend struct detail::ChiralWilsonCoefficients;
  friend struct detail::ChiralContact;

  enum class FermionCategory {
    neutrino,
    chargedLepton,
    upQuark,
    downQuark
  };
  struct FermionInfo {
    Generation generation;
    FermionCategory category;
  };
  // Maps the existing GRIFFIN particle identifiers; rejects unsupported ones.
  // The shared mapping is implemented in SMEFTInput.cc.
  static FermionInfo lookupFermion(int type);

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
  // Canonical representatives of exact operator-identical permutations.
  // Every supplied index is validated before any storage is accessed.
  static std::size_t cllIndex(Generation p, Generation r,
                             Generation s, Generation t);
  static std::size_t ceeIndex(Generation p, Generation r,
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

  // Four-index backing arrays; Cll/Cee use canonical representative slots.
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
