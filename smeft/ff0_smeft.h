/*
  ff0_smeft.h: tree-level dimension-six shifts of the Zff form factors

  The functions in this file are the SMEFT counterparts of az0(), vz0(), and
  z0() in ff0.h.  They return only the term linear in dimension-six Wilson
  coefficients.  In particular, azD6Tree() and vzD6Tree() do not include the
  SM tree-level form factors themselves.

  The direct current shift is evaluated in {alpha, MW, MZ}.  In GmuAlphaMZ the
  functions additionally propagate the linear displacement of the MW derived
  from {G_mu, alpha, M_Z}.  This follows the derived-MW organization without
  mutating the SM input or generating products of Wilson coefficients.

  The historical three-argument API is retained and means GmuAlphaMZ.  New
  code that may compare input schemes should pass EWInputScheme explicitly.
*/

#ifndef GRIFFIN_SMEFT_FF0_H
#define GRIFFIN_SMEFT_FF0_H

#include "classes.h"
#include "smeft/EWInputScheme.h"
#include "smeft/SMEFTInput.h"

namespace griffin {
namespace smeft {

// Linear dimension-six tree-level shift of the axial-vector Zff factor.
double azD6Tree(int type, const inval& smInput, const Input& eftInput);
double azD6Tree(int type, const inval& smInput, const Input& eftInput,
                EWInputScheme scheme);

// Linear dimension-six tree-level shift of the vector Zff factor.
double vzD6Tree(int type, const inval& smInput, const Input& eftInput);
double vzD6Tree(int type, const inval& smInput, const Input& eftInput,
                EWInputScheme scheme);

// Select the vector or axial-vector shift using VEC/AXV.
double zD6Tree(int type, int form, const inval& smInput,
               const Input& eftInput);
double zD6Tree(int type, int form, const inval& smInput,
               const Input& eftInput, EWInputScheme scheme);

// Strongly typed overloads complement the legacy integer API.
inline double azD6Tree(Fermion type, const inval& smInput,
                       const Input& eftInput)
{
  return azD6Tree(index(type),smInput,eftInput);
}

inline double azD6Tree(Fermion type, const inval& smInput,
                       const Input& eftInput, EWInputScheme scheme)
{
  return azD6Tree(index(type),smInput,eftInput,scheme);
}

inline double vzD6Tree(Fermion type, const inval& smInput,
                       const Input& eftInput)
{
  return vzD6Tree(index(type),smInput,eftInput);
}

inline double vzD6Tree(Fermion type, const inval& smInput,
                       const Input& eftInput, EWInputScheme scheme)
{
  return vzD6Tree(index(type),smInput,eftInput,scheme);
}

inline double zD6Tree(Fermion type, Current form, const inval& smInput,
                      const Input& eftInput)
{
  return zD6Tree(index(type),index(form),smInput,eftInput);
}

inline double zD6Tree(Fermion type, Current form, const inval& smInput,
                      const Input& eftInput, EWInputScheme scheme)
{
  return zD6Tree(index(type),index(form),smInput,eftInput,scheme);
}

} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_FF0_H
