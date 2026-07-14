/*
  SMEFTFlavor.h: shared mapping from GRIFFIN fermions to SMEFT flavor labels

  GRIFFIN's historical fermion identifiers encode particle species, while
  Warsaw-basis Wilson coefficients use a separate generation index.  Keeping
  this conversion in one internal helper prevents form factors and matrix
  elements from developing different flavor conventions.
*/

#ifndef GRIFFIN_SMEFT_FLAVOR_H
#define GRIFFIN_SMEFT_FLAVOR_H

#include <stdexcept>

#include "classes.h"
#include "smeft/SMEFTInput.h"

namespace griffin {
namespace smeft {
namespace detail {

enum class FermionFamily {
  neutrino,
  chargedLepton,
  upQuark,
  downQuark
};

struct FermionInfo {
  Generation generation;
  FermionFamily family;
};

inline FermionInfo fermionInfo(int type)
{
  switch(type)
  {
    case NUE: return {Generation::first,FermionFamily::neutrino};
    case NUM: return {Generation::second,FermionFamily::neutrino};
    case NUT: return {Generation::third,FermionFamily::neutrino};
    case ELE: return {Generation::first,FermionFamily::chargedLepton};
    case MUO: return {Generation::second,FermionFamily::chargedLepton};
    case TAU: return {Generation::third,FermionFamily::chargedLepton};
    case UQU: return {Generation::first,FermionFamily::upQuark};
    case CQU: return {Generation::second,FermionFamily::upQuark};
    case DQU: return {Generation::first,FermionFamily::downQuark};
    case SQU: return {Generation::second,FermionFamily::downQuark};
    case BQU: return {Generation::third,FermionFamily::downQuark};
    default:
      throw std::invalid_argument("unsupported fermion type in SMEFT calculation");
  }
}

} // namespace detail
} // namespace smeft
} // namespace griffin

#endif // GRIFFIN_SMEFT_FLAVOR_H
