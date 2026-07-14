/* EWPOZ_smeft.cc: linearized SMEFT contributions to Z-pole observables */

#include "smeft/EWPOZ_smeft.h"

#include <cmath>
#include <stdexcept>

#include "ff0.h"
#include "smeft/ff0_smeft.h"

namespace griffin {
namespace smeft {

Cplx SW_SMEFTLO::resD6Tree() const
{
  const double absQf = std::fabs(Qf[ftyp]);
  if(absQf == 0.0)
    throw std::domain_error("SW is undefined for a neutral fermion");

  const double a0 = az0(ftyp,*ival);
  const double v0 = vz0(ftyp,*ival);
  const double da = azD6Tree(ftyp,*ival,*eftInput,scheme);
  const double dv = vzD6Tree(ftyp,*ival,*eftInput,scheme);

  // Linear term in (1-v/a)/(4|Q|); products of SMEFT shifts are omitted.
  return -(dv/a0-v0*da/(a0*a0))/(4.0*absQf);
}

Cplx SW_SMEFTLO::result() const
{
  // Check the domain before invoking SW_SMLO, whose neutral-fermion ratio
  // would otherwise divide by zero.
  if(std::fabs(Qf[ftyp]) == 0.0)
    throw std::domain_error("SW is undefined for a neutral fermion");
  return SW_SMLO::result()+resD6Tree();
}

Cplx FA_SMEFTLO::resD6Tree() const
{
  // Linear term in (a0 + delta a)^2.
  return 2.0*az0(ftyp,*ival)*azD6Tree(ftyp,*ival,*eftInput,scheme);
}

Cplx FA_SMEFTLO::result() const
{
  return FA_SMLO::result()+resD6Tree();
}

Cplx FV_SMEFTLO::resD6Tree() const
{
  // FV_SMLO is algebraically vz0^2; keep only 2*vz0*delta(vz).
  return 2.0*vz0(ftyp,*ival)*vzD6Tree(ftyp,*ival,*eftInput,scheme);
}

Cplx FV_SMEFTLO::result() const
{
  return FV_SMLO::result()+resD6Tree();
}

} // namespace smeft
} // namespace griffin
