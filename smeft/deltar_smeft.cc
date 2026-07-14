/*
  deltar_smeft.cc: linear dimension-six contribution to Delta r

  With dimensionless Wilson coefficients defined by C_i Q_i/Lambda^2,

    Delta r_D6 = epsilonD6 [
        2 cw/sw CphiWB + cw^2/(2 sw^2) CphiD + Delta_mu ].

  In GmuAlphaMZ, this raw input-transformation coefficient uses the tree-level
  relation

    sw^2 cw^2 = pi alpha/(sqrt(2) G_mu M_Z^2),

  rather than an independent M_W.  derivedMWShiftD6() converts it into the
  final loop-improved derived-MW displacement, including the linear response
  of the SM Delta-r fixed point.  In AlphaMWMZ, the weak angle is defined
  directly by cw=M_W/M_Z.

  The relation follows from Eqs. (2.11), (2.20), and the tree-level muon-decay
  matching in arXiv:2305.03763, translated to dimensionless coefficients.
*/

#include "smeft/deltar_smeft.h"

#include <cmath>
#include <stdexcept>

namespace griffin {
namespace smeft {

Cplx dr_SMEFTLO::resD6Tree() const
{
  return detail::deltaRSMEFTTree(scheme,*ival,*eftInput);
}

Cplx dr_SMEFTLO::result() const
{
  return dr_SMNLO::result()+resD6Tree();
}

} // namespace smeft
} // namespace griffin
