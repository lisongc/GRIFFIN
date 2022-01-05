#include <iostream>
//#include <math.h>
#include "ff.h"
//#include "classes.h"
//#include "Cplx.h"
//#include "oneloop.h"
using namespace std;

extern const double Qf[7] = {-1., 0, +0.6666666666666667,
                             -0.3333333333333333, -0.3333333333333333, +0.6666666666666666667, -0.333333333333333333333333};
extern const double I3f[7] = {-0.5, +0.5, +0.5, -0.5, -0.5, 0.5, -0.5};

// tree-level axial-vector Z vertex factors
double az0(int type, const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)),
           mz = ival->get(MZ),
           mw = ival->get(MW),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw);

    switch (type)
    {
    case LEP:
    case DQU:
    case BQU:
    case SQU:
        return (-el / (4 * sw * cw));
    case NEU:
    case UQU:
    case CQU:
        return (+el / (4 * sw * cw));
    }
    return 0;
}

// tree-level vector Z vertex factors
double vz0(int type, const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)),
           mz = ival->get(MZ),
           mw = ival->get(MW),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw);
    return az0(type, input) * (1 - 4 * fabs(Qf[type]) * (1 - mw * mw / (mz * mz)));
}

double z0(int type, int formt, const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)),
           mz = ival->get(MZ),
           mw = ival->get(MW),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw);
    if (formt == VEC)
    {
        return vz0(type, input);
    }
    else
    {
        return az0(type, input);
    }
}

double g0(int type, int formt, const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al));
    if (formt == VEC)
    {
        return -el * Qf[type];
    }
    else
    {
        return 0;
    }
}

//one-loop Z-boson (renormalized) self-energies and derivatives
Cplx SigZ1(const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)),
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz;
    return (Cplx(0, 1) * ((el * el * mzs * (63 - 120 * (sw * sw) + 160 * (sw * sw * sw * sw)) * B0Im(mzs, 0, 0))) / (288. * (cw * cw) * (Pi * Pi) * (sw * sw)));
}
//one-loop gamma-gamma renormalized self energy
Cplx SigG1(const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)),
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz;
    return (-(delal * mzs) + (el * el * (96 * mts - 108 * mws - 31 * mzs)) / (432. * (Pi * Pi)) - (el * el * (2 * mts + mzs) * A0fin(mts)) / (9. * mts * (Pi * Pi)) +
            (el * el * (4 * mws + 3 * mzs) * A0fin(mws)) / (16. * mws * (Pi * Pi)) + (Cplx(0, 0.5555555555555556) * (el * el) * mzs * B0Im(mzs, 0, 0)) / (Pi * Pi) +
            (el * el * (2 * mts + mzs) * B0Refin(mzs, mts, mts)) / (9. * (Pi * Pi)) - (el * el * (4 * mws + 3 * mzs) * B0Refin(mzs, mws, mws)) / (16. * (Pi * Pi)));
}
Cplx SigZ1p(const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = 0.303,
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz;
    return (Cplx(0, 1) * ((el * el * (63. - 120. * (sw * sw) + 160 * (sw * sw * sw * sw)) * B0Im(mzs, 0, 0))) / (288. * (cw * cw) * (Pi * Pi) * (sw * sw)));
}
Cplx SigZ1p2(const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = 0.303,
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz;
    return ((el * el * (9 * (cw * cw * cw * cw * cw * cw) * mhs * (mz * mz * mz * mz * mz * mz) * ((mhs - 4 * mzs) * (mhs - 4 * mzs)) * ((-4 * mts + mzs) * (-4 * mts + mzs)) * (-240 * (mws * mws * mws) + 13 * (mz * mz * mz * mz * mz * mz) + 172 * (mws * mws) * mzs - 98 * mws * (mzs * mzs)) + 36 * mws * mzs * ((-4 * mts + mzs) * (-4 * mts + mzs)) * ((-4 * mws + mzs) * (-4 * mws + mzs)) * (-(mhs * mhs * (mz * mz * mz * mz * mz * mz)) + 4 * (mz * mz * mz * mz * mz * mz) * (mzs * mzs) + mhs * (-8 * (mz * mz * mz * mz * mz * mz) * mzs + 12 * (mzs * mzs * mzs * mzs))) + 6 * (cw * cw * cw * cw) * mhs * (mz * mz * mz * mz * mz * mz) * ((mhs - 4 * mzs) * (mhs - 4 * mzs)) * ((-4 * mts + mzs) * (-4 * mts + mzs)) * (48 * (mws * mws * mws) + mz * mz * mz * mz * mz * mz - 20 * (mws * mws) * mzs - 2 * mws * (mzs * mzs)) * (sw * sw) + cw * cw * mhs * (-9 * (mhs * mhs * mhs) * (mz * mz * mz * mz * mz * mz) * ((-4 * mts + mzs) * (-4 * mts + mzs)) * ((-4 * mws + mzs) * (-4 * mws + mzs)) - 4 * mhs * mzs * (48 * (mts * mts * mts) * (mz * mz * mz * mz * mz * mz) * ((-4 * mws + mzs) * (-4 * mws + mzs)) * (-9 - 48 * (sw * sw) + 64 * (sw * sw * sw * sw)) + 8 * (mts * mts) * (136 * (mzs * mzs * mzs * mzs * mzs * mzs) + 2880 * (mws * mws * mws) * (mz * mz * mz * mz * mz * mz) * (sw * sw * sw * sw) + mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * (-589 + 1152 * (sw * sw) - 1548 * (sw * sw * sw * sw)) + 16 * (mws * mws) * mzs * (136 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (-589 + 1152 * (sw * sw) - 1629 * (sw * sw * sw * sw))) + 8 * mws * (mzs * mzs) * (-136 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (589 - 1152 * (sw * sw) + 1575 * (sw * sw * sw * sw)))) - 4 * mts * mzs * (136 * (mzs * mzs * mzs * mzs * mzs * mzs) + 2880 * (mws * mws * mws) * (mz * mz * mz * mz * mz * mz) * (sw * sw * sw * sw) + mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * (-643 + 1152 * (sw * sw) - 1548 * (sw * sw * sw * sw)) + 16 * (mws * mws) * mzs * (136 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (-643 + 1152 * (sw * sw) - 1629 * (sw * sw * sw * sw))) + 8 * mws * (mzs * mzs) * (-136 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (643 - 1152 * (sw * sw) + 1575 * (sw * sw * sw * sw)))) + mzs * mzs * (68 * (mzs * mzs * mzs * mzs * mzs * mzs) + 1440 * (mws * mws * mws) * (mz * mz * mz * mz * mz * mz) * (sw * sw * sw * sw) + mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * (-335 + 576 * (sw * sw) - 774 * (sw * sw * sw * sw)) + 8 * (mws * mws) * mzs * (136 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (-670 + 1152 * (sw * sw) - 1629 * (sw * sw * sw * sw))) + 4 * mws * (mzs * mzs) * (-136 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (670 - 1152 * (sw * sw) + 1575 * (sw * sw * sw * sw))))) + mhs * mhs * (24 * (mts * mts * mts) * (mz * mz * mz * mz * mz * mz) * ((-4 * mws + mzs) * (-4 * mws + mzs)) * (-9 - 48 * (sw * sw) + 64 * (sw * sw * sw * sw)) + 4 * (mts * mts) * (640 * (mzs * mzs * mzs * mzs * mzs * mzs) + 2880 * (mws * mws * mws) * (mz * mz * mz * mz * mz * mz) * (sw * sw * sw * sw) + mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * (-895 + 1152 * (sw * sw) - 1548 * (sw * sw * sw * sw)) + 16 * (mws * mws) * mzs * (640 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (-895 + 1152 * (sw * sw) - 1629 * (sw * sw * sw * sw))) + 8 * mws * (mzs * mzs) * (-640 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (895 - 1152 * (sw * sw) + 1575 * (sw * sw * sw * sw)))) - 2 * mts * mzs * (640 * (mzs * mzs * mzs * mzs * mzs * mzs) + 2880 * (mws * mws * mws) * (mz * mz * mz * mz * mz * mz) * (sw * sw * sw * sw) + mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * (-949 + 1152 * (sw * sw) - 1548 * (sw * sw * sw * sw)) + 16 * (mws * mws) * mzs * (640 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (-949 + 1152 * (sw * sw) - 1629 * (sw * sw * sw * sw))) + 8 * mws * (mzs * mzs) * (-640 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (949 - 1152 * (sw * sw) + 1575 * (sw * sw * sw * sw)))) + mzs * mzs * (160 * (mzs * mzs * mzs * mzs * mzs * mzs) + 720 * (mws * mws * mws) * (mz * mz * mz * mz * mz * mz) * (sw * sw * sw * sw) + mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * (-244 + 288 * (sw * sw) - 387 * (sw * sw * sw * sw)) + 4 * (mws * mws) * mzs * (640 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (-976 + 1152 * (sw * sw) - 1629 * (sw * sw * sw * sw))) + 2 * mws * (mzs * mzs) * (-640 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (976 - 1152 * (sw * sw) + 1575 * (sw * sw * sw * sw))))) - 16 * (mzs * mzs) * (-24 * (mts * mts * mts) * (mz * mz * mz * mz * mz * mz) * ((-4 * mws + mzs) * (-4 * mws + mzs)) * (-9 - 48 * (sw * sw) + 64 * (sw * sw * sw * sw)) - 4 * (mts * mts) * (2880 * (mws * mws * mws) * (mz * mz * mz * mz * mz * mz) * (sw * sw * sw * sw) - mz * mz * mz * mz * mz * mz * (32 * (mz * mz * mz * mz * mz * mz) + mz * mz * mz * mz * mz * mz * (475 - 1152 * (sw * sw) + 1548 * (sw * sw * sw * sw))) + 8 * mws * (32 * (mzs * mzs * mzs * mzs * mzs) + mz * mz * mz * mz * mz * mz * (mzs * mzs) * (475 - 1152 * (sw * sw) + 1575 * (sw * sw * sw * sw))) - 16 * (mws * mws) * (32 * (mzs * mzs * mzs * mzs) + mz * mz * mz * mz * mz * mz * mzs * (475 - 1152 * (sw * sw) + 1629 * (sw * sw * sw * sw)))) - 2 * mts * mzs * (32 * (mzs * mzs * mzs * mzs * mzs * mzs) - 2880 * (mws * mws * mws) * (mz * mz * mz * mz * mz * mz) * (sw * sw * sw * sw) + mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * (529 - 1152 * (sw * sw) + 1548 * (sw * sw * sw * sw)) - 8 * mws * (32 * (mzs * mzs * mzs * mzs * mzs) + mz * mz * mz * mz * mz * mz * (mzs * mzs) * (529 - 1152 * (sw * sw) + 1575 * (sw * sw * sw * sw))) + 16 * (mws * mws) * (32 * (mzs * mzs * mzs * mzs) + mz * mz * mz * mz * mz * mz * mzs * (529 - 1152 * (sw * sw) + 1629 * (sw * sw * sw * sw)))) + mzs * mzs * (8 * (mzs * mzs * mzs * mzs * mzs * mzs) - 720 * (mws * mws * mws) * (mz * mz * mz * mz * mz * mz) * (sw * sw * sw * sw) + mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * (139 - 288 * (sw * sw) + 387 * (sw * sw * sw * sw)) - 2 * mws * (32 * (mzs * mzs * mzs * mzs * mzs) + mz * mz * mz * mz * mz * mz * (mzs * mzs) * (556 - 1152 * (sw * sw) + 1575 * (sw * sw * sw * sw))) + 4 * (mws * mws) * (32 * (mzs * mzs * mzs * mzs) + mz * mz * mz * mz * mz * mz * mzs * (556 - 1152 * (sw * sw) + 1629 * (sw * sw * sw * sw)))))))) /
                (576. * (cw * cw * cw * cw) * mhs * (mz * mz * mz * mz * mz * mz) * ((mhs - 4 * mzs) * (mhs - 4 * mzs)) * (mzs * mzs) * ((-4 * mts + mzs) * (-4 * mts + mzs)) *
                 ((-4 * mws + mzs) * (-4 * mws + mzs)) * (Pi * Pi) * (sw * sw)) -
            (el * el * (-12 * mws * mzs * (-2 * (mhs * mhs) * (mz * mz * mz * mz * mz * mz) + 12 * mhs * (mz * mz * mz * mz * mz * mz) * mzs + 10 * (mz * mz * mz * mz * mz * mz) * (mzs * mzs) - 24 * (mzs * mzs * mzs * mzs * mzs)) + cw * cw * mhs * (6 * (mhs * mhs * mhs) * (mz * mz * mz * mz * mz * mz) - 120 * (mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz) - 54 * (mhs * mhs) * (mz * mz * mz * mz * mz * mz) * mzs + mhs * (-26 * (mz * mz * mz * mz * mz * mz) * (mzs * mzs) + 176 * (mzs * mzs * mzs * mzs * mzs)))) * A0fin(mhs)) /
                (192. * (cw * cw * cw * cw) * mhs * (mz * mz * mz * mz * mz * mz) * ((mhs - 4 * mzs) * (mhs - 4 * mzs)) * (mzs * mzs * mzs) * (Pi * Pi) * (sw * sw)) -
            (el * el * mts * (8 * mzs * (sw * sw) * (3 - 4 * (sw * sw)) + mts * (-9 - 48 * (sw * sw) + 64 * (sw * sw * sw * sw))) * A0fin(mts)) /
                (24. * (cw * cw) * (mzs * mzs) * ((-4 * mts + mzs) * (-4 * mts + mzs)) * (Pi * Pi) * (sw * sw)) +
            (el * el * mws * (cw * cw * cw * cw * (60 * mws - 33 * mzs) + 2 * (cw * cw) * (-4 * mws + mzs) * (sw * sw) + (-20 * mws + 7 * mzs) * (sw * sw * sw * sw)) * A0fin(mws)) /
                (16. * (cw * cw) * (mzs * mzs) * ((-4 * mws + mzs) * (-4 * mws + mzs)) * (Pi * Pi) * (sw * sw)) +
            (el * el * (-12 * mws * mzs * (-2 * (mhs * mhs) * (mz * mz * mz * mz * mz * mz) + 10 * mhs * (mz * mz * mz * mz * mz * mz) * mzs - 4 * (mz * mz * mz * mz * mz * mz) * (mzs * mzs)) + cw * cw * mhs * (6 * (mhs * mhs * mhs) * (mz * mz * mz * mz * mz * mz) - 80 * (mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz) + 32 * (mzs * mzs * mzs * mzs * mzs * mzs) - 2 * (mhs * mhs) * (5 * (mz * mz * mz * mz * mz * mz) * mzs + 19 * (mzs * mzs * mzs * mzs)) + mhs * (44 * (mz * mz * mz * mz * mz * mz) * (mzs * mzs) + 64 * (mzs * mzs * mzs * mzs * mzs)))) * A0fin(mzs)) /
                (192. * (cw * cw * cw * cw) * mhs * (mz * mz * mz * mz * mz * mz) * ((mhs - 4 * mzs) * (mhs - 4 * mzs)) * (mzs * mzs * mzs) * (Pi * Pi) * (sw * sw)) +
            (el * el * (-12 * mws * mzs * (-2 * (mhs * mhs * mhs) * (mz * mz * mz * mz * mz * mz) + 4 * (mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz) + mhs * mhs * (-23 * (mz * mz * mz * mz * mz * mz) * mzs + 37 * (mzs * mzs * mzs * mzs)) + mhs * (44 * (mz * mz * mz * mz * mz * mz) * (mzs * mzs) - 68 * (mzs * mzs * mzs * mzs * mzs))) + cw * cw * mhs * (6 * (mhs * mhs * mhs * mhs) * (mz * mz * mz * mz * mz * mz) - 60 * (mhs * mhs * mhs) * (mz * mz * mz * mz * mz * mz) * mzs + 198 * (mhs * mhs) * (mz * mz * mz * mz * mz * mz) * (mzs * mzs) + 48 * (mz * mz * mz * mz * mz * mz) * (mzs * mzs * mzs * mzs) + 4 * mhs * (11 * (mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz * mz) - 68 * (mzs * mzs * mzs * mzs * mzs * mzs)))) * B0Refin(mzs, mhs, mzs)) /
                (192. * (cw * cw * cw * cw) * mhs * (mz * mz * mz * mz * mz * mz) * ((mhs - 4 * mzs) * (mhs - 4 * mzs)) * (mzs * mzs * mzs) * (Pi * Pi) * (sw * sw)) +
            (el * el * (mts * mts) * (8 * mzs * (sw * sw) * (3 - 4 * (sw * sw)) + mts * (-9 - 48 * (sw * sw) + 64 * (sw * sw * sw * sw))) * B0Refin(mzs, mts, mts)) /
                (24. * (cw * cw) * (mzs * mzs) * ((-4 * mts + mzs) * (-4 * mts + mzs)) * (Pi * Pi) * (sw * sw)) -
            (el * el * (mws * mws) * (cw * cw * cw * cw * (60 * mws - 33 * mzs) + 2 * (cw * cw) * (-4 * mws + mzs) * (sw * sw) + (-20 * mws + 7 * mzs) * (sw * sw * sw * sw)) *
             B0Refin(mzs, mws, mws)) /
                (16. * (cw * cw) * (mzs * mzs) * ((-4 * mws + mzs) * (-4 * mws + mzs)) * (Pi * Pi) * (sw * sw)));

} //one-loop Z, Z', and G vertex form factors
// Z, G, Z' factors at one-loop SM
Cplx Z1(const int ftyp, const int formt, const inval &input)
//Cplx Z_SMNLO::result(void) const
{
    const inval *ival;
    ival = &input;
    // Cplx zf;
    double //el = sqrt(4 * Pi * ival->get(al)),
        el = 0.303,
        mz = ival->get(MZ),
        mw = ival->get(MW),
        mt = ival->get(MT),
        mh = ival->get(MH),
        delal = ival->get(Delal),
        cw = mw / mz,
        sw = sqrt(1 - cw * cw),
        mws = mw * mw,
        mhs = mh * mh,
        mts = mt * mt,
        mzs = mz * mz;
    switch (ftyp)
    {

    case LEP:
        switch (formt)
        {
        case VEC:
            return ((el * el * el * (24 * mhs * mws * (4 * mws - 3 * mzs) * (mzs * mzs) + 24 * mws * (mzs * mzs * mzs) * (-4 * mws + 3 * mzs) + mhs * mhs * mzs * (-44 * (mws * mws) + 41 * mws * mzs - 12 * (mzs * mzs)) + mhs * mhs * mhs * (8 * (mws * mws) - 8 * mws * mzs + 3 * (mzs * mzs))) * A0fin(mhs)) /
                        (1536. * mhs * Power(mws, 2.5) * (mhs - 4 * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (-32 * (mws * mws) * (4 * mws - 3 * mzs) * ((mws - mzs) * (mws - mzs)) * (mzs * mzs) - 36 * (mts * mts * mts) * (2 * mws - 3 * mzs) * (mzs * mzs * mzs) + 2 * mts * mws * mzs * (128 * (mws * mws * mws * mws) - 352 * (mws * mws * mws) * mzs + 332 * (mws * mws) * (mzs * mzs) - 135 * mws * (mzs * mzs * mzs) + 27 * (mzs * mzs * mzs * mzs)) + mts * mts * (2560 * (mws * mws * mws * mws * mws) - 7424 * (mws * mws * mws * mws) * mzs + 7576 * (mws * mws * mws) * (mzs * mzs) - 2946 * (mws * mws) * (mzs * mzs * mzs) + 252 * mws * (mzs * mzs * mzs * mzs) - 27 * (mzs * mzs * mzs * mzs * mzs))) * A0fin(mts)) /
                        (2304. * mts * Power(mws, 2.5) * mzs * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (480 * (mws * mws * mws * mws * mws) - 552 * (mws * mws * mws * mws) * mzs - 340 * (mws * mws * mws) * (mzs * mzs) + 538 * (mws * mws) * (mzs * mzs * mzs) + 2 * mws * (mhs - 53 * mzs) * (mzs * mzs * mzs) - 3 * (mzs * mzs * mzs * mzs) * (mhs + mzs)) * A0fin(mws)) /
                        (1536. * Power(mws, 2.5) * mzs * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs * mws * mzs * (-8 * (mws * mws) + 16 * mws * mzs - 9 * (mzs * mzs)) - mhs * (480 * (mws * mws * mws * mws * mws) - 1560 * (mws * mws * mws * mws) * mzs + 1992 * (mws * mws * mws) * (mzs * mzs) - 1183 * (mws * mws) * (mzs * mzs * mzs) + 247 * mws * (mzs * mzs * mzs * mzs) + 3 * (mzs * mzs * mzs * mzs * mzs)) + 4 * mzs * (480 * (mws * mws * mws * mws * mws) - 1560 * (mws * mws * mws * mws) * mzs + 2012 * (mws * mws * mws) * (mzs * mzs) - 1226 * (mws * mws) * (mzs * mzs * mzs) + 274 * mws * (mzs * mzs * mzs * mzs) + 3 * (mzs * mzs * mzs * mzs * mzs))) * A0fin(mzs)) /
                        (1536. * Power(mws, 2.5) * mzs * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * ((-9 * (el * el) * (mhs * mhs * (4 * mws - 3 * mzs) * (mzs * mzs) + 2 * mhs * (240 * (mws * mws * mws * mws) - 76 * (mws * mws * mws) * mzs - 392 * (mws * mws) * (mzs * mzs) + 293 * mws * (mzs * mzs * mzs) - 6 * (mzs * mzs * mzs * mzs)) + 4 * mzs * (-480 * (mws * mws * mws * mws) + 152 * (mws * mws * mws) * mzs + 784 * (mws * mws) * (mzs * mzs) - 590 * mws * (mzs * mzs * mzs) + 15 * (mzs * mzs * mzs * mzs)))) / (mzs * (-mhs + 4 * mzs)) - (2 * (el * el * (8 * (mzs * mzs) * (320 * (mws * mws * mws) - 632 * (mws * mws) * mzs + 474 * mws * (mzs * mzs) - 135 * (mzs * mzs * mzs)) + mts * mts * (-7680 * (mws * mws * mws) + 14592 * (mws * mws) * mzs - 8136 * mws * (mzs * mzs) + 702 * (mzs * mzs * mzs)) + mts * mzs * (-7168 * (mws * mws * mws) + 14272 * (mws * mws) * mzs - 11604 * mws * (mzs * mzs) + 3807 * (mzs * mzs * mzs))) + 864 * delal * mws * (4 * mts - mzs) * mzs * (4 * (mws * mws) - 7 * mws * mzs + 3 * (mzs * mzs)) * (Pi * Pi))) / (mzs * (-4 * mts + mzs)) +
                           Cplx(0, 960) * (el * el) * mws * (8 * mws - 5 * mzs) * (-mws + mzs) * B0Im(mzs, 0, 0) +
                           Cplx(0, 54) * (el * el) * (80 * (mws * mws * mws) - 184 * (mws * mws) * mzs + 144 * mws * (mzs * mzs) - 45 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (13824. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (3 * (el * el * el) * (2 * mws - 3 * mzs) * (mzs * mzs) * B0Refin(mws, 0, 0)) / (128. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs - 4 * mhs * mws + 12 * (mws * mws)) * (2 * mws - 3 * mzs) * (mzs * mzs) * B0Refin(mws, mhs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (-mts + mws) * (mts + 2 * mws) * (2 * mws - 3 * mzs) * (mzs * mzs) * B0Refin(mws, mts, 0)) / (256. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (2 * mws - 3 * mzs) * (4 * mws - mzs) * mzs * (12 * (mws * mws) + 20 * mws * mzs + mzs * mzs) * B0Refin(mws, mzs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (240 * (mws * mws * mws * mws) - 792 * (mws * mws * mws) * mzs + 1028 * (mws * mws) * (mzs * mzs) - 632 * mws * (mzs * mzs * mzs) + 135 * (mzs * mzs * mzs * mzs)) *
                     B0Refin(mzs, 0, 0)) /
                        (768. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs * mhs * (-8 * (mws * mws) + 16 * mws * mzs - 9 * (mzs * mzs)) + 12 * (mzs * mzs * mzs) * (12 * (mws * mws) - 29 * mws * mzs + 21 * (mzs * mzs)) - 4 * mhs * (mzs * mzs) * (32 * (mws * mws) - 70 * mws * mzs + 45 * (mzs * mzs)) + mhs * mhs * mzs * (52 * (mws * mws) - 107 * mws * mzs + 63 * (mzs * mzs))) *
                     B0Refin(mzs, mhs, mzs)) /
                        (1536. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) -
                    (el * el * el * (mws * (mzs * mzs) * (-128 * (mws * mws * mws) + 352 * (mws * mws) * mzs - 308 * mws * (mzs * mzs) + 75 * (mzs * mzs * mzs)) + 2 * (mts * mts) * (1280 * (mws * mws * mws * mws) - 3712 * (mws * mws * mws) * mzs + 3788 * (mws * mws) * (mzs * mzs) - 1437 * mws * (mzs * mzs * mzs) + 63 * (mzs * mzs * mzs * mzs)) + mts * mzs * (256 * (mws * mws * mws * mws) - 704 * (mws * mws * mws) * mzs + 616 * (mws * mws) * (mzs * mzs) - 204 * mws * (mzs * mzs * mzs) + 81 * (mzs * mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) /
                        (2304. * Power(mws, 1.5) * mzs * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (480 * (mws * mws * mws * mws) - 552 * (mws * mws * mws) * mzs - 412 * (mws * mws) * (mzs * mzs) + 626 * mws * (mzs * mzs * mzs) - 43 * (mzs * mzs * mzs * mzs)) *
                     B0Refin(mzs, mws, mws)) /
                        (1536. * sqrt(mws) * mzs * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (Cplx(0, 0.015625) * (el * el * el) * mzs * ((mws + mzs) * (mws + mzs)) * C0Im(mzs, 0, 0, 0, 0, mws)) / (sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.03125) * (el * el * el) * mzs * (-16 * (mws * mws * mws) + 36 * (mws * mws) * mzs - 30 * mws * (mzs * mzs) + 9 * (mzs * mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * Power(mws, 1.5) * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) /
                        (32. * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * mzs * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) /
                        (64. * sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (16 * (mws * mws * mws) - 36 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) *
                     C0Re(mzs, 0, 0, 0, 0, mzs)) /
                        (32. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        case AXV:
            return (-(el * el * el * mzs * (mhs * mhs * (19 * mws - 4 * mzs) * mzs - 24 * mhs * mws * (mzs * mzs) + 24 * mws * (mzs * mzs * mzs) + mhs * mhs * mhs * (-4 * mws + mzs)) * A0fin(mhs)) /
                        (1536. * mhs * Power(mws, 2.5) * (mhs - 4 * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (-32 * (mws * mws) * ((mws - mzs) * (mws - mzs)) * (mzs * mzs) + 36 * (mts * mts * mts) * (2 * mws - mzs) * (mzs * mzs) + 2 * mts * mws * mzs * (32 * (mws * mws * mws) - 88 * (mws * mws) * mzs + 65 * mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) + mts * mts * (640 * (mws * mws * mws * mws) - 1184 * (mws * mws * mws) * mzs + 622 * (mws * mws) * (mzs * mzs) - 96 * mws * (mzs * mzs * mzs) + 9 * (mzs * mzs * mzs * mzs))) *
                     A0fin(mts)) /
                        (2304. * mts * Power(mws, 2.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (168 * (mws * mws * mws * mws) - 76 * (mws * mws * mws) * mzs - 158 * (mws * mws) * (mzs * mzs) - 2 * mws * (mhs - 41 * mzs) * (mzs * mzs) + mzs * mzs * mzs * (mhs + mzs)) * A0fin(mws)) / (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs * mws * mzs * (-4 * mws + 3 * mzs) + mhs * (-360 * (mws * mws * mws * mws) + 948 * (mws * mws * mws) * mzs - 769 * (mws * mws) * (mzs * mzs) + 201 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs) - 4 * mzs * (-360 * (mws * mws * mws * mws) + 948 * (mws * mws * mws) * mzs - 782 * (mws * mws) * (mzs * mzs) + 210 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs)) * A0fin(mzs)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * ((2 * (el * el) * (6 * (mts * mts) * (320 * (mws * mws) - 272 * mws * mzs + 39 * (mzs * mzs)) - 8 * (mzs * mzs) * (80 * (mws * mws) - 98 * mws * mzs + 45 * (mzs * mzs)) + mts * mzs * (1792 * (mws * mws) - 2368 * mws * mzs + 1269 * (mzs * mzs)))) / (-4 * mts + mzs) +
                           (9 * (el * el) * (mhs * mhs * (mzs * mzs) + 2 * mhs * (60 * (mws * mws * mws) + 108 * (mws * mws) * mzs - 117 * mws * (mzs * mzs) + 8 * (mzs * mzs * mzs)) - 4 * mzs * (120 * (mws * mws * mws) + 216 * (mws * mws) * mzs - 234 * mws * (mzs * mzs) + 17 * (mzs * mzs * mzs)))) / (mhs - 4 * mzs) +
                           1728 * delal * mws * (mws - mzs) * mzs * (Pi * Pi) + Cplx(0, 54) * (el * el) * mzs * (56 * (mws * mws) - 96 * mws * mzs + 35 * (mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (13824. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (3 * (el * el * el) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, 0, 0)) / (128. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (mhs * mhs - 4 * mhs * mws + 12 * (mws * mws)) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mhs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mts - mws) * (mts + 2 * mws) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mts, 0)) / (256. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (2 * mws - mzs) * (4 * mws - mzs) * mzs * (12 * (mws * mws) + 20 * mws * mzs + mzs * mzs) * B0Refin(mws, mzs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (-168 * (mws * mws * mws) + 496 * (mws * mws) * mzs - 412 * mws * (mzs * mzs) + 105 * (mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (768. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (mhs * mhs * (29 * mws - 21 * mzs) * mzs + 12 * (11 * mws - 7 * mzs) * (mzs * mzs * mzs) + mhs * mhs * mhs * (-4 * mws + 3 * mzs) + mhs * (-88 * mws * (mzs * mzs) + 60 * (mzs * mzs * mzs))) * B0Refin(mzs, mhs, mzs)) / (1536. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) + (el * el * el * (mws * (mzs * mzs) * (32 * (mws * mws) - 40 * mws * mzs + 17 * (mzs * mzs)) + mts * mzs * (-64 * (mws * mws * mws) + 80 * (mws * mws) * mzs - 88 * mws * (mzs * mzs) + 27 * (mzs * mzs * mzs)) + mts * mts * (-640 * (mws * mws * mws) + 1184 * (mws * mws) * mzs - 550 * mws * (mzs * mzs) + 42 * (mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) / (2304. * Power(mws, 1.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (168 * (mws * mws * mws) + 44 * (mws * mws) * mzs - 130 * mws * (mzs * mzs) + 17 * (mzs * mzs * mzs)) * B0Refin(mzs, mws, mws)) /
                        (1536. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (Cplx(0, 0.015625) * (el * el * el) * mzs * ((mws + mzs) * (mws + mzs)) * C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (Cplx(0, 0.03125) * (el * el * el) * (mzs * mzs) * (12 * (mws * mws) - 18 * mws * mzs + 7 * (mzs * mzs)) *
                     C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * Power(mws, 1.5) * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * mzs * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) / (64. * sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (mzs * mzs) * (12 * (mws * mws) - 18 * mws * mzs + 7 * (mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (32. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        }
        break;

    case DQU:
    case SQU:

        switch (formt)

        {
        case VEC:
            return ((el * el * el * (24 * mhs * mws * (4 * mws - mzs) * (mzs * mzs) + 24 * mws * (mzs * mzs * mzs) * (-4 * mws + mzs) + mhs * mhs * mzs * (-44 * (mws * mws) + 3 * mws * mzs - 4 * (mzs * mzs)) + mhs * mhs * mhs * (8 * (mws * mws) + mzs * mzs)) * A0fin(mhs)) /
                        (4608. * mhs * Power(mws, 2.5) * (mhs - 4 * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (-32 * (mws * mws) * ((mws - mzs) * (mws - mzs)) * (4 * mws - mzs) * (mzs * mzs) + 36 * (mts * mts * mts) * (mzs * mzs * mzs) * (2 * mws + mzs) + 2 * mts * mws * mzs * (128 * (mws * mws * mws * mws) - 288 * (mws * mws * mws) * mzs + 156 * (mws * mws) * (mzs * mzs) - 5 * mws * (mzs * mzs * mzs) + 9 * (mzs * mzs * mzs * mzs)) + mts * mts * (2560 * (mws * mws * mws * mws * mws) - 6144 * (mws * mws * mws * mws) * mzs + 5208 * (mws * mws * mws) * (mzs * mzs) - 1702 * (mws * mws) * (mzs * mzs * mzs) + 60 * mws * (mzs * mzs * mzs * mzs) - 9 * (mzs * mzs * mzs * mzs * mzs))) * A0fin(mts)) /
                        (6912. * mts * Power(mws, 2.5) * mzs * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (480 * (mws * mws * mws * mws * mws) - 312 * (mws * mws * mws * mws) * mzs - 444 * (mws * mws * mws) * (mzs * mzs) + 414 * (mws * mws) * (mzs * mzs * mzs) - mzs * mzs * mzs * mzs * (mhs + mzs) - 2 * mws * (mzs * mzs * mzs) * (mhs + 43 * mzs)) * A0fin(mws)) /
                        (4608. * Power(mws, 2.5) * mzs * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (-3 * (mhs * mhs) * mws * mzs * (8 * (mws * mws) - 8 * mws * mzs + 3 * (mzs * mzs)) - mhs * (160 * (mws * mws * mws * mws * mws) - 280 * (mws * mws * mws * mws) * mzs + 168 * (mws * mws * mws) * (mzs * mzs) - 295 * (mws * mws) * (mzs * mzs * mzs) + 55 * mws * (mzs * mzs * mzs * mzs) + 3 * (mzs * mzs * mzs * mzs * mzs)) + 4 * mzs * (160 * (mws * mws * mws * mws * mws) - 280 * (mws * mws * mws * mws) * mzs + 228 * (mws * mws * mws) * (mzs * mzs) - 346 * (mws * mws) * (mzs * mzs * mzs) + 82 * mws * (mzs * mzs * mzs * mzs) + 3 * (mzs * mzs * mzs * mzs * mzs))) * A0fin(mzs)) /
                        (13824. * Power(mws, 2.5) * mzs * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * ((3 * (el * el) * (3 * (mhs * mhs) * (mzs * mzs) * (-4 * mws + mzs) - 2 * mhs * (720 * (mws * mws * mws * mws) + 68 * (mws * mws * mws) * mzs - 516 * (mws * mws) * (mzs * mzs) + 255 * mws * (mzs * mzs * mzs) + 4 * (mzs * mzs * mzs * mzs)) + 4 * mzs * (1440 * (mws * mws * mws * mws) + 136 * (mws * mws * mws) * mzs - 1032 * (mws * mws) * (mzs * mzs) + 522 * mws * (mzs * mzs * mzs) + 5 * (mzs * mzs * mzs * mzs)))) / (mzs * (-mhs + 4 * mzs)) + (2 * (el * el * (mts * mzs * (7168 * (mws * mws * mws) - 10688 * (mws * mws) * mzs + 6868 * mws * (mzs * mzs) - 1269 * (mzs * mzs * mzs)) + 6 * (mts * mts) * (1280 * (mws * mws * mws) - 1792 * (mws * mws) * mzs + 812 * mws * (mzs * mzs) - 39 * (mzs * mzs * mzs)) + 8 * (mzs * mzs) * (-320 * (mws * mws * mws) + 472 * (mws * mws) * mzs - 278 * mws * (mzs * mzs) + 45 * (mzs * mzs * mzs))) - 864 * delal * mws * (4 * mts - mzs) * mzs * (4 * (mws * mws) - 5 * mws * mzs + mzs * mzs) * (Pi * Pi))) / (mzs * (-4 * mts + mzs)) +
                           Cplx(0, 960) * (el * el) * mws * (8 * mws - 5 * mzs) * (-mws + mzs) * B0Im(mzs, 0, 0) -
                           Cplx(0, 6) * (el * el) * (4 * mws - mzs) * (16 * (mws * mws) + 64 * mws * mzs - 35 * (mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (41472. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (mzs * mzs) * (2 * mws + mzs) * B0Refin(mws, 0, 0)) / (128. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (mhs * mhs - 4 * mhs * mws + 12 * (mws * mws)) * (mzs * mzs) * (2 * mws + mzs) * B0Refin(mws, mhs, mws)) /
                        (4608. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mts - mws) * (mts + 2 * mws) * (mzs * mzs) * (2 * mws + mzs) * B0Refin(mws, mts, 0)) / (768. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (4 * mws - mzs) * mzs * (2 * mws + mzs) * (12 * (mws * mws) + 20 * mws * mzs + mzs * mzs) * B0Refin(mws, mzs, mws)) /
                        (4608. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (64 * (mws * mws * mws * mws) + 176 * (mws * mws * mws) * mzs - 336 * (mws * mws) * (mzs * mzs) + 320 * mws * (mzs * mzs * mzs) - 35 * (mzs * mzs * mzs * mzs)) *
                     B0Refin(mzs, 0, 0)) /
                        (6912. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs * mhs * (-8 * (mws * mws) + 8 * mws * mzs - 3 * (mzs * mzs)) + 12 * (mzs * mzs * mzs) * (12 * (mws * mws) - 7 * mws * mzs + 7 * (mzs * mzs)) - 4 * mhs * (mzs * mzs) * (32 * (mws * mws) - 26 * mws * mzs + 15 * (mzs * mzs)) + mhs * mhs * mzs * (52 * (mws * mws) - 49 * mws * mzs + 21 * (mzs * mzs))) *
                     B0Refin(mzs, mhs, mzs)) /
                        (4608. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) -
                    (el * el * el * (mws * (mzs * mzs) * (-128 * (mws * mws * mws) + 288 * (mws * mws) * mzs - 228 * mws * (mzs * mzs) + 41 * (mzs * mzs * mzs)) + 2 * (mts * mts) * (1280 * (mws * mws * mws * mws) - 3072 * (mws * mws * mws) * mzs + 2604 * (mws * mws) * (mzs * mzs) - 887 * mws * (mzs * mzs * mzs) + 21 * (mzs * mzs * mzs * mzs)) + mts * mzs * (256 * (mws * mws * mws * mws) - 576 * (mws * mws * mws) * mzs + 456 * (mws * mws) * (mzs * mzs) - 28 * mws * (mzs * mzs * mzs) + 27 * (mzs * mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) / (6912. * Power(mws, 1.5) * mzs * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (160 * (mws * mws * mws * mws) - 72 * (mws * mws * mws) * mzs - 108 * (mws * mws) * (mzs * mzs) + 122 * mws * (mzs * mzs * mzs) - 3 * (mzs * mzs * mzs * mzs)) *
                     B0Refin(mzs, mws, mws)) /
                        (1536. * sqrt(mws) * mzs * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (Cplx(0, 0.005208333333333333) * (el * el * el) * (4 * mws - mzs) * ((mws + mzs) * (mws + mzs)) * C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.0011574074074074073) * (el * el * el) * mzs *
                     (-16 * (mws * mws * mws) + 12 * (mws * mws) * mzs - 30 * mws * (mzs * mzs) + 7 * (mzs * mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * Power(mws, 1.5) * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) /
                        (32. * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (4 * mws - mzs) * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) /
                        (192. * sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * mzs *
                     (16 * (mws * mws * mws) - 12 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) /
                        (864. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;

        case AXV:
            return (-(el * el * el * mzs * (mhs * mhs * (19 * mws - 4 * mzs) * mzs - 24 * mhs * mws * (mzs * mzs) + 24 * mws * (mzs * mzs * mzs) + mhs * mhs * mhs * (-4 * mws + mzs)) * A0fin(mhs)) /
                        (1536. * mhs * Power(mws, 2.5) * (mhs - 4 * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (-32 * (mws * mws) * ((mws - mzs) * (mws - mzs)) * (mzs * mzs) + 36 * (mts * mts * mts) * (2 * mws - mzs) * (mzs * mzs) + 2 * mts * mws * mzs * (32 * (mws * mws * mws) - 88 * (mws * mws) * mzs + 65 * mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) + mts * mts * (640 * (mws * mws * mws * mws) - 1184 * (mws * mws * mws) * mzs + 622 * (mws * mws) * (mzs * mzs) - 96 * mws * (mzs * mzs * mzs) + 9 * (mzs * mzs * mzs * mzs))) *
                     A0fin(mts)) /
                        (2304. * mts * Power(mws, 2.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (136 * (mws * mws * mws * mws) - 60 * (mws * mws * mws) * mzs - 94 * (mws * mws) * (mzs * mzs) - 2 * mws * (mhs - 17 * mzs) * (mzs * mzs) + mzs * mzs * mzs * (mhs + mzs)) * A0fin(mws)) / (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs * mws * mzs * (-4 * mws + 3 * mzs) + mhs * (-40 * (mws * mws * mws * mws) + 108 * (mws * mws * mws) * mzs - 49 * (mws * mws) * (mzs * mzs) + mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs) - 4 * mzs * (-40 * (mws * mws * mws * mws) + 108 * (mws * mws * mws) * mzs - 62 * (mws * mws) * (mzs * mzs) + 10 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs)) * A0fin(mzs)) / (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * ((2 * (el * el) * (6 * (mts * mts) * (320 * (mws * mws) - 272 * mws * mzs + 39 * (mzs * mzs)) - 8 * (mzs * mzs) * (80 * (mws * mws) - 98 * mws * mzs + 45 * (mzs * mzs)) + mts * mzs * (1792 * (mws * mws) - 2368 * mws * mzs + 1269 * (mzs * mzs)))) / (-4 * mts + mzs) +
                           (9 * (el * el) * (mhs * mhs * (mzs * mzs) + 2 * mhs * (60 * (mws * mws * mws) + 88 * (mws * mws) * mzs - 87 * mws * (mzs * mzs) - 2 * (mzs * mzs * mzs)) + 4 * mzs * (-120 * (mws * mws * mws) - 176 * (mws * mws) * mzs + 174 * mws * (mzs * mzs) + 3 * (mzs * mzs * mzs)))) / (mhs - 4 * mzs) +
                           1728 * delal * mws * (mws - mzs) * mzs * (Pi * Pi) - Cplx(0, 18) * (el * el) * (16 * (mws * mws * mws) + 4 * mws * (mzs * mzs) - 5 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (13824. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (3 * (el * el * el) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, 0, 0)) / (128. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (mhs * mhs - 4 * mhs * mws + 12 * (mws * mws)) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mhs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mts - mws) * (mts + 2 * mws) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mts, 0)) / (256. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (2 * mws - mzs) * (4 * mws - mzs) * mzs * (12 * (mws * mws) + 20 * mws * mzs + mzs * mzs) * B0Refin(mws, mzs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (16 * (mws * mws * mws * mws) - 16 * (mws * mws * mws) * mzs + 44 * (mws * mws) * (mzs * mzs) - 28 * mws * (mzs * mzs * mzs) + 5 * (mzs * mzs * mzs * mzs)) *
                     B0Refin(mzs, 0, 0)) /
                        (768. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (mhs * mhs * (29 * mws - 21 * mzs) * mzs + 12 * (11 * mws - 7 * mzs) * (mzs * mzs * mzs) + mhs * mhs * mhs * (-4 * mws + 3 * mzs) + mhs * (-88 * mws * (mzs * mzs) + 60 * (mzs * mzs * mzs))) * B0Refin(mzs, mhs, mzs)) / (1536. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) + (el * el * el * (mws * (mzs * mzs) * (32 * (mws * mws) - 40 * mws * mzs + 17 * (mzs * mzs)) + mts * mzs * (-64 * (mws * mws * mws) + 80 * (mws * mws) * mzs - 88 * mws * (mzs * mzs) + 27 * (mzs * mzs * mzs)) + mts * mts * (-640 * (mws * mws * mws) + 1184 * (mws * mws) * mzs - 550 * mws * (mzs * mzs) + 42 * (mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) / (2304. * Power(mws, 1.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (168 * (mws * mws * mws) + 44 * (mws * mws) * mzs - 130 * mws * (mzs * mzs) + 17 * (mzs * mzs * mzs)) * B0Refin(mzs, mws, mws)) /
                        (1536. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (Cplx(0, 0.005208333333333333) * (el * el * el) * (4 * mws - mzs) * ((mws + mzs) * (mws + mzs)) *
                     C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * (mzs * mzs) * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * Power(mws, 1.5) * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) /
                        (32. * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (4 * mws - mzs) * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) /
                        (192. * sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (mzs * mzs) * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) /
                        (96. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        }
        break;

    case BQU:
        switch (formt)
        {
        case VEC:
            return ((el * el * el * (24 * mhs * mws * (4 * mws - mzs) * (mzs * mzs) + 24 * mws * (mzs * mzs * mzs) * (-4 * mws + mzs) + mhs * mhs * mzs * (-44 * (mws * mws) + 3 * mws * mzs - 4 * (mzs * mzs)) + mhs * mhs * mhs * (8 * (mws * mws) + mzs * mzs)) * A0fin(mhs)) /
                        (4608. * mhs * Power(mws, 2.5) * (mhs - 4 * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * ((9 * mws * (-mws + mzs) * (2 * (mts * mts * mts) * (2 * mws + mzs) + mts * mts * mzs * (2 * mws + mzs) - 6 * mts * mws * (2 * (mws * mws) + 3 * mws * mzs + mzs * mzs) + 4 * (mws * mws) * (2 * (mws * mws) + 5 * mws * mzs + 2 * (mzs * mzs)))) / ((mts - mws) * (mts - mws)) + (-32 * (mws * mws) * ((mws - mzs) * (mws - mzs)) * (4 * mws - mzs) * (mzs * mzs) + 36 * (mts * mts * mts) * (mzs * mzs * mzs) * (2 * mws + mzs) + 2 * mts * mws * mzs * (128 * (mws * mws * mws * mws) - 288 * (mws * mws * mws) * mzs + 156 * (mws * mws) * (mzs * mzs) - 5 * mws * (mzs * mzs * mzs) + 9 * (mzs * mzs * mzs * mzs)) + mts * mts * (2560 * (mws * mws * mws * mws * mws) - 6144 * (mws * mws * mws * mws) * mzs + 5208 * (mws * mws * mws) * (mzs * mzs) - 1702 * (mws * mws) * (mzs * mzs * mzs) + 60 * mws * (mzs * mzs * mzs * mzs) - 9 * (mzs * mzs * mzs * mzs * mzs))) / (mts * mzs * (-4 * mts + mzs))) * A0fin(mts)) /
                        (6912. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (12 * (mts * mts * mts) * mws * mzs * (2 * (mws * mws) - mws * mzs - mzs * mzs) + mts * mts * (480 * (mws * mws * mws * mws * mws) - 360 * (mws * mws * mws * mws) * mzs - 468 * (mws * mws * mws) * (mzs * mzs) + 462 * (mws * mws) * (mzs * mzs * mzs) - mzs * mzs * mzs * mzs * (mhs + mzs) - 2 * mws * (mzs * mzs * mzs) * (mhs + 31 * mzs)) + mws * mws * (480 * (mws * mws * mws * mws * mws) - 312 * (mws * mws * mws * mws) * mzs - 444 * (mws * mws * mws) * (mzs * mzs) + 414 * (mws * mws) * (mzs * mzs * mzs) - mzs * mzs * mzs * mzs * (mhs + mzs) - 2 * mws * (mzs * mzs * mzs) * (mhs + 43 * mzs)) + 2 * mts * mws * (-480 * (mws * mws * mws * mws * mws) + 324 * (mws * mws * mws * mws) * mzs + 480 * (mws * mws * mws) * (mzs * mzs) - 441 * (mws * mws) * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs * (mhs + mzs) + mws * (mzs * mzs * mzs) * (2 * mhs + 65 * mzs))) * A0fin(mws)) /
                        (4608. * ((mts - mws) * (mts - mws)) * Power(mws, 2.5) * mzs * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (-3 * (mhs * mhs) * mws * mzs * (8 * (mws * mws) - 8 * mws * mzs + 3 * (mzs * mzs)) - mhs * (160 * (mws * mws * mws * mws * mws) - 280 * (mws * mws * mws * mws) * mzs + 168 * (mws * mws * mws) * (mzs * mzs) - 295 * (mws * mws) * (mzs * mzs * mzs) + 55 * mws * (mzs * mzs * mzs * mzs) + 3 * (mzs * mzs * mzs * mzs * mzs)) + 4 * mzs * (160 * (mws * mws * mws * mws * mws) - 280 * (mws * mws * mws * mws) * mzs + 228 * (mws * mws * mws) * (mzs * mzs) - 346 * (mws * mws) * (mzs * mzs * mzs) + 82 * mws * (mzs * mzs * mzs * mzs) + 3 * (mzs * mzs * mzs * mzs * mzs))) * A0fin(mzs)) /
                        (13824. * Power(mws, 2.5) * mzs * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * ((3 * (el * el) * (3 * (mhs * mhs) * (mts - mws) * (4 * mws - mzs) * (mzs * mzs) + 4 * mzs * (27 * (mts * mts) * (6 * mws - 5 * mzs) * (mzs * mzs) - mts * (1440 * (mws * mws * mws * mws) + 136 * (mws * mws * mws) * mzs - 762 * (mws * mws) * (mzs * mzs) + 441 * mws * (mzs * mzs * mzs) + 5 * (mzs * mzs * mzs * mzs)) + mws * (1440 * (mws * mws * mws * mws) + 136 * (mws * mws * mws) * mzs - 1032 * (mws * mws) * (mzs * mzs) + 522 * mws * (mzs * mzs * mzs) + 5 * (mzs * mzs * mzs * mzs))) + mhs * (-27 * (mts * mts) * (6 * mws - 5 * mzs) * (mzs * mzs) - 2 * mws * (720 * (mws * mws * mws * mws) + 68 * (mws * mws * mws) * mzs - 516 * (mws * mws) * (mzs * mzs) + 255 * mws * (mzs * mzs * mzs) + 4 * (mzs * mzs * mzs * mzs)) + mts * (1440 * (mws * mws * mws * mws) + 136 * (mws * mws * mws) * mzs - 762 * (mws * mws) * (mzs * mzs) + 429 * mws * (mzs * mzs * mzs) + 8 * (mzs * mzs * mzs * mzs))))) / ((-mts + mws) * mzs * (-mhs + 4 * mzs)) +
                           (2 * (el * el * (mts * mzs * (7168 * (mws * mws * mws) - 10688 * (mws * mws) * mzs + 6868 * mws * (mzs * mzs) - 1269 * (mzs * mzs * mzs)) + 6 * (mts * mts) * (1280 * (mws * mws * mws) - 1792 * (mws * mws) * mzs + 812 * mws * (mzs * mzs) - 39 * (mzs * mzs * mzs)) + 8 * (mzs * mzs) * (-320 * (mws * mws * mws) + 472 * (mws * mws) * mzs - 278 * mws * (mzs * mzs) + 45 * (mzs * mzs * mzs))) -
                                 864 * delal * mws * (4 * mts - mzs) * mzs * (4 * (mws * mws) - 5 * mws * mzs + mzs * mzs) * (Pi * Pi))) /
                               (mzs * (-4 * mts + mzs)) +
                           Cplx(0, 960) * (el * el) * mws * (8 * mws - 5 * mzs) * (-mws + mzs) * B0Im(mzs, 0, 0) +
                           Cplx(0, 30) * (el * el) * (16 * (mws * mws * mws) - 12 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (41472. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (mzs * mzs) * (2 * mws + mzs) * B0Refin(mws, 0, 0)) / (128. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (mhs * mhs - 4 * mhs * mws + 12 * (mws * mws)) * (mzs * mzs) * (2 * mws + mzs) * B0Refin(mws, mhs, mws)) /
                        (4608. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mts - mws) * (mts + 2 * mws) * (mzs * mzs) * (2 * mws + mzs) * B0Refin(mws, mts, 0)) / (768. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (4 * mws - mzs) * mzs * (2 * mws + mzs) * (12 * (mws * mws) + 20 * mws * mzs + mzs * mzs) * B0Refin(mws, mzs, mws)) /
                        (4608. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (80 * (mws * mws * mws * mws) - 140 * (mws * mws * mws) * mzs + 102 * (mws * mws) * (mzs * mzs) - 266 * mws * (mzs * mzs * mzs) + 35 * (mzs * mzs * mzs * mzs)) *
                     B0Refin(mzs, 0, 0)) /
                        (6912. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs * mhs * (-8 * (mws * mws) + 8 * mws * mzs - 3 * (mzs * mzs)) + 12 * (mzs * mzs * mzs) * (12 * (mws * mws) - 7 * mws * mzs + 7 * (mzs * mzs)) - 4 * mhs * (mzs * mzs) * (32 * (mws * mws) - 26 * mws * mzs + 15 * (mzs * mzs)) + mhs * mhs * mzs * (52 * (mws * mws) - 49 * mws * mzs + 21 * (mzs * mzs))) *
                     B0Refin(mzs, mhs, mzs)) /
                        (4608. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * el * el * (288 * (mts * mts * mts) * ((mws - mzs) * (mws - mzs)) * mzs + mws * (mzs * mzs) * (272 * (mws * mws * mws) - 252 * (mws * mws) * mzs - 6 * mws * (mzs * mzs) + 13 * (mzs * mzs * mzs)) + mts * mzs * (-832 * (mws * mws * mws * mws) + 360 * (mws * mws * mws) * mzs + 480 * (mws * mws) * (mzs * mzs) - 80 * mws * (mzs * mzs * mzs) - 63 * (mzs * mzs * mzs * mzs)) + mts * mts * (-2560 * (mws * mws * mws * mws) + 6432 * (mws * mws * mws) * mzs - 5280 * (mws * mws) * (mzs * mzs) + 1486 * mws * (mzs * mzs * mzs) + 30 * (mzs * mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) / (6912. * Power(mws, 1.5) * mzs * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (-160 * (mws * mws * mws * mws * mws) + 72 * (mws * mws * mws * mws) * mzs + 6 * mts * (mzs * mzs * mzs) * (2 * mts + mzs) + 12 * (mws * mws * mws) * mzs * (2 * mts + 9 * mzs) + 3 * mws * (mzs * mzs) * (-12 * (mts * mts) - 10 * mts * mzs + mzs * mzs) + 2 * (mws * mws) * (12 * (mts * mts) * mzs - 61 * (mzs * mzs * mzs))) * B0Refin(mzs, mws, mws)) / (1536. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (Cplx(0, 0.0011574074074074073) * (el * el * el) * mzs * (-16 * (mws * mws * mws) + 12 * (mws * mws) * mzs - 30 * mws * (mzs * mzs) + 7 * (mzs * mzs * mzs)) *
                     C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (mts * mts * mts * (2 * mws - mzs) + mts * mts * (4 * mws - mzs) * mzs + 4 * (mws * mws * mws) * (mws + 2 * mzs) + mts * mws * (-6 * (mws * mws) - 5 * mws * mzs + 4 * (mzs * mzs))) * C0Re(0, mzs, 0, mts, mws, mws)) / (128. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (16 * (mws * mws * mws) - 12 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) /
                        (864. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (4 * (mts * mts * mts) * (mws - mzs) + mts * mts * (10 * mws - mzs) * mzs + 2 * mws * (4 * mws - mzs) * ((mws + mzs) * (mws + mzs)) - 4 * mts * mws * (3 * (mws * mws) + 2 * mws * mzs + mzs * mzs)) * C0Re(mzs, 0, 0, mts, mts, mws)) / (384. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        case AXV:
            return (-(el * el * el * mzs * (mhs * mhs * (19 * mws - 4 * mzs) * mzs - 24 * mhs * mws * (mzs * mzs) + 24 * mws * (mzs * mzs * mzs) + mhs * mhs * mhs * (-4 * mws + mzs)) * A0fin(mhs)) /
                        (1536. * mhs * Power(mws, 2.5) * (mhs - 4 * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * ((3 * mws * (-mws + mzs) * (2 * (mts * mts * mts) * (2 * mws + mzs) + mts * mts * mzs * (2 * mws + mzs) - 6 * mts * mws * (2 * (mws * mws) + 3 * mws * mzs + mzs * mzs) + 4 * (mws * mws) * (2 * (mws * mws) + 5 * mws * mzs + 2 * (mzs * mzs)))) / ((mts - mws) * (mts - mws)) + (-32 * (mws * mws) * ((mws - mzs) * (mws - mzs)) * (mzs * mzs) + 36 * (mts * mts * mts) * (2 * mws - mzs) * (mzs * mzs) + 2 * mts * mws * mzs * (32 * (mws * mws * mws) - 88 * (mws * mws) * mzs + 65 * mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) + mts * mts * (640 * (mws * mws * mws * mws) - 1184 * (mws * mws * mws) * mzs + 622 * (mws * mws) * (mzs * mzs) - 96 * mws * (mzs * mzs * mzs) + 9 * (mzs * mzs * mzs * mzs))) / (mts * (-4 * mts + mzs))) * A0fin(mts)) / (2304. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (120 * (mws * mws * mws * mws) - 160 * (mws * mws * mws) * mzs + 84 * mws * ((mws - mzs) * (mws - mzs)) * mzs + (4 * mws * (2 * (mts * mts) + 3 * mts * mws - 2 * (mws * mws)) * ((mws - mzs) * (mws - mzs)) * mzs) / ((mts - mws) * (mts - mws)) - 2 * mhs * mws * (mzs * mzs) + 94 * (mws * mws) * (mzs * mzs) + (6 * mws * (2 * (mts * mts) + 3 * mts * mws - 2 * (mws * mws)) * (mws - mzs) * (mzs * mzs)) / ((mts - mws) * (mts - mws)) + mhs * (mzs * mzs * mzs) - 38 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs + (4 * mws * (mws - mzs) * (2 * mws + mzs) * (mts * mts + mts * mws - 2 * mws * (mws + 2 * mzs))) / (mts - mws)) * A0fin(mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs * mws * mzs * (-4 * mws + 3 * mzs) + mhs * (-40 * (mws * mws * mws * mws) + 108 * (mws * mws * mws) * mzs - 49 * (mws * mws) * (mzs * mzs) + mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs) - 4 * mzs * (-40 * (mws * mws * mws * mws) + 108 * (mws * mws * mws) * mzs - 62 * (mws * mws) * (mzs * mzs) + 10 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs)) * A0fin(mzs)) / (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * ((2 * (el * el) * (6 * (mts * mts) * (320 * (mws * mws) - 272 * mws * mzs + 39 * (mzs * mzs)) - 8 * (mzs * mzs) * (80 * (mws * mws) - 98 * mws * mzs + 45 * (mzs * mzs)) + mts * mzs * (1792 * (mws * mws) - 2368 * mws * mzs + 1269 * (mzs * mzs)))) / (-4 * mts + mzs) +
                           (9 * (el * el) * (mhs * mhs * (mts - mws) * (mzs * mzs) + mhs * (-3 * (mts * mts) * (6 * mws - 5 * mzs) * mzs + mts * (120 * (mws * mws * mws) + 206 * (mws * mws) * mzs - 183 * mws * (mzs * mzs) - 4 * (mzs * mzs * mzs)) + 2 * mws * (-60 * (mws * mws * mws) - 88 * (mws * mws) * mzs + 87 * mws * (mzs * mzs) + 2 * (mzs * mzs * mzs))) + 4 * mzs * (3 * (mts * mts) * (6 * mws - 5 * mzs) * mzs + mws * (120 * (mws * mws * mws) + 176 * (mws * mws) * mzs - 174 * mws * (mzs * mzs) - 3 * (mzs * mzs * mzs)) + mts * (-120 * (mws * mws * mws) - 206 * (mws * mws) * mzs + 183 * mws * (mzs * mzs) + 3 * (mzs * mzs * mzs))))) / ((mts - mws) * (mhs - 4 * mzs)) +
                           1728 * delal * mws * (mws - mzs) * mzs * (Pi * Pi) + Cplx(0, 90) * (el * el) * mzs * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * B0Im(mzs, 0, 0))) /
                        (13824. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (3 * (el * el * el) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, 0, 0)) / (128. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (mhs * mhs - 4 * mhs * mws + 12 * (mws * mws)) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mhs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mts - mws) * (mts + 2 * mws) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mts, 0)) / (256. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (2 * mws - mzs) * (4 * mws - mzs) * mzs * (12 * (mws * mws) + 20 * mws * mzs + mzs * mzs) * B0Refin(mws, mzs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (-20 * (mws * mws * mws) + 70 * (mws * mws) * mzs - 34 * mws * (mzs * mzs) + 5 * (mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (768. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (mhs * mhs * (29 * mws - 21 * mzs) * mzs + 12 * (11 * mws - 7 * mzs) * (mzs * mzs * mzs) + mhs * mhs * mhs * (-4 * mws + 3 * mzs) + mhs * (-88 * mws * (mzs * mzs) + 60 * (mzs * mzs * mzs))) * B0Refin(mzs, mhs, mzs)) / (1536. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) + (el * el * el * (96 * (mts * mts * mts) * ((mws - mzs) * (mws - mzs)) + mws * mzs * (48 * (mws * mws * mws) + 44 * (mws * mws) * mzs - 118 * mws * (mzs * mzs) + 35 * (mzs * mzs * mzs)) + mts * mts * (-544 * (mws * mws * mws) + 1160 * (mws * mws) * mzs - 646 * mws * (mzs * mzs) + 66 * (mzs * mzs * mzs)) + mts * (-192 * (mws * mws * mws * mws) - 136 * (mws * mws * mws) * mzs + 392 * (mws * mws) * (mzs * mzs) - 124 * mws * (mzs * mzs * mzs) + 15 * (mzs * mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) / (2304. * Power(mws, 1.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (12 * (mts * mts) * (2 * (mws * mws) - 3 * mws * mzs + mzs * mzs) + 6 * mts * (4 * (mws * mws * mws) - 5 * mws * (mzs * mzs) + mzs * mzs * mzs) - mws * (168 * (mws * mws * mws) + 44 * (mws * mws) * mzs - 130 * mws * (mzs * mzs) + 17 * (mzs * mzs * mzs))) * B0Refin(mzs, mws, mws)) /
                        (1536. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * (mzs * mzs) * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el *
                     (mts * mts * mts * (2 * mws - mzs) + mts * mts * (4 * mws - mzs) * mzs + 4 * (mws * mws * mws) * (mws + 2 * mzs) +
                      mts * mws * (-6 * (mws * mws) - 5 * mws * mzs + 4 * (mzs * mzs))) *
                     C0Re(0, mzs, 0, mts, mws, mws)) /
                        (128. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (mzs * mzs) * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (96. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (4 * (mts * mts * mts) * (mws - mzs) + mts * mts * (10 * mws - mzs) * mzs + 2 * mws * (4 * mws - mzs) * ((mws + mzs) * (mws + mzs)) - 4 * mts * mws * (3 * (mws * mws) + 2 * mws * mzs + mzs * mzs)) * C0Re(mzs, 0, 0, mts, mts, mws)) / (384. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        }
        break;
    case UQU:
    case CQU:
        switch (formt)
        {
        case VEC:
            return ((el * el * el * (-24 * mhs * mws * (8 * mws - 5 * mzs) * (mzs * mzs) + 24 * mws * (8 * mws - 5 * mzs) * (mzs * mzs * mzs) + mhs * mhs * mhs * (-16 * (mws * mws) + 12 * mws * mzs - 5 * (mzs * mzs)) + mhs * mhs * mzs * (88 * (mws * mws) - 63 * mws * mzs + 20 * (mzs * mzs))) * A0fin(mhs)) /
                        (4608. * mhs * Power(mws, 2.5) * (mhs - 4 * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (32 * (mws * mws) * (8 * mws - 5 * mzs) * ((mws - mzs) * (mws - mzs)) * (mzs * mzs) + 36 * (mts * mts * mts) * (2 * mws - 5 * mzs) * (mzs * mzs * mzs) - 2 * mts * mws * mzs * (256 * (mws * mws * mws * mws) - 672 * (mws * mws * mws) * mzs + 576 * (mws * mws) * (mzs * mzs) - 205 * mws * (mzs * mzs * mzs) + 45 * (mzs * mzs * mzs * mzs)) + mts * mts * (-5120 * (mws * mws * mws * mws * mws) + 14208 * (mws * mws * mws * mws) * mzs - 13968 * (mws * mws * mws) * (mzs * mzs) + 5270 * (mws * mws) * (mzs * mzs * mzs) - 408 * mws * (mzs * mzs * mzs * mzs) + 45 * (mzs * mzs * mzs * mzs * mzs))) * A0fin(mts)) /
                        (6912. * mts * Power(mws, 2.5) * mzs * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (-960 * (mws * mws * mws * mws * mws) + 984 * (mws * mws * mws * mws) * mzs + 732 * (mws * mws * mws) * (mzs * mzs) - 1014 * (mws * mws) * (mzs * mzs * mzs) - 2 * mws * (mhs - 101 * mzs) * (mzs * mzs * mzs) + 5 * (mzs * mzs * mzs * mzs) * (mhs + mzs)) * A0fin(mws)) /
                        (4608. * Power(mws, 2.5) * mzs * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (3 * (mhs * mhs) * mws * mzs * (16 * (mws * mws) - 28 * mws * mzs + 15 * (mzs * mzs)) + mhs * (1280 * (mws * mws * mws * mws * mws) - 3680 * (mws * mws * mws * mws) * mzs + 4368 * (mws * mws * mws) * (mzs * mzs) - 2711 * (mws * mws) * (mzs * mzs * mzs) + 539 * mws * (mzs * mzs * mzs * mzs) + 15 * (mzs * mzs * mzs * mzs * mzs)) - 4 * mzs * (1280 * (mws * mws * mws * mws * mws) - 3680 * (mws * mws * mws * mws) * mzs + 4488 * (mws * mws * mws) * (mzs * mzs) - 2930 * (mws * mws) * (mzs * mzs * mzs) + 674 * mws * (mzs * mzs * mzs * mzs) + 15 * (mzs * mzs * mzs * mzs * mzs))) * A0fin(mzs)) /
                        (13824. * Power(mws, 2.5) * mzs * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * ((3 * (el * el) * (3 * (mhs * mhs) * (8 * mws - 5 * mzs) * (mzs * mzs) + 4 * mzs * (-2880 * (mws * mws * mws * mws) + 712 * (mws * mws * mws) * mzs + 3792 * (mws * mws) * (mzs * mzs) - 2700 * mws * (mzs * mzs * mzs) + 5 * (mzs * mzs * mzs * mzs)) + 2 * mhs * (1440 * (mws * mws * mws * mws) - 356 * (mws * mws * mws) * mzs - 1896 * (mws * mws) * (mzs * mzs) + 1338 * mws * (mzs * mzs * mzs) + 5 * (mzs * mzs * mzs * mzs)))) / (mzs * (-mhs + 4 * mzs)) +
                           (2 * (el * el * (8 * (mzs * mzs) * (640 * (mws * mws * mws) - 1184 * (mws * mws) * mzs + 850 * mws * (mzs * mzs) - 225 * (mzs * mzs * mzs)) - 6 * (mts * mts) * (2560 * (mws * mws * mws) - 4544 * (mws * mws) * mzs + 2440 * mws * (mzs * mzs) - 195 * (mzs * mzs * mzs)) + mts * mzs * (-14336 * (mws * mws * mws) + 26752 * (mws * mws) * mzs - 20840 * mws * (mzs * mzs) + 6345 * (mzs * mzs * mzs))) +
                                 864 * delal * mws * (4 * mts - mzs) * mzs * (8 * (mws * mws) - 13 * mws * mzs + 5 * (mzs * mzs)) * (Pi * Pi))) /
                               (mzs * (-4 * mts + mzs)) +
                           Cplx(0, 1920) * (el * el) * mws * (8 * mws - 5 * mzs) * (mws - mzs) * B0Im(mzs, 0, 0) -
                           Cplx(0, 6) * (el * el) * (568 * (mws * mws * mws) - 1344 * (mws * mws) * mzs + 966 * mws * (mzs * mzs) - 325 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (41472. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (2 * mws - 5 * mzs) * (mzs * mzs) * B0Refin(mws, 0, 0)) / (128. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (mhs * mhs - 4 * mhs * mws + 12 * (mws * mws)) * (2 * mws - 5 * mzs) * (mzs * mzs) * B0Refin(mws, mhs, mws)) /
                        (4608. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mts - mws) * (mts + 2 * mws) * (2 * mws - 5 * mzs) * (mzs * mzs) * B0Refin(mws, mts, 0)) / (768. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (2 * mws - 5 * mzs) * (4 * mws - mzs) * mzs * (12 * (mws * mws) + 20 * mws * mzs + mzs * mzs) * B0Refin(mws, mzs, mws)) /
                        (4608. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (568 * (mws * mws * mws * mws) - 1912 * (mws * mws * mws) * mzs + 2454 * (mws * mws) * (mzs * mzs) - 1624 * mws * (mzs * mzs * mzs) + 325 * (mzs * mzs * mzs * mzs)) *
                     B0Refin(mzs, 0, 0)) /
                        (6912. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs * mzs * (-104 * (mws * mws) + 185 * mws * mzs - 105 * (mzs * mzs)) + mhs * mhs * mhs * (16 * (mws * mws) - 28 * mws * mzs + 15 * (mzs * mzs)) - 12 * (mzs * mzs * mzs) * (24 * (mws * mws) - 47 * mws * mzs + 35 * (mzs * mzs)) + 4 * mhs * (mzs * mzs) * (64 * (mws * mws) - 118 * mws * mzs + 75 * (mzs * mzs))) *
                     B0Refin(mzs, mhs, mzs)) /
                        (4608. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * el * el * (mws * (mzs * mzs) * (-256 * (mws * mws * mws) + 672 * (mws * mws) * mzs - 576 * mws * (mzs * mzs) + 133 * (mzs * mzs * mzs)) + 2 * (mts * mts) * (2560 * (mws * mws * mws * mws) - 7104 * (mws * mws * mws) * mzs + 6984 * (mws * mws) * (mzs * mzs) - 2599 * mws * (mzs * mzs * mzs) + 105 * (mzs * mzs * mzs * mzs)) + mts * mzs * (512 * (mws * mws * mws * mws) - 1344 * (mws * mws * mws) * mzs + 1152 * (mws * mws) * (mzs * mzs) - 320 * mws * (mzs * mzs * mzs) + 135 * (mzs * mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) /
                        (6912. * Power(mws, 1.5) * mzs * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (320 * (mws * mws * mws * mws) - 312 * (mws * mws * mws) * mzs - 260 * (mws * mws) * (mzs * mzs) + 374 * mws * (mzs * mzs * mzs) - 23 * (mzs * mzs * mzs * mzs)) *
                     B0Refin(mzs, mws, mws)) /
                        (1536. * sqrt(mws) * mzs * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (Cplx(0, 0.005208333333333333) * (el * el * el) * ((mws + mzs) * (mws + mzs)) * (2 * mws + mzs) * C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.0011574074074074073) * (el * el * el) * mzs *
                     (128 * (mws * mws * mws) - 240 * (mws * mws) * mzs + 204 * mws * (mzs * mzs) - 65 * (mzs * mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) /
                        (32. * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * ((mws + mzs) * (mws + mzs)) * (2 * mws + mzs) * C0Re(mzs, 0, 0, 0, 0, mws)) /
                        (192. * sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * mzs *
                     (-128 * (mws * mws * mws) + 240 * (mws * mws) * mzs - 204 * mws * (mzs * mzs) + 65 * (mzs * mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) /
                        (864. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        case AXV:
            return ((el * el * el * mzs * (mhs * mhs * (19 * mws - 4 * mzs) * mzs - 24 * mhs * mws * (mzs * mzs) + 24 * mws * (mzs * mzs * mzs) + mhs * mhs * mhs * (-4 * mws + mzs)) * A0fin(mhs)) /
                        (1536. * mhs * Power(mws, 2.5) * (mhs - 4 * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (32 * (mws * mws) * ((mws - mzs) * (mws - mzs)) * (mzs * mzs) + 36 * (mts * mts * mts) * (mzs * mzs) * (-2 * mws + mzs) - 2 * mts * mws * mzs * (32 * (mws * mws * mws) - 88 * (mws * mws) * mzs + 65 * mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) + mts * mts * (-640 * (mws * mws * mws * mws) + 1184 * (mws * mws * mws) * mzs - 622 * (mws * mws) * (mzs * mzs) + 96 * mws * (mzs * mzs * mzs) - 9 * (mzs * mzs * mzs * mzs))) *
                     A0fin(mts)) /
                        (2304. * mts * Power(mws, 2.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (152 * (mws * mws * mws * mws) - 68 * (mws * mws * mws) * mzs - 126 * (mws * mws) * (mzs * mzs) - 2 * mws * (mhs - 29 * mzs) * (mzs * mzs) + mzs * mzs * mzs * (mhs + mzs)) * A0fin(mws)) / (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs * mws * (4 * mws - 3 * mzs) * mzs + mhs * (160 * (mws * mws * mws * mws) - 408 * (mws * mws * mws) * mzs + 289 * (mws * mws) * (mzs * mzs) - 61 * mws * (mzs * mzs * mzs) - mzs * mzs * mzs * mzs) + 4 * mzs * (-160 * (mws * mws * mws * mws) + 408 * (mws * mws * mws) * mzs - 302 * (mws * mws) * (mzs * mzs) + 70 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs)) * A0fin(mzs)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * ((-2 * (el * el) * (6 * (mts * mts) * (320 * (mws * mws) - 272 * mws * mzs + 39 * (mzs * mzs)) - 8 * (mzs * mzs) * (80 * (mws * mws) - 98 * mws * mzs + 45 * (mzs * mzs)) + mts * mzs * (1792 * (mws * mws) - 2368 * mws * mzs + 1269 * (mzs * mzs)))) / (-4 * mts + mzs) -
                           (9 * (el * el) * (mhs * mhs * (mzs * mzs) + 2 * mhs * (60 * (mws * mws * mws) + 96 * (mws * mws) * mzs - 98 * mws * (mzs * mzs) + mzs * mzs * mzs) - 4 * mzs * (120 * (mws * mws * mws) + 192 * (mws * mws) * mzs - 196 * mws * (mzs * mzs) + 3 * (mzs * mzs * mzs)))) / (mhs - 4 * mzs) -
                           1728 * delal * mws * (mws - mzs) * mzs * (Pi * Pi) + Cplx(0, 18) * (el * el) * (8 * (mws * mws * mws) - 64 * (mws * mws) * mzs + 106 * mws * (mzs * mzs) - 35 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (13824. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (3 * (el * el * el) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, 0, 0)) / (128. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs - 4 * mhs * mws + 12 * (mws * mws)) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mhs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (-mts + mws) * (mts + 2 * mws) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mts, 0)) / (256. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (2 * mws - mzs) * (4 * mws - mzs) * mzs * (12 * (mws * mws) + 20 * mws * mzs + mzs * mzs) * B0Refin(mws, mzs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (8 * (mws * mws * mws * mws) - 72 * (mws * mws * mws) * mzs + 210 * (mws * mws) * (mzs * mzs) - 160 * mws * (mzs * mzs * mzs) + 35 * (mzs * mzs * mzs * mzs)) *
                     B0Refin(mzs, 0, 0)) /
                        (768. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (mhs * mhs * mhs * (4 * mws - 3 * mzs) + 4 * mhs * (22 * mws - 15 * mzs) * (mzs * mzs) + 12 * (mzs * mzs * mzs) * (-11 * mws + 7 * mzs) + mhs * mhs * mzs * (-29 * mws + 21 * mzs)) * B0Refin(mzs, mhs, mzs)) / (1536. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * el * el * (mws * (mzs * mzs) * (-32 * (mws * mws) + 40 * mws * mzs - 17 * (mzs * mzs)) + mts * mzs * (64 * (mws * mws * mws) - 80 * (mws * mws) * mzs + 88 * mws * (mzs * mzs) - 27 * (mzs * mzs * mzs)) + 2 * (mts * mts) * (320 * (mws * mws * mws) - 592 * (mws * mws) * mzs + 275 * mws * (mzs * mzs) - 21 * (mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) /
                        (2304. * Power(mws, 1.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (168 * (mws * mws * mws) + 44 * (mws * mws) * mzs - 130 * mws * (mzs * mzs) + 17 * (mzs * mzs * mzs)) * B0Refin(mzs, mws, mws)) /
                        (1536. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (Cplx(0, 0.005208333333333333) * (el * el * el) * ((mws + mzs) * (mws + mzs)) * (2 * mws + mzs) *
                     C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * (mzs * mzs) * (16 * (mws * mws) - 20 * mws * mzs + 7 * (mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) /
                        (32. * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * ((mws + mzs) * (mws + mzs)) * (2 * mws + mzs) * C0Re(mzs, 0, 0, 0, 0, mws)) /
                        (192. * sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (mzs * mzs) * (16 * (mws * mws) - 20 * mws * mzs + 7 * (mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) /
                        (96. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        }
        break;
    case NEU:
        switch (formt)
        {
        case VEC:
            return ((el * el * el * mzs * (mhs * mhs * (19 * mws - 4 * mzs) * mzs - 24 * mhs * mws * (mzs * mzs) + 24 * mws * (mzs * mzs * mzs) + mhs * mhs * mhs * (-4 * mws + mzs)) * A0fin(mhs)) /
                        (1536. * mhs * Power(mws, 2.5) * (mhs - 4 * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (32 * (mws * mws) * ((mws - mzs) * (mws - mzs)) * (mzs * mzs) + 36 * (mts * mts * mts) * (mzs * mzs) * (-2 * mws + mzs) - 2 * mts * mws * mzs * (32 * (mws * mws * mws) - 88 * (mws * mws) * mzs + 65 * mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) + mts * mts * (-640 * (mws * mws * mws * mws) + 1184 * (mws * mws * mws) * mzs - 622 * (mws * mws) * (mzs * mzs) + 96 * mws * (mzs * mzs * mzs) - 9 * (mzs * mzs * mzs * mzs))) *
                     A0fin(mts)) /
                        (2304. * mts * Power(mws, 2.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (120 * (mws * mws * mws * mws) - 52 * (mws * mws * mws) * mzs - 62 * (mws * mws) * (mzs * mzs) - 2 * mws * (mhs - 5 * mzs) * (mzs * mzs) + mzs * mzs * mzs * (mhs + mzs)) * A0fin(mws)) / (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (mhs * mhs * mws * (4 * mws - 3 * mzs) - mhs * (48 * (mws * mws * mws) - 49 * (mws * mws) * mzs + 21 * mws * (mzs * mzs) + mzs * mzs * mzs) + 4 * mzs * (48 * (mws * mws * mws) - 62 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) + mzs * mzs * mzs)) * A0fin(mzs)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * ((-2 * (el * el) * (6 * (mts * mts) * (320 * (mws * mws) - 272 * mws * mzs + 39 * (mzs * mzs)) - 8 * (mzs * mzs) * (80 * (mws * mws) - 98 * mws * mzs + 45 * (mzs * mzs)) + mts * mzs * (1792 * (mws * mws) - 2368 * mws * mzs + 1269 * (mzs * mzs)))) / (-4 * mts + mzs) -
                           (9 * (el * el) * (mhs * mhs * (mzs * mzs) + 2 * mhs * (60 * (mws * mws * mws) + 84 * (mws * mws) * mzs - 84 * mws * (mzs * mzs) - mzs * mzs * mzs) + 4 * mzs * (-120 * (mws * mws * mws) - 168 * (mws * mws) * mzs + 168 * mws * (mzs * mzs) + mzs * mzs * mzs))) / (mhs - 4 * mzs) -
                           1728 * delal * mws * (mws - mzs) * mzs * (Pi * Pi) + Cplx(0, 54) * (el * el) * (8 * (mws * mws * mws) + 8 * (mws * mws) * mzs - 6 * mws * (mzs * mzs) - 5 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (13824. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (3 * (el * el * el) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, 0, 0)) / (128. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs - 4 * mhs * mws + 12 * (mws * mws)) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mhs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (-mts + mws) * (mts + 2 * mws) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mts, 0)) / (256. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (2 * mws - mzs) * (4 * mws - mzs) * mzs * (12 * (mws * mws) + 20 * mws * mzs + mzs * mzs) * B0Refin(mws, mzs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (24 * (mws * mws * mws * mws) - 2 * (mws * mws) * (mzs * mzs) - 16 * mws * (mzs * mzs * mzs) + 15 * (mzs * mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (768. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (mhs * mhs * mhs * (4 * mws - 3 * mzs) + 4 * mhs * (22 * mws - 15 * mzs) * (mzs * mzs) + 12 * (mzs * mzs * mzs) * (-11 * mws + 7 * mzs) + mhs * mhs * mzs * (-29 * mws + 21 * mzs)) * B0Refin(mzs, mhs, mzs)) / (1536. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * el * el * (mws * (mzs * mzs) * (-32 * (mws * mws) + 40 * mws * mzs - 17 * (mzs * mzs)) + mts * mzs * (64 * (mws * mws * mws) - 80 * (mws * mws) * mzs + 88 * mws * (mzs * mzs) - 27 * (mzs * mzs * mzs)) + 2 * (mts * mts) * (320 * (mws * mws * mws) - 592 * (mws * mws) * mzs + 275 * mws * (mzs * mzs) - 21 * (mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) /
                        (2304. * Power(mws, 1.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (168 * (mws * mws * mws) + 44 * (mws * mws) * mzs - 130 * mws * (mzs * mzs) + 17 * (mzs * mzs * mzs)) * B0Refin(mzs, mws, mws)) /
                        (1536. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (Cplx(0, 0.015625) * (el * el * el) * (2 * mws - mzs) * ((mws + mzs) * (mws + mzs)) *
                     C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.03125) * (el * el * el) * (mzs * mzs * mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) / (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (2 * mws - mzs) * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) / (64. * sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (mzs * mzs * mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (32. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        case AXV:
            return ((el * el * el * mzs * (mhs * mhs * (19 * mws - 4 * mzs) * mzs - 24 * mhs * mws * (mzs * mzs) + 24 * mws * (mzs * mzs * mzs) + mhs * mhs * mhs * (-4 * mws + mzs)) * A0fin(mhs)) /
                        (1536. * mhs * Power(mws, 2.5) * (mhs - 4 * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (32 * (mws * mws) * ((mws - mzs) * (mws - mzs)) * (mzs * mzs) + 36 * (mts * mts * mts) * (mzs * mzs) * (-2 * mws + mzs) - 2 * mts * mws * mzs * (32 * (mws * mws * mws) - 88 * (mws * mws) * mzs + 65 * mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) + mts * mts * (-640 * (mws * mws * mws * mws) + 1184 * (mws * mws * mws) * mzs - 622 * (mws * mws) * (mzs * mzs) + 96 * mws * (mzs * mzs * mzs) - 9 * (mzs * mzs * mzs * mzs))) *
                     A0fin(mts)) /
                        (2304. * mts * Power(mws, 2.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (120 * (mws * mws * mws * mws) - 52 * (mws * mws * mws) * mzs - 62 * (mws * mws) * (mzs * mzs) - 2 * mws * (mhs - 5 * mzs) * (mzs * mzs) + mzs * mzs * mzs * (mhs + mzs)) * A0fin(mws)) / (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (mhs * mhs * mws * (4 * mws - 3 * mzs) - mhs * (48 * (mws * mws * mws) - 49 * (mws * mws) * mzs + 21 * mws * (mzs * mzs) + mzs * mzs * mzs) + 4 * mzs * (48 * (mws * mws * mws) - 62 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) + mzs * mzs * mzs)) * A0fin(mzs)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * ((-2 * (el * el) * (6 * (mts * mts) * (320 * (mws * mws) - 272 * mws * mzs + 39 * (mzs * mzs)) - 8 * (mzs * mzs) * (80 * (mws * mws) - 98 * mws * mzs + 45 * (mzs * mzs)) + mts * mzs * (1792 * (mws * mws) - 2368 * mws * mzs + 1269 * (mzs * mzs)))) / (-4 * mts + mzs) -
                           (9 * (el * el) * (mhs * mhs * (mzs * mzs) + 2 * mhs * (60 * (mws * mws * mws) + 84 * (mws * mws) * mzs - 84 * mws * (mzs * mzs) - mzs * mzs * mzs) + 4 * mzs * (-120 * (mws * mws * mws) - 168 * (mws * mws) * mzs + 168 * mws * (mzs * mzs) + mzs * mzs * mzs))) / (mhs - 4 * mzs) -
                           1728 * delal * mws * (mws - mzs) * mzs * (Pi * Pi) + Cplx(0, 54) * (el * el) * (8 * (mws * mws * mws) + 8 * (mws * mws) * mzs - 6 * mws * (mzs * mzs) - 5 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (13824. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (3 * (el * el * el) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, 0, 0)) / (128. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (mhs * mhs - 4 * mhs * mws + 12 * (mws * mws)) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mhs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (-mts + mws) * (mts + 2 * mws) * (2 * mws - mzs) * (mzs * mzs) * B0Refin(mws, mts, 0)) / (256. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (2 * mws - mzs) * (4 * mws - mzs) * mzs * (12 * (mws * mws) + 20 * mws * mzs + mzs * mzs) * B0Refin(mws, mzs, mws)) /
                        (1536. * Power(mws, 2.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) -
                    (el * el * el * (24 * (mws * mws * mws * mws) - 2 * (mws * mws) * (mzs * mzs) - 16 * mws * (mzs * mzs * mzs) + 15 * (mzs * mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (768. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (mhs * mhs * mhs * (4 * mws - 3 * mzs) + 4 * mhs * (22 * mws - 15 * mzs) * (mzs * mzs) + 12 * (mzs * mzs * mzs) * (-11 * mws + 7 * mzs) + mhs * mhs * mzs * (-29 * mws + 21 * mzs)) * B0Refin(mzs, mhs, mzs)) / (1536. * Power(mws, 1.5) * Power(-mws + mzs, 2.5) * (-mhs + 4 * mzs) * (Pi * Pi)) +
                    (el * el * el * (mws * (mzs * mzs) * (-32 * (mws * mws) + 40 * mws * mzs - 17 * (mzs * mzs)) + mts * mzs * (64 * (mws * mws * mws) - 80 * (mws * mws) * mzs + 88 * mws * (mzs * mzs) - 27 * (mzs * mzs * mzs)) + 2 * (mts * mts) * (320 * (mws * mws * mws) - 592 * (mws * mws) * mzs + 275 * mws * (mzs * mzs) - 21 * (mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) /
                        (2304. * Power(mws, 1.5) * (-4 * mts + mzs) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (el * el * el * (168 * (mws * mws * mws) + 44 * (mws * mws) * mzs - 130 * mws * (mzs * mzs) + 17 * (mzs * mzs * mzs)) * B0Refin(mzs, mws, mws)) /
                        (1536. * sqrt(mws) * Power(-mws + mzs, 2.5) * (Pi * Pi)) +
                    (Cplx(0, 0.015625) * (el * el * el) * (2 * mws - mzs) * ((mws + mzs) * (mws + mzs)) *
                     C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.03125) * (el * el * el) * (mzs * mzs * mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) / (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (2 * mws - mzs) * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) / (64. * sqrt(mws) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (mzs * mzs * mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (32. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        }
        break;
    }
}

//Cplx Zp_SMNLO::result(void) const
Cplx Z1p(const int ftyp, const int formt, const inval &input)
{
    const inval *ival;
    ival = &input;
    double //el = sqrt(4 * Pi * ival->get(al)),
        el = 0.303,
        mz = ival->get(MZ),
        mw = ival->get(MW),
        mt = ival->get(MT),
        mh = ival->get(MH),
        delal = ival->get(Delal),
        cw = mw / mz,
        sw = sqrt(1 - cw * cw),
        mws = mw * mw,
        mhs = mh * mh,
        mts = mt * mt,
        mzs = mz * mz;
    switch (ftyp)
    {
    case LEP:
        switch (formt)
        {
        case VEC:
            return (-(el * el * el * mts * (8 * mws - 5 * mzs) * A0fin(mts)) / (6. * sqrt(mws) * (mzs * mzs) * (-4 * mts + mzs) * sqrt(-mws + mzs) * (Pi * Pi)) +
                    (el * el * el * (96 * (mws * mws * mws * mws) - 40 * (mws * mws * mws) * mzs - 44 * (mws * mws) * (mzs * mzs) + 2 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs) * A0fin(mws)) /
                        (64. * Power(mws, 1.5) * (4 * mws - mzs) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (3 * (el * el * el) * (16 * (mws * mws * mws) - 36 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) * A0fin(mzs)) /
                        (128. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * ((64 * mws * (8 * mws - 5 * mzs) * (mws - mzs) * (-2 * (mts * mts) - 3 * mts * mzs + mzs * mzs)) / (-4 * mts + mzs) + (-1152 * (mws * mws * mws * mws * mws) + 672 * (mws * mws * mws * mws) * mzs + 1120 * (mws * mws * mws) * (mzs * mzs) - 680 * (mws * mws) * (mzs * mzs * mzs) + 184 * mws * (mzs * mzs * mzs * mzs) - 27 * (mzs * mzs * mzs * mzs * mzs)) / (4 * mws - mzs) + Cplx(0, 6) * mzs * (-48 * (mws * mws * mws) + 112 * (mws * mws) * mzs - 88 * mws * (mzs * mzs) + 27 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (768. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (48 * (mws * mws * mws) - 112 * (mws * mws) * mzs + 88 * mws * (mzs * mzs) - 27 * (mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (128. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (mts * mts) * (8 * mws - 5 * mzs) * B0Refin(mzs, mts, mts)) / (6. * sqrt(mws) * (mzs * mzs) * (-4 * mts + mzs) * sqrt(-mws + mzs) * (Pi * Pi)) -
                    (el * el * el * sqrt(mws) * (24 * (mws * mws * mws) - 10 * (mws * mws) * mzs - 9 * mws * (mzs * mzs) + mzs * mzs * mzs) * B0Refin(mzs, mws, mws)) /
                        (16. * (4 * mws - mzs) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (Cplx(0, 0.03125) * (el * el * el) * sqrt(mws) * (mws + mzs) * C0Im(mzs, 0, 0, 0, 0, mws)) / (Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.03125) * (el * el * el) * (16 * (mws * mws * mws) - 36 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (mws + mzs) * C0Re(0, mzs, 0, 0, mws, mws)) /
                        (16. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * sqrt(mws) * (mws + mzs) * C0Re(mzs, 0, 0, 0, 0, mws)) / (32. * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (16 * (mws * mws * mws) - 36 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) /
                        (32. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        case AXV:
            return ((el * el * el * (16 * (mws * mws * mws) - 2 * mws * (mzs * mzs) + mzs * mzs * mzs) * A0fin(mws)) /
                        (64. * Power(mws, 1.5) * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (3 * (el * el * el) * (12 * (mws * mws) - 18 * mws * mzs + 7 * (mzs * mzs)) * A0fin(mzs)) / (128. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (16 * (mws * mws * mws) - 96 * (mws * mws) * mzs + 48 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs) + Cplx(0, 2) * (128 * (mws * mws * mws) - 256 * (mws * mws) * mzs + 140 * mws * (mzs * mzs) - 21 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (256. * Power(mws, 1.5) * (4 * mws - mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (32 * (mws * mws) - 56 * mws * mzs + 21 * (mzs * mzs)) * B0Refin(mzs, 0, 0)) / (128. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * Power(mws, 1.5) * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (8. * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (Cplx(0, 0.03125) * (el * el * el) * sqrt(mws) * (mws + mzs) * C0Im(mzs, 0, 0, 0, 0, mws)) / (Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.03125) * (el * el * el) * mzs * (12 * (mws * mws) - 18 * mws * mzs + 7 * (mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (mws + mzs) * C0Re(0, mzs, 0, 0, mws, mws)) /
                        (16. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * sqrt(mws) * (mws + mzs) * C0Re(mzs, 0, 0, 0, 0, mws)) / (32. * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * mzs * (12 * (mws * mws) - 18 * mws * mzs + 7 * (mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (32. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        }
        break;
    case UQU:
    case CQU:
        switch (formt)
        {
        case VEC:
            return ((el * el * el * mts * (8 * mws - 5 * mzs) * A0fin(mts)) / (9. * sqrt(mws) * (mzs * mzs) * (-4 * mts + mzs) * sqrt(-mws + mzs) * (Pi * Pi)) -
                    (el * el * el * (192 * (mws * mws * mws * mws) - 80 * (mws * mws * mws) * mzs - 76 * (mws * mws) * (mzs * mzs) + 8 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs) * A0fin(mws)) /
                        (192. * Power(mws, 1.5) * (4 * mws - mzs) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (128 * (mws * mws * mws) - 240 * (mws * mws) * mzs + 204 * mws * (mzs * mzs) - 65 * (mzs * mzs * mzs)) * A0fin(mzs)) /
                        (1152. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * ((384 * mws * (8 * mws - 5 * mzs) * (-mws + mzs) * (-2 * (mts * mts) - 3 * mts * mzs + mzs * mzs)) / (-4 * mts + mzs) + (6912 * (mws * mws * mws * mws * mws) - 4672 * (mws * mws * mws * mws) * mzs - 5360 * (mws * mws * mws) * (mzs * mzs) + 2400 * (mws * mws) * (mzs * mzs * mzs) - 398 * mws * (mzs * mzs * mzs * mzs) + 65 * (mzs * mzs * mzs * mzs * mzs)) / (4 * mws - mzs) + Cplx(0, 6) * mzs * (104 * (mws * mws * mws) - 264 * (mws * mws) * mzs + 198 * mws * (mzs * mzs) - 65 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (6912. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (104 * (mws * mws * mws) - 264 * (mws * mws) * mzs + 198 * mws * (mzs * mzs) - 65 * (mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (1152. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (mts * mts) * (8 * mws - 5 * mzs) * B0Refin(mzs, mts, mts)) / (9. * sqrt(mws) * (mzs * mzs) * (-4 * mts + mzs) * sqrt(-mws + mzs) * (Pi * Pi)) +
                    (el * el * el * sqrt(mws) * (24 * (mws * mws * mws) - 8 * (mws * mws) * mzs - 8 * mws * (mzs * mzs) + mzs * mzs * mzs) * B0Refin(mzs, mws, mws)) /
                        (24. * (4 * mws - mzs) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * sqrt(mws) * (mws + mzs) * (2 * mws + mzs) * C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (Cplx(0, 0.0011574074074074073) * (el * el * el) *
                     (128 * (mws * mws * mws) - 240 * (mws * mws) * mzs + 204 * mws * (mzs * mzs) - 65 * (mzs * mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * Power(mws, 1.5) * (mws + mzs) * C0Re(0, mzs, 0, 0, mws, mws)) /
                        (16. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * sqrt(mws) * (mws + mzs) * (2 * mws + mzs) * C0Re(mzs, 0, 0, 0, 0, mws)) /
                        (96. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (128 * (mws * mws * mws) - 240 * (mws * mws) * mzs + 204 * mws * (mzs * mzs) - 65 * (mzs * mzs * mzs)) *
                     C0Re(mzs, 0, 0, 0, 0, mzs)) /
                        (864. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        case AXV:
            return (-(el * el * el * (32 * (mws * mws * mws) + 12 * (mws * mws) * mzs + mzs * mzs * mzs) * A0fin(mws)) /
                        (192. * Power(mws, 1.5) * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (16 * (mws * mws) - 20 * mws * mzs + 7 * (mzs * mzs)) * A0fin(mzs)) / (128. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (mzs * (48 * (mws * mws * mws) + 112 * (mws * mws) * mzs - 50 * mws * (mzs * mzs) + 7 * (mzs * mzs * mzs)) + Cplx(0, 2) * (32 * (mws * mws * mws * mws) - 168 * (mws * mws * mws) * mzs + 288 * (mws * mws) * (mzs * mzs) - 146 * mws * (mzs * mzs * mzs) + 21 * (mzs * mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) / (768. * Power(mws, 1.5) * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (8 * (mws * mws * mws) - 40 * (mws * mws) * mzs + 62 * mws * (mzs * mzs) - 21 * (mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (384. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (8. * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * sqrt(mws) * (mws + mzs) * (2 * mws + mzs) * C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * mzs * (16 * (mws * mws) - 20 * mws * mzs + 7 * (mzs * mzs)) *
                     C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * Power(mws, 1.5) * (mws + mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (16. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * sqrt(mws) * (mws + mzs) * (2 * mws + mzs) * C0Re(mzs, 0, 0, 0, 0, mws)) / (96. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * mzs * (16 * (mws * mws) - 20 * mws * mzs + 7 * (mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (96. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
        }
        break;
    case DQU:
    case SQU:
        switch (formt)
        {
        case VEC:
            return (-(el * el * el * mts * (8 * mws - 5 * mzs) * A0fin(mts)) / (18. * sqrt(mws) * (mzs * mzs) * (-4 * mts + mzs) * sqrt(-mws + mzs) * (Pi * Pi)) +
                    (el * el * el * (24 * (mws * mws * mws) - 4 * (mws * mws) * mzs - 6 * mws * (mzs * mzs) + mzs * mzs * mzs) * A0fin(mws)) /
                        (192. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (16 * (mws * mws * mws) - 12 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) * A0fin(mzs)) /
                        (1152. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (-864 * (mws * mws * mws * mws) + 416 * (mws * mws * mws) * mzs + 816 * (mws * mws) * (mzs * mzs) - 24 * mws * (mzs * mzs * mzs) + 7 * (mzs * mzs * mzs * mzs) + (192 * mws * (8 * mws - 5 * mzs) * (mws - mzs) * (-2 * (mts * mts) - 3 * mts * mzs + mzs * mzs)) / (-4 * mts + mzs) + Cplx(0, 6) * (4 * mws - mzs) * mzs * (8 * (mws * mws) + 8 * mws * mzs - 7 * (mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (6912. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (4 * mws - mzs) * (8 * (mws * mws) + 8 * mws * mzs - 7 * (mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (1152. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (mts * mts) * (8 * mws - 5 * mzs) * B0Refin(mzs, mts, mts)) / (18. * sqrt(mws) * (mzs * mzs) * (-4 * mts + mzs) * sqrt(-mws + mzs) * (Pi * Pi)) -
                    (el * el * el * sqrt(mws) * (6 * (mws * mws) + mws * mzs - mzs * mzs) * B0Refin(mzs, mws, mws)) / (48. * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * sqrt(mws) * (4 * mws - mzs) * (mws + mzs) * C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.0011574074074074073) * (el * el * el) *
                     (16 * (mws * mws * mws) - 12 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (mws + mzs) * C0Re(0, mzs, 0, 0, mws, mws)) /
                        (16. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * sqrt(mws) * (4 * mws - mzs) * (mws + mzs) * C0Re(mzs, 0, 0, 0, 0, mws)) /
                        (96. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (16 * (mws * mws * mws) - 12 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) *
                     C0Re(mzs, 0, 0, 0, 0, mzs)) /
                        (864. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        case AXV:
            return ((el * el * el * (16 * (mws * mws * mws) + 24 * (mws * mws) * mzs + 6 * mws * (mzs * mzs) - mzs * mzs * mzs) * A0fin(mws)) /
                        (192. * Power(mws, 1.5) * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * A0fin(mzs)) / (128. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (-(mzs * (112 * (mws * mws * mws) + 8 * (mws * mws) * mzs - 4 * mws * (mzs * mzs) + mzs * mzs * mzs)) - Cplx(0, 2) * (64 * (mws * mws * mws * mws) - 48 * (mws * mws * mws) * mzs + 24 * (mws * mws) * (mzs * mzs) - 16 * mws * (mzs * mzs * mzs) + 3 * (mzs * mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) / (768. * Power(mws, 1.5) * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (16 * (mws * mws * mws) - 8 * (mws * mws) * mzs + 4 * mws * (mzs * mzs) - 3 * (mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (384. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * Power(mws, 1.5) * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (8. * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * sqrt(mws) * (4 * mws - mzs) * (mws + mzs) * C0Im(mzs, 0, 0, 0, 0, mws)) /
                        (mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * mzs * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) *
                     C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (mws + mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (16. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * sqrt(mws) * (4 * mws - mzs) * (mws + mzs) * C0Re(mzs, 0, 0, 0, 0, mws)) / (96. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * mzs * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (96. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));

            break;
        }
        break;
    case BQU:
        switch (formt)
        {
        case VEC:
            return ((el * el * el * (32 * mts * mws * (8 * mws - 5 * mzs) * (mws - mzs) + (3 * mzs * (4 * (mts * mts * mts * mts * mts * mts) * (2 * mws + mzs) - 2 * (mts * mts * mts * mts * mts) * (8 * (mws * mws) + 5 * mws * mzs - 4 * (mzs * mzs)) - 2 * (mws * mws * mws * mws) * mzs * (2 * (mws * mws) + 9 * mws * mzs + 7 * (mzs * mzs)) + mts * mts * mts * mws * (64 * (mws * mws * mws) + 84 * (mws * mws) * mzs + 5 * mws * (mzs * mzs) - 6 * (mzs * mzs * mzs)) + mts * mts * mts * mts * (-16 * (mws * mws * mws) - 4 * (mws * mws) * mzs + mws * (mzs * mzs) + 4 * (mzs * mzs * mzs)) + mts * (mws * mws) * (16 * (mws * mws * mws * mws) + 70 * (mws * mws * mws) * mzs + 99 * (mws * mws) * (mzs * mzs) + 9 * mws * (mzs * mzs * mzs) - 8 * (mzs * mzs * mzs * mzs)) + mts * mts * mws * (-56 * (mws * mws * mws * mws) - 140 * (mws * mws * mws) * mzs - 71 * (mws * mws) * (mzs * mzs) + 31 * mws * (mzs * mzs * mzs) + 5 * (mzs * mzs * mzs * mzs)))) / ((mts * mts + mws * mws + mts * (-2 * mws + mzs)) * (mts * mts - 2 * mts * mws + mws * (mws + mzs)))) * A0fin(mts)) /
                        (576. * Power(mws, 1.5) * (mzs * mzs) * (-4 * mts + mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (2 * (mts * mts * mts * mts * mts) * mzs * (2 * mws + mzs) + 3 * (mts * mts * mts * mts) * (16 * (mws * mws * mws) - 8 * (mws * mws) * mzs - 12 * mws * (mzs * mzs) + mzs * mzs * mzs) + mts * mts * mts * (-192 * (mws * mws * mws * mws) + 104 * (mws * mws * mws) * mzs + 108 * (mws * mws) * (mzs * mzs) - 12 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs) + 2 * (mws * mws * mws) * (24 * (mws * mws * mws * mws) + 20 * (mws * mws * mws) * mzs - 10 * (mws * mws) * (mzs * mzs) - 5 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs) + mts * mts * mws * (288 * (mws * mws * mws * mws) - 112 * (mws * mws * mws) * mzs - 136 * (mws * mws) * (mzs * mzs) - 11 * mws * (mzs * mzs * mzs) + 4 * (mzs * mzs * mzs * mzs)) - 2 * mts * mws * (96 * (mws * mws * mws * mws * mws) + 6 * (mws * mws * mws * mws) * mzs - 65 * (mws * mws * mws) * (mzs * mzs) - 13 * (mws * mws) * (mzs * mzs * mzs) + 7 * mws * (mzs * mzs * mzs * mzs) - mzs * mzs * mzs * mzs * mzs)) * A0fin(mws)) /
                        (384. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (mts * mts + mws * mws + mts * (-2 * mws + mzs)) * (mts * mts - 2 * mts * mws + mws * (mws + mzs)) *
                         (Pi * Pi)) +
                    (el * el * el * (16 * (mws * mws * mws) - 12 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) * A0fin(mzs)) /
                        (1152. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * ((192 * mws * (8 * mws - 5 * mzs) * (mws - mzs) * (-2 * (mts * mts) - 3 * mts * mzs + mzs * mzs)) / (-4 * mts + mzs) - (72 * (mts * mts * mts * mts * mts * mts) * (5 * mws - 2 * mzs) * (mzs * mzs) + mts * mts * mts * mts * mts * (-3456 * (mws * mws * mws * mws) + 1664 * (mws * mws * mws) * mzs + 1824 * (mws * mws) * (mzs * mzs) - 150 * mws * (mzs * mzs * mzs) - 125 * (mzs * mzs * mzs * mzs)) + mws * mws * mws * mzs * (864 * (mws * mws * mws * mws * mws) + 448 * (mws * mws * mws * mws) * mzs - 1232 * (mws * mws * mws) * (mzs * mzs) - 792 * (mws * mws) * (mzs * mzs * mzs) + 17 * mws * (mzs * mzs * mzs * mzs) - 7 * (mzs * mzs * mzs * mzs * mzs)) + 4 * (mts * mts * mts * mts) * (3456 * (mws * mws * mws * mws * mws) - 2312 * (mws * mws * mws * mws) * mzs - 2412 * (mws * mws * mws) * (mzs * mzs) + 951 * (mws * mws) * (mzs * mzs * mzs) + 17 * mws * (mzs * mzs * mzs * mzs) + 3 * (mzs * mzs * mzs * mzs * mzs)) - mts * mts * mts * (20736 * (mws * mws * mws * mws * mws * mws) - 9984 * (mws * mws * mws * mws * mws) * mzs - 19008 * (mws * mws * mws * mws) * (mzs * mzs) + 2144 * (mws * mws * mws) * (mzs * mzs * mzs) + 1332 * (mws * mws) * (mzs * mzs * mzs * mzs) + 129 * mws * (mzs * mzs * mzs * mzs * mzs) + 7 * (mzs * mzs * mzs * mzs * mzs * mzs)) + mts * mts * mws * (13824 * (mws * mws * mws * mws * mws * mws) + 1984 * (mws * mws * mws * mws * mws) * mzs - 21176 * (mws * mws * mws * mws) * (mzs * mzs) - 5372 * (mws * mws * mws) * (mzs * mzs * mzs) + 4532 * (mws * mws) * (mzs * mzs * mzs * mzs) - 55 * mws * (mzs * mzs * mzs * mzs * mzs) + 26 * (mzs * mzs * mzs * mzs * mzs * mzs)) + mts * mws * (-3456 * (mws * mws * mws * mws * mws * mws * mws) - 5248 * (mws * mws * mws * mws * mws * mws) * mzs + 5728 * (mws * mws * mws * mws * mws) * (mzs * mzs) + 7766 * (mws * mws * mws * mws) * (mzs * mzs * mzs) + 137 * (mws * mws * mws) * (mzs * mzs * mzs * mzs) - 793 * (mws * mws) * (mzs * mzs * mzs * mzs * mzs) + 31 * mws * (mzs * mzs * mzs * mzs * mzs * mzs) - 7 * (mzs * mzs * mzs * mzs * mzs * mzs * mzs))) / ((-4 * mts + mzs) * (mts * mts + mws * mws + mts * (-2 * mws + mzs)) * (mts * mts - 2 * mts * mws + mws * (mws + mzs))) + Cplx(0, 6) * mzs * (-16 * (mws * mws * mws) + 12 * (mws * mws) * mzs - 30 * mws * (mzs * mzs) + 7 * (mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (6912. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (16 * (mws * mws * mws) - 12 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (1152. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (96 * (mts * mts * mts * mts * mts) * (mws - mzs) * mzs + 2 * (mts * mts * mts) * mws * (512 * (mws * mws * mws) - 976 * (mws * mws) * mzs + 212 * mws * (mzs * mzs) - 45 * (mzs * mzs * mzs)) + 6 * (mws * mws) * (mzs * mzs) * (8 * (mws * mws * mws) + 10 * (mws * mws) * mzs + mws * (mzs * mzs) - mzs * mzs * mzs) - 12 * mts * (mws * mws) * mzs * (16 * (mws * mws * mws) + 30 * (mws * mws) * mzs + 6 * mws * (mzs * mzs) - mzs * mzs * mzs) - 2 * (mts * mts * mts * mts) * (256 * (mws * mws * mws) - 368 * (mws * mws) * mzs + 28 * mws * (mzs * mzs) + 3 * (mzs * mzs * mzs)) + mts * mts * mws * (-512 * (mws * mws * mws * mws) + 800 * (mws * mws * mws) * mzs + 872 * (mws * mws) * (mzs * mzs) - 236 * mws * (mzs * mzs * mzs) - 15 * (mzs * mzs * mzs * mzs))) * B0Refin(mzs, mts, mts)) /
                        (1152. * Power(mws, 1.5) * (mzs * mzs) * (-4 * mts + mzs) * Power(-mws + mzs, 1.5) * (mts * mts - 2 * mts * mws + mws * (mws + mzs)) * (Pi * Pi)) -
                    (el * el * el * (3 * (mts * mts * mts * mts) * mzs * (-2 * mws + mzs) + 4 * (mws * mws * mws * mws) * (6 * (mws * mws) + mws * mzs - mzs * mzs) + 3 * (mts * mts * mts) * mzs * (2 * (mws * mws) - 5 * mws * mzs + mzs * mzs) + mts * mts * mws * (24 * (mws * mws * mws) + 10 * (mws * mws) * mzs - mws * (mzs * mzs) - 9 * (mzs * mzs * mzs)) + mts * (mws * mws) * (-48 * (mws * mws * mws) + 10 * (mws * mws) * mzs + 9 * mws * (mzs * mzs) - 4 * (mzs * mzs * mzs))) * B0Refin(mzs, mws, mws)) /
                        (192. * Power(mws, 1.5) * (mzs * mzs) * Power(-mws + mzs, 1.5) * (mts * mts + mws * mws + mts * (-2 * mws + mzs)) * (Pi * Pi)) -
                    (Cplx(0, 0.0011574074074074073) * (el * el * el) * (16 * (mws * mws * mws) - 12 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) *
                     C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (mts * mts * mts * (4 * mws - 2 * mzs) + mts * mts * (6 * mws - mzs) * mzs + 8 * (mws * mws * mws) * (mws + mzs) - 2 * mts * mws * (6 * (mws * mws) + 3 * mws * mzs - 2 * (mzs * mzs))) * C0Re(0, mzs, 0, mts, mws, mws)) /
                        (128. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (16 * (mws * mws * mws) - 12 * (mws * mws) * mzs + 30 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) /
                        (864. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (8 * (mts * mts * mts) * (mws - mzs) + mts * mts * (16 * mws - mzs) * mzs + 4 * (mws * mws) * (4 * (mws * mws) + 3 * mws * mzs - mzs * mzs) - 4 * mts * mws * (6 * (mws * mws) + 2 * mws * mzs + mzs * mzs)) * C0Re(mzs, 0, 0, mts, mts, mws)) / (384. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;

        case AXV:
            return (
                (el * el * el * (4 * (mts * mts * mts * mts * mts * mts) * (2 * mws + mzs) - 2 * (mts * mts * mts * mts * mts) * (8 * (mws * mws) + 5 * mws * mzs - 4 * (mzs * mzs)) - 2 * (mws * mws * mws * mws) * mzs * (2 * (mws * mws) + 9 * mws * mzs + 7 * (mzs * mzs)) + mts * mts * mts * mws * (64 * (mws * mws * mws) + 84 * (mws * mws) * mzs + 5 * mws * (mzs * mzs) - 6 * (mzs * mzs * mzs)) + mts * mts * mts * mts * (-16 * (mws * mws * mws) - 4 * (mws * mws) * mzs + mws * (mzs * mzs) + 4 * (mzs * mzs * mzs)) + mts * (mws * mws) * (16 * (mws * mws * mws * mws) + 70 * (mws * mws * mws) * mzs + 99 * (mws * mws) * (mzs * mzs) + 9 * mws * (mzs * mzs * mzs) - 8 * (mzs * mzs * mzs * mzs)) + mts * mts * mws * (-56 * (mws * mws * mws * mws) - 140 * (mws * mws * mws) * mzs - 71 * (mws * mws) * (mzs * mzs) + 31 * mws * (mzs * mzs * mzs) + 5 * (mzs * mzs * mzs * mzs))) * A0fin(mts)) /
                    (192. * Power(mws, 1.5) * mzs * (-4 * mts + mzs) * Power(-mws + mzs, 1.5) * (mts * mts + mws * mws + mts * (-2 * mws + mzs)) *
                     (mts * mts - 2 * mts * mws + mws * (mws + mzs)) * (Pi * Pi)) -
                (el * el * el * (-2 * (mts * mts * mts * mts * mts) * (8 * (mws * mws) + 2 * mws * mzs - mzs * mzs) + mts * mts * mts * mts * (32 * (mws * mws * mws) + 32 * (mws * mws) * mzs - 40 * mws * (mzs * mzs) + 3 * (mzs * mzs * mzs)) - 2 * (mws * mws * mws) * (16 * (mws * mws * mws * mws) + 40 * (mws * mws * mws) * mzs + 30 * (mws * mws) * (mzs * mzs) + 5 * mws * (mzs * mzs * mzs) - mzs * mzs * mzs * mzs) + mts * mts * mts * (32 * (mws * mws * mws * mws) - 88 * (mws * mws * mws) * mzs + 36 * (mws * mws) * (mzs * mzs) - 8 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs) + mts * mts * mws * (-128 * (mws * mws * mws * mws) + 16 * (mws * mws * mws) * mzs + 44 * (mws * mws) * (mzs * mzs) - 35 * mws * (mzs * mzs * mzs) + 4 * (mzs * mzs * mzs * mzs)) + 2 * mts * mws * (56 * (mws * mws * mws * mws * mws) + 62 * (mws * mws * mws * mws) * mzs - 15 * (mws * mws * mws) * (mzs * mzs) - 7 * (mws * mws) * (mzs * mzs * mzs) - 7 * mws * (mzs * mzs * mzs * mzs) + mzs * mzs * mzs * mzs * mzs)) * A0fin(mws)) /
                    (384. * Power(mws, 1.5) * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (mts * mts + mws * mws + mts * (-2 * mws + mzs)) *
                     (mts * mts - 2 * mts * mws + mws * (mws + mzs)) * (Pi * Pi)) +
                (el * el * el * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * A0fin(mzs)) / (128. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                (el * el * el * ((-8 * (mts * mts * mts * mts * mts * mts) * (20 * (mws * mws) - 13 * mws * mzs + 2 * (mzs * mzs)) + 3 * (mts * mts * mts * mts * mts) * (64 * (mws * mws * mws) - 56 * (mws * mws) * mzs + 26 * mws * (mzs * mzs) - 7 * (mzs * mzs * mzs)) + mts * mts * mts * mts * (832 * (mws * mws * mws * mws) - 400 * (mws * mws * mws) * mzs - 92 * (mws * mws) * (mzs * mzs) + 60 * mws * (mzs * mzs * mzs) - 4 * (mzs * mzs * mzs * mzs)) + mws * mws * mws * mzs * (112 * (mws * mws * mws * mws) + 120 * (mws * mws * mws) * mzs + 4 * (mws * mws) * (mzs * mzs) - 3 * mws * (mzs * mzs * mzs) + mzs * mzs * mzs * mzs) + mts * mts * mws * (1632 * (mws * mws * mws * mws * mws) + 1144 * (mws * mws * mws * mws) * mzs - 652 * (mws * mws * mws) * (mzs * mzs) - 88 * (mws * mws) * (mzs * mzs * mzs) + 49 * mws * (mzs * mzs * mzs * mzs) - 6 * (mzs * mzs * mzs * mzs * mzs)) + mts * mts * mts * (-2048 * (mws * mws * mws * mws * mws) + 160 * (mws * mws * mws * mws) * mzs + 384 * (mws * mws * mws) * (mzs * mzs) - 24 * (mws * mws) * (mzs * mzs * mzs) - 21 * mws * (mzs * mzs * mzs * mzs) + mzs * mzs * mzs * mzs * mzs) + mts * mws * (-448 * (mws * mws * mws * mws * mws * mws) - 952 * (mws * mws * mws * mws * mws) * mzs - 110 * (mws * mws * mws * mws) * (mzs * mzs) + 125 * (mws * mws * mws) * (mzs * mzs * mzs) + 3 * (mws * mws) * (mzs * mzs * mzs * mzs) - 5 * mws * (mzs * mzs * mzs * mzs * mzs) + mzs * mzs * mzs * mzs * mzs * mzs)) / ((4 * mws - mzs) * (-4 * mts + mzs) * (mts * mts + mws * mws + mts * (-2 * mws + mzs)) * (mts * mts - 2 * mts * mws + mws * (mws + mzs))) - Cplx(0, 6) * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * B0Im(mzs, 0, 0))) / (768. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                (el * el * el * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * B0Refin(mzs, 0, 0)) / (128. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                (el * el * el * (32 * (mts * mts * mts * mts * mts) * (mws - mzs) - 2 * (mts * mts * mts * mts) * (16 * (mws * mws) - 44 * mws * mzs + mzs * mzs) - 6 * (mts * mts * mts) * mws * (16 * (mws * mws) + 12 * mws * mzs + 5 * (mzs * mzs)) + mts * mts * mws * (160 * (mws * mws * mws) + 120 * (mws * mws) * mzs + 28 * mws * (mzs * mzs) - 5 * (mzs * mzs * mzs)) + 2 * (mws * mws) * mzs * (8 * (mws * mws * mws) + 10 * (mws * mws) * mzs + mws * (mzs * mzs) - mzs * mzs * mzs) - 4 * mts * (mws * mws) * (16 * (mws * mws * mws) + 30 * (mws * mws) * mzs + 6 * mws * (mzs * mzs) - mzs * mzs * mzs)) * B0Refin(mzs, mts, mts)) /
                    (384. * Power(mws, 1.5) * mzs * (-4 * mts + mzs) * Power(-mws + mzs, 1.5) * (mts * mts - 2 * mts * mws + mws * (mws + mzs)) * (Pi * Pi)) -
                (el * el * el * (8 * (mws * mws * mws * mws * mws) * (2 * mws + mzs) - mts * mts * mts * mts * (8 * (mws * mws) - 6 * mws * mzs + mzs * mzs) + mts * (mws * mws * mws) * (-40 * (mws * mws) - 2 * mws * mzs + 9 * (mzs * mzs)) + mts * mts * mts * (8 * (mws * mws * mws) - 22 * (mws * mws) * mzs + 9 * mws * (mzs * mzs) - mzs * mzs * mzs) + mts * mts * mws * (24 * (mws * mws * mws) + 10 * (mws * mws) * mzs - 13 * mws * (mzs * mzs) + 3 * (mzs * mzs * mzs))) * B0Refin(mzs, mws, mws)) /
                    (64. * Power(mws, 1.5) * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (mts * mts + mws * mws + mts * (-2 * mws + mzs)) * (Pi * Pi)) -
                (Cplx(0, 0.010416666666666666) * (el * el * el) * mzs * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                    (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                (el * el * el *
                 (mts * mts * mts * (4 * mws - 2 * mzs) + mts * mts * (6 * mws - mzs) * mzs + 8 * (mws * mws * mws) * (mws + mzs) -
                  2 * mts * mws * (6 * (mws * mws) + 3 * mws * mzs - 2 * (mzs * mzs))) *
                 C0Re(0, mzs, 0, mts, mws, mws)) /
                    (128. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                (el * el * el * mzs * (4 * (mws * mws) - 2 * mws * mzs + mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (96. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                (el * el * el * (8 * (mts * mts * mts) * (mws - mzs) + mts * mts * (16 * mws - mzs) * mzs + 4 * (mws * mws) * (4 * (mws * mws) + 3 * mws * mzs - mzs * mzs) - 4 * mts * mws * (6 * (mws * mws) + 2 * mws * mzs + mzs * mzs)) * C0Re(mzs, 0, 0, mts, mts, mws)) / (384. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        }
        break;
    case NEU:
        switch (formt)
        {
        case VEC:
            return (-(el * el * el * (12 * (mws * mws) + 4 * mws * mzs - mzs * mzs) * A0fin(mws)) / (64. * Power(mws, 1.5) * (4 * mws - mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (3 * (el * el * el) * mzs * A0fin(mzs)) / (128. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (mzs * (48 * (mws * mws * mws) - 8 * (mws * mws) * mzs - 2 * mws * (mzs * mzs) + mzs * mzs * mzs) + Cplx(0, 2) * (32 * (mws * mws * mws * mws) - 8 * (mws * mws * mws) * mzs - 8 * (mws * mws) * (mzs * mzs) - 10 * mws * (mzs * mzs * mzs) + 3 * (mzs * mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) / (256. * Power(mws, 1.5) * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (8 * (mws * mws * mws) - 2 * mws * (mzs * mzs) - 3 * (mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (128. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (8. * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.03125) * (el * el * el) * sqrt(mws) * (2 * mws - mzs) * (mws + mzs) * C0Im(mzs, 0, 0, 0, 0, mws)) / (mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (Cplx(0, 0.03125) * (el * el * el) * (mzs * mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) / (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * Power(mws, 1.5) * (mws + mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (16. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * sqrt(mws) * (2 * mws - mzs) * (mws + mzs) * C0Re(mzs, 0, 0, 0, 0, mws)) / (32. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (mzs * mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (32. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));
            break;
        case AXV:
            return (-(el * el * el * (12 * (mws * mws) + 4 * mws * mzs - mzs * mzs) * A0fin(mws)) / (64. * Power(mws, 1.5) * (4 * mws - mzs) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (3 * (el * el * el) * mzs * A0fin(mzs)) / (128. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (mzs * (48 * (mws * mws * mws) - 8 * (mws * mws) * mzs - 2 * mws * (mzs * mzs) + mzs * mzs * mzs) + Cplx(0, 2) * (32 * (mws * mws * mws * mws) - 8 * (mws * mws * mws) * mzs - 8 * (mws * mws) * (mzs * mzs) - 10 * mws * (mzs * mzs * mzs) + 3 * (mzs * mzs * mzs * mzs)) * B0Im(mzs, 0, 0))) / (256. * Power(mws, 1.5) * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * (8 * (mws * mws * mws) - 2 * mws * (mzs * mzs) - 3 * (mzs * mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (128. * Power(mws, 1.5) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * Power(mws, 1.5) * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (8. * (4 * mws - mzs) * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (Cplx(0, 0.03125) * (el * el * el) * sqrt(mws) * (2 * mws - mzs) * (mws + mzs) * C0Im(mzs, 0, 0, 0, 0, mws)) / (mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (Cplx(0, 0.03125) * (el * el * el) * (mzs * mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) / (Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * Power(mws, 1.5) * (mws + mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (16. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) -
                    (el * el * el * sqrt(mws) * (2 * mws - mzs) * (mws + mzs) * C0Re(mzs, 0, 0, 0, 0, mws)) / (32. * mzs * Power(-mws + mzs, 1.5) * (Pi * Pi)) +
                    (el * el * el * (mzs * mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (32. * Power(mws, 1.5) * Power(-mws + mzs, 1.5) * (Pi * Pi)));

            break;
        }
        break;
    }
}
//Cplx G_SMNLO::result(void) const
Cplx G1(const int ftyp, const int formt, const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)),
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz;

    switch (ftyp)
    {
    case LEP:
        switch (formt)
        {
        case VEC:
            return (-(el * el * el * (2 * mws + mzs) * A0fin(mws)) / (64. * mws * (mws - mzs) * (Pi * Pi)) -
                    (5 * (el * el * el) * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) * A0fin(mzs)) / (128. * mws * (mws - mzs) * mzs * (Pi * Pi)) +
                    (el * el * el * (8 * (mws * mws) - 26 * mws * mzs + 5 * (mzs * mzs) + Cplx(0, 10) * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (256. * mws * (mws - mzs) * (Pi * Pi)) +
                    (5 * (el * el * el) * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (128. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (64. * (mws - mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.0625) * (el * el * el) * mzs * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) / (mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mws * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * (-mws + mzs) * (Pi * Pi)) +
                    (el * el * el * mzs * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (16. * mws * (mws - mzs) * (Pi * Pi)));
            break;
        case AXV:
            return (-(el * el * el * (2 * mws + mzs) * A0fin(mws)) / (64. * mws * (mws - mzs) * (Pi * Pi)) -
                    (5 * (el * el * el) * (4 * mws - 3 * mzs) * A0fin(mzs)) / (128. * mws * (mws - mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.00390625) * (el * el * el) * mzs * (Cplx(0, 1) * (10 * mws + 3 * mzs) + 10 * (4 * mws - 3 * mzs) * B0Im(mzs, 0, 0))) / (mws * (mws - mzs) * (Pi * Pi)) +
                    (5 * (el * el * el) * (4 * mws - 3 * mzs) * mzs * B0Refin(mzs, 0, 0)) / (128. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (64. * (mws - mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.0625) * (el * el * el) * (4 * mws - 3 * mzs) * (mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) / (mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mws * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * (-mws + mzs) * (Pi * Pi)) +
                    (el * el * el * (4 * mws - 3 * mzs) * (mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (16. * mws * (mws - mzs) * (Pi * Pi)));
            break;
        }
        break;
    case UQU:
    case CQU:
        switch (formt)
        {
        case VEC:
            return ((el * el * el * A0fin(mws)) / (48. * (mws - mzs) * (Pi * Pi)) + (5 * (el * el * el) * (32 * (mws * mws) - 40 * mws * mzs + 17 * (mzs * mzs)) * A0fin(mzs)) / (1728. * mws * (mws - mzs) * mzs * (Pi * Pi)) + (el * el * el * (-32 * (mws * mws) + 238 * mws * mzs - 17 * (mzs * mzs) - Cplx(0, 2) * (142 * (mws * mws) - 227 * mws * mzs + 85 * (mzs * mzs)) * B0Im(mzs, 0, 0))) / (3456. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (-142 * mws + 85 * mzs) * B0Refin(mzs, 0, 0)) / (1728. * mws * (Pi * Pi)) +
                    (el * el * el * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (64. * (-mws + mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * ((mws + mzs) * (mws + mzs)) * C0Im(mzs, 0, 0, 0, 0, mws)) / ((mws - mzs) * (Pi * Pi)) -
                    (Cplx(0, 0.004629629629629629) * (el * el * el) * mzs * (32 * (mws * mws) - 40 * mws * mzs + 17 * (mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mws * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) / (96. * (mws - mzs) * (Pi * Pi)) -
                    (el * el * el * mzs * (32 * (mws * mws) - 40 * mws * mzs + 17 * (mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (216. * mws * (mws - mzs) * (Pi * Pi)));
            break;
        case AXV:
            return ((el * el * el * A0fin(mws)) / (48. * (mws - mzs) * (Pi * Pi)) + (5 * (el * el * el) * (8 * mws - 5 * mzs) * A0fin(mzs)) / (576. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (mzs * (58 * mws + 5 * mzs) + Cplx(0, 2) * (6 * (mws * mws) - 31 * mws * mzs + 25 * (mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (1152. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (6 * mws - 25 * mzs) * B0Refin(mzs, 0, 0)) / (576. * mws * (Pi * Pi)) +
                    (el * el * el * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (64. * (-mws + mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.010416666666666666) * (el * el * el) * ((mws + mzs) * (mws + mzs)) * C0Im(mzs, 0, 0, 0, 0, mws)) / ((mws - mzs) * (Pi * Pi)) -
                    (Cplx(0, 0.013888888888888888) * (el * el * el) * (8 * mws - 5 * mzs) * (mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) / (mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mws * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) / (96. * (mws - mzs) * (Pi * Pi)) -
                    (el * el * el * (8 * mws - 5 * mzs) * (mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (72. * mws * (mws - mzs) * (Pi * Pi)));

            break;
        }
        break;
    case DQU:
    case SQU:
        switch (formt)
        {
        case VEC:
            return (-(el * el * el * (2 * mws - 3 * mzs) * A0fin(mws)) / (192. * mws * (mws - mzs) * (Pi * Pi)) -
                    (5 * (el * el * el) * (8 * (mws * mws) - 4 * mws * mzs + 5 * (mzs * mzs)) * A0fin(mzs)) / (3456. * mws * (mws - mzs) * mzs * (Pi * Pi)) +
                    (el * el * el * (8 * (mws * mws) - 418 * mws * mzs + 5 * (mzs * mzs) - Cplx(0, 2) * (32 * (mws * mws) + 128 * mws * mzs - 25 * (mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (6912. * mws * (mws - mzs) * (Pi * Pi)) -
                    (el * el * el * (32 * (mws * mws) + 128 * mws * mzs - 25 * (mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (3456. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (64. * (mws - mzs) * (Pi * Pi)) -
                    (Cplx(0, 0.020833333333333332) * (el * el * el) * ((mws + mzs) * (mws + mzs)) * C0Im(mzs, 0, 0, 0, 0, mws)) / ((mws - mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.0023148148148148147) * (el * el * el) * mzs * (8 * (mws * mws) - 4 * mws * mzs + 5 * (mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mws * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * (-mws + mzs) * (Pi * Pi)) -
                    (el * el * el * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) / (48. * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mzs * (8 * (mws * mws) - 4 * mws * mzs + 5 * (mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (432. * mws * (mws - mzs) * (Pi * Pi)));
            break;
        case AXV:
            return (-(el * el * el * (2 * mws - 3 * mzs) * A0fin(mws)) / (192. * mws * (mws - mzs) * (Pi * Pi)) -
                    (5 * (el * el * el) * (4 * mws - mzs) * A0fin(mzs)) / (1152. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (-(mzs * (134 * mws + mzs)) - Cplx(0, 2) * (24 * (mws * mws) + 16 * mws * mzs + 5 * (mzs * mzs)) * B0Im(mzs, 0, 0))) /
                        (2304. * mws * (mws - mzs) * (Pi * Pi)) -
                    (el * el * el * (24 * (mws * mws) + 16 * mws * mzs + 5 * (mzs * mzs)) * B0Refin(mzs, 0, 0)) /
                        (1152. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (64. * (mws - mzs) * (Pi * Pi)) -
                    (Cplx(0, 0.020833333333333332) * (el * el * el) * ((mws + mzs) * (mws + mzs)) * C0Im(mzs, 0, 0, 0, 0, mws)) / ((mws - mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.006944444444444444) * (el * el * el) * (4 * mws - mzs) * (mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) / (mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mws * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * (-mws + mzs) * (Pi * Pi)) -
                    (el * el * el * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) / (48. * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (4 * mws - mzs) * (mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (144. * mws * (mws - mzs) * (Pi * Pi)));

            break;
        }
        break;
    case NEU:
        switch (formt)

        {
        case VEC:
            return ((el * el * el * mzs * A0fin(mws)) / (32. * mws * (-mws + mzs) * (Pi * Pi)) +
                    (el * el * el * (4 * mzs + Cplx(0, 1) * (2 * mws + 3 * mzs) * B0Im(mzs, 0, 0))) / (64. * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (2 * mws + 3 * mzs) * B0Refin(mzs, 0, 0)) / (64. * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (64. * (-mws + mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.03125) * (el * el * el) * ((mws + mzs) * (mws + mzs)) * C0Im(mzs, 0, 0, 0, 0, mws)) / ((mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mws * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) / (32. * (mws - mzs) * (Pi * Pi)));
            break;
        case AXV:
            return ((el * el * el * mzs * A0fin(mws)) / (32. * mws * (-mws + mzs) * (Pi * Pi)) +
                    (el * el * el * (4 * mzs + Cplx(0, 1) * (2 * mws + 3 * mzs) * B0Im(mzs, 0, 0))) / (64. * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (2 * mws + 3 * mzs) * B0Refin(mzs, 0, 0)) / (64. * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (2 * mws + mzs) * B0Refin(mzs, mws, mws)) / (64. * (-mws + mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.03125) * (el * el * el) * ((mws + mzs) * (mws + mzs)) * C0Im(mzs, 0, 0, 0, 0, mws)) / ((mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mws * (mws + 2 * mzs) * C0Re(0, mzs, 0, 0, mws, mws)) / (32. * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * ((mws + mzs) * (mws + mzs)) * C0Re(mzs, 0, 0, 0, 0, mws)) / (32. * (mws - mzs) * (Pi * Pi)));
            break;
        }
        break;

    case BQU:
        switch (formt)

        {
        case VEC:
            return ((el * el * el * (2 * (mts * mts * mts) + mts * mts * mzs - 6 * mts * mws * (mws + mzs) + 4 * (mws * mws) * (mws + 2 * mzs)) * A0fin(mts)) /
                        (384. * ((mts - mws) * (mts - mws)) * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (-2 * (mts * mts * mts) - 4 * (mws * mws * mws) + mts * mws * (6 * mws - 19 * mzs) + 10 * (mts * mts) * mzs + 6 * (mws * mws) * mzs) * A0fin(mws)) /
                        (384. * ((mts - mws) * (mts - mws)) * mws * (mws - mzs) * (Pi * Pi)) -
                    (5 * (el * el * el) * (8 * (mws * mws) - 4 * mws * mzs + 5 * (mzs * mzs)) * A0fin(mzs)) / (3456. * mws * (mws - mzs) * mzs * (Pi * Pi)) -
                    (el * el * el * (81 * (mts * mts) * mzs + mts * (-8 * (mws * mws) + 283 * mws * mzs - 5 * (mzs * mzs)) + mws * (8 * (mws * mws) - 418 * mws * mzs + 5 * (mzs * mzs)) - Cplx(0, 10) * (mts - mws) * (8 * (mws * mws) - 4 * mws * mzs + 5 * (mzs * mzs)) * B0Im(mzs, 0, 0))) / (6912. * (mts - mws) * mws * (mws - mzs) * (Pi * Pi)) +
                    (5 * (el * el * el) * (8 * (mws * mws) - 4 * mws * mzs + 5 * (mzs * mzs)) * B0Refin(mzs, 0, 0)) / (3456. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (2 * (mts * mts) + mts * (2 * mws + mzs) - 2 * mws * (2 * mws + 3 * mzs)) * B0Refin(mzs, mts, mts)) / (192. * mws * (mws - mzs) * (Pi * Pi)) -
                    (el * el * el * (2 * (mts * mts) + mts * (2 * mws + mzs) - 2 * mws * (2 * mws + mzs)) * B0Refin(mzs, mws, mws)) / (128. * mws * (mws - mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.0023148148148148147) * (el * el * el) * mzs * (8 * (mws * mws) - 4 * mws * mzs + 5 * (mzs * mzs)) * C0Im(mzs, 0, 0, 0, 0, mzs)) /
                        (mws * (mws - mzs) * (Pi * Pi)) -
                    (el * el * el * (mts * mts * mts + mts * mts * mzs + 2 * (mws * mws) * (mws + 2 * mzs) - mts * mws * (3 * mws + 2 * mzs)) *
                     C0Re(0, mzs, 0, mts, mws, mws)) /
                        (64. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mzs * (8 * (mws * mws) - 4 * mws * mzs + 5 * (mzs * mzs)) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (432. * mws * (mws - mzs) * (Pi * Pi)) -
                    (el * el * el * (mts * mts * mts + mts * mts * mzs + 2 * mws * ((mws + mzs) * (mws + mzs)) - mts * mws * (3 * mws + 2 * mzs)) * C0Re(mzs, 0, 0, mts, mts, mws)) /
                        (96. * mws * (mws - mzs) * (Pi * Pi)));
            break;
        case AXV:
            return ((el * el * el * (2 * (mts * mts * mts) + mts * mts * mzs - 6 * mts * mws * (mws + mzs) + 4 * (mws * mws) * (mws + 2 * mzs)) * A0fin(mts)) /
                        (384. * ((mts - mws) * (mts - mws)) * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (-2 * (mts * mts * mts) - 4 * (mws * mws * mws) + mts * mws * (6 * mws - 19 * mzs) + 10 * (mts * mts) * mzs + 6 * (mws * mws) * mzs) * A0fin(mws)) /
                        (384. * ((mts - mws) * (mts - mws)) * mws * (mws - mzs) * (Pi * Pi)) -
                    (5 * (el * el * el) * (4 * mws - mzs) * A0fin(mzs)) / (1152. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * mzs * (-27 * (mts * mts) - mts * (89 * mws + mzs) + mws * (134 * mws + mzs) + Cplx(0, 10) * (mts - mws) * (4 * mws - mzs) * B0Im(mzs, 0, 0))) /
                        (2304. * (mts - mws) * mws * (mws - mzs) * (Pi * Pi)) +
                    (5 * (el * el * el) * (4 * mws - mzs) * mzs * B0Refin(mzs, 0, 0)) / (1152. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (2 * (mts * mts) + mts * (2 * mws + mzs) - 2 * mws * (2 * mws + 3 * mzs)) * B0Refin(mzs, mts, mts)) / (192. * mws * (mws - mzs) * (Pi * Pi)) -
                    (el * el * el * (2 * (mts * mts) + mts * (2 * mws + mzs) - 2 * mws * (2 * mws + mzs)) * B0Refin(mzs, mws, mws)) / (128. * mws * (mws - mzs) * (Pi * Pi)) +
                    (Cplx(0, 0.006944444444444444) * (el * el * el) * (4 * mws - mzs) * (mzs * mzs) * C0Im(mzs, 0, 0, 0, 0, mzs)) / (mws * (mws - mzs) * (Pi * Pi)) -
                    (el * el * el * (mts * mts * mts + mts * mts * mzs + 2 * (mws * mws) * (mws + 2 * mzs) - mts * mws * (3 * mws + 2 * mzs)) * C0Re(0, mzs, 0, mts, mws, mws)) /
                        (64. * mws * (mws - mzs) * (Pi * Pi)) +
                    (el * el * el * (4 * mws - mzs) * (mzs * mzs) * C0Re(mzs, 0, 0, 0, 0, mzs)) / (144. * mws * (mws - mzs) * (Pi * Pi)) -
                    (el * el * el * (mts * mts * mts + mts * mts * mzs + 2 * mws * ((mws + mzs) * (mws + mzs)) - mts * mws * (3 * mws + 2 * mzs)) * C0Re(mzs, 0, 0, mts, mts, mws)) /
                        (96. * mws * (mws - mzs) * (Pi * Pi)));

            break;
        }
        break;
    }
}

//one-loop gauge-boson mediated boxes contribution
Cplx B1(const int ftyp, const int iff, const int off, const int GB, const inval &input, double s)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)), cost = 0.5,
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           gamz = ival->get(GamZ),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz,
           elq = el * el * el * el,
           swq = sw * sw * sw * sw;

    switch (ftyp)
    {
    case LEP:
        if (iff == VEC && off == VEC && GB == GG)

            //gamma mgamma
            return (-((el * el * el * el * (2 * cost * (Pi * Pi) + 2 * (cost * cost * cost) * (Pi * Pi) + B0Refin(-((1 - cost) * mzs) / 2., 0, 0) - B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + clog1((1 - cost) / (1 + cost)) - 2 * (cost * cost) * clog1((1 - cost) / (1 + cost)) + cost * cost * cost * cost * clog1((1 - cost) / (1 + cost)) - clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) + 2 * (cost * cost) * clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) - cost * cost * cost * cost * clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) + 2 * cost * (clog1(mzs) * clog1(mzs)) + 2 * (cost * cost * cost) * (clog1(mzs) * clog1(mzs)) + 2 * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * cost * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + 2 * (cost * cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * cost * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * (cost * cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + 2 * (cost * cost) * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + cost * cost * Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, -4) + Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, -2) + cost * cost * cost * cost * Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, -2) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, -1) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1) + 2 * cost * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - 2 * (cost * cost * cost) * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, 2) + cost * cost * cost * cost * Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, 2) + cost * cost * Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, 4))) /
                      (8. * ((-1 + cost) * (-1 + cost)) * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi))));
        if (iff == VEC && off == VEC && GB == GZ)
            //gamma Z
            return (-(elq * ((cost * ((4 * mws - 3 * mzs) * (4 * mws - 3 * mzs)) - mzs * mzs) * (Pi * Pi) +
                             6 * (cost * ((4 * mws - 3 * mzs) * (4 * mws - 3 * mzs)) - mzs * mzs) * (B0Refin(mzs, mzs, 0) + B0Im(mzs, mzs, 0) * Cplx(0, 1)) +
                             6 * (1 - cost) * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) *
                                 (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) -
                             24 * (1 + cost) * (mws - mzs) * (2 * mws - mzs) * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) -
                             3 * ((-1 + cost) * (-1 + cost)) * mzs * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) *
                                 (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) + C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) * Cplx(0, 1)) +
                             12 * ((1 + cost) * (1 + cost)) * (mws - mzs) * (2 * mws - mzs) * mzs *
                                 (C0Re(-((1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) + C0Im(-((1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) * Cplx(0, 1)) -
                             24 * (1 + cost) * (mws - mzs) * (2 * mws - mzs) * li2((-1 + cost) / (1 + cost)) -
                             6 * (-1 + cost) * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) * li2((1 + cost) / (-1 + cost)) +
                             6 * (-1 + cost) * ((4 * mws - 3 * mzs) * (4 * mws - 3 * mzs)) * log(-2 / (-1 + cost)) -
                             3 * (cost * ((4 * mws - 3 * mzs) * (4 * mws - 3 * mzs)) - mzs * mzs) * (log(2 / (1 + cost)) * log(2 / (1 + cost))) -
                             3 * (-1 + cost) * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) * (log((1 - cost) / (1 + cost)) * log((1 - cost) / (1 + cost))) +
                             6 * log(2 / (1 + cost)) * ((1 + cost) * ((4 * mws - 3 * mzs) * (4 * mws - 3 * mzs)) + (-1 + cost) * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) * log((1 - cost) / (1 + cost)) - 2 * (cost * ((4 * mws - 3 * mzs) * (4 * mws - 3 * mzs)) - mzs * mzs) * log(1 - s / (mzs - I * gamz * mz))) -
                             3 * (1 - cost) * log((1 - cost) / (1 + cost)) * ((1 + cost) * ((4 * mws - 3 * mzs) * (4 * mws - 3 * mzs)) + 4 * (8 * (mws * mws) - 12 * mws * mzs + 5 * (mzs * mzs)) * log(1 - s / (mzs - I * gamz * mz))))) /
                    (192. * (-1 + cost) * (1 + cost) * mws * (mws - mzs) * mzs * (Pi * Pi)));
        //ZZ + WW
        if (iff == VEC && off == VEC && GB == ZZWW)
            return (-(-(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-((1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) + cost * ((1 - 4 * (sw * sw) + 8 * swq) * (1 - 4 * (sw * sw) + 8 * swq))) *
                       (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (256. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) +
                        B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (32. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) +
                              B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (1 - 8 * (sw * sw) - 32 * (sw * sw * sw * sw * sw * sw) + 32 * (sw * sw * sw * sw * sw * sw * sw * sw) + 24 * swq) * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (256. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) *
                       (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (2 - 16 * (sw * sw) + 32 * swq - 3 * cost * ((1 - 4 * (sw * sw) + 8 * swq) * (1 - 4 * (sw * sw) + 8 * swq)) + cost * cost * cost * ((1 - 4 * (sw * sw) + 8 * swq) * (1 - 4 * (sw * sw) + 8 * swq))) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) *
                           (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq *
                       ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) +
                        C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (64. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * (1 - 8 * (sw * sw) - 32 * (sw * sw * sw * sw * sw * sw) + 32 * (sw * sw * sw * sw * sw * sw * sw * sw) + 24 * swq) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) +
                        C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) *
                       (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (2 - 16 * (sw * sw) + 32 * swq - 3 * cost * ((1 - 4 * (sw * sw) + 8 * swq) * (1 - 4 * (sw * sw) + 8 * swq)) + cost * cost * cost * ((1 - 4 * (sw * sw) + 8 * swq) * (1 - 4 * (sw * sw) + 8 * swq))) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) *
                           (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq *
                       ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(((-1 - cost) * mzs) / 2., 0, 0, 0, 0, mzs) +
                        C0Im(((-1 - cost) * mzs) / 2., 0, 0, 0, 0, mzs) * Cplx(0, 1))) /
                          (64. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * elq * (2 * mws + cost * mzs) *
                       (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) +
                        C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) * Cplx(0, 1))) /
                          (64. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * (2 + cost) * elq * (1 - 8 * (sw * sw) - 32 * (sw * sw * sw * sw * sw * sw) + 32 * (sw * sw * sw * sw * sw * sw * sw * sw) + 24 * swq) *
                       (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) +
                        C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * ((8 * (-3 + cost) * (-3 + cost * cost) * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw)))) / ((-1 + cost) * (-1 + cost)) + ((1 + cost) * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw)))) / (sw * sw * sw * sw) + 16 * (1 + cost) * swq) *
                       D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * (cw * cw * cw * cw) * (Pi * Pi)) +
                      (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost - cost * cost + cost * cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) /
                          (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (-9 - 3 * cost + 3 * (cost * cost) + cost * cost * cost - 8 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 64 * (-5 + cost * (-2 + cost * (2 + cost))) * (sw * sw * sw * sw * sw * sw) + 64 * (-5 + cost * (-2 + cost * (2 + cost))) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (2 + cost) * (-7 + cost + 2 * (cost * cost)) * swq) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == VEC && off == AXV && GB == ZZWW)
            //ZZ + WW
            return (-(-(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) / (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (256. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (1 - 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (256. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((2 + cost) * elq * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (256. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-1 + cost) * (2 + cost) * elq * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((2 + cost) * elq * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (256. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * elq * (2 * mws + cost * mzs) * (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) + C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) * Cplx(0, 1))) / (64. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-1 + cost) * (2 + cost) * elq * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) + C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((1 + cost) * elq * mzs * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost) * (1 + cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) /
                          (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((3 + cost) * (-3 + cost * cost) * elq * mzs * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));

        if (iff == VEC && off == AXV && GB == GZ)
            //gamma Z
            return (-(elq * (-1 + 4 * (sw * sw)) * (-(Pi * Pi) + cost * (Pi * Pi) - 6 * B0Refin(mz * mz, 0, mz * mz) + 6 * cost * B0Refin(mz * mz, 0, mz * mz) + 6 * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) - 6 * cost * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) - 3 * (mz * mz) * C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + 6 * cost * (mz * mz) * C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) - 3 * (cost * cost) * (mz * mz) * C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + (-1 + cost) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -6) + mz * mz * C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, -3) + cost * cost * (mz * mz) * C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, -3) + (-1 + cost) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, 6) + cost * (mz * mz) * C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, 6) + 6 * li2((1 + cost) / (-1 + cost)) - 6 * cost * li2((1 + cost) / (-1 + cost)) - 6 * log(-2 / (-1 + cost)) + 6 * cost * log(-2 / (-1 + cost)) + 6 * log(2 / (1 + cost)) + 6 * cost * log(2 / (1 + cost)) - 3 * log((1 - cost) / (1 + cost)) + 3 * (cost * cost) * log((1 - cost) / (1 + cost)) + 3 * (log(mz * mz) * log(mz * mz)) - 3 * cost * (log(mz * mz) * log(mz * mz)) - 6 * log(mz * mz) * log(-((-1 + cost) * (mz * mz)) / 2.) + 6 * cost * log(mz * mz) * log(-((-1 + cost) * (mz * mz)) / 2.) + 3 * (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) - 3 * cost * (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) + 12 * log(mz * mz) * log(1 - s / (mzs - I * gamz * mz)) - 12 * cost * log(mz * mz) * log(1 - s / (mzs - I * gamz * mz)) - 12 * log(-((-1 + cost) * (mz * mz)) / 2.) * log(1 - s / (mzs - I * gamz * mz)) + 12 * cost * log(-((-1 + cost) * (mz * mz)) / 2.) * log(1 - s / (mzs - I * gamz * mz)))) / (192. * (-1 + cost * cost) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        if (iff == AXV && off == VEC && GB == GZ)
            //gamma Z
            return (-(elq * (-1 + 4 * (sw * sw)) * (-(Pi * Pi) + cost * (Pi * Pi) - 6 * B0Refin(mz * mz, 0, mz * mz) + 6 * cost * B0Refin(mz * mz, 0, mz * mz) + 6 * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) - 6 * cost * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) - 3 * (mz * mz) * C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + 6 * cost * (mz * mz) * C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) - 3 * (cost * cost) * (mz * mz) * C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + (-1 + cost) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -6) + mz * mz * C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, -3) + cost * cost * (mz * mz) * C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, -3) + (-1 + cost) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, 6) + cost * (mz * mz) * C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, 6) + 6 * li2((1 + cost) / (-1 + cost)) - 6 * cost * li2((1 + cost) / (-1 + cost)) - 6 * log(-2 / (-1 + cost)) + 6 * cost * log(-2 / (-1 + cost)) + 6 * log(2 / (1 + cost)) + 6 * cost * log(2 / (1 + cost)) - 3 * log((1 - cost) / (1 + cost)) + 3 * (cost * cost) * log((1 - cost) / (1 + cost)) + 3 * (log(mz * mz) * log(mz * mz)) - 3 * cost * (log(mz * mz) * log(mz * mz)) - 6 * log(mz * mz) * log(-((-1 + cost) * (mz * mz)) / 2.) + 6 * cost * log(mz * mz) * log(-((-1 + cost) * (mz * mz)) / 2.) + 3 * (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) - 3 * cost * (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) + 12 * log(mz * mz) * log(1 - s / (mzs - I * gamz * mz)) - 12 * cost * log(mz * mz) * log(1 - s / (mzs - I * gamz * mz)) - 12 * log(-((-1 + cost) * (mz * mz)) / 2.) * log(1 - s / (mzs - I * gamz * mz)) + 12 * cost * log(-((-1 + cost) * (mz * mz)) / 2.) * log(1 - s / (mzs - I * gamz * mz)))) / (192. * (-1 + cost * cost) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        // ZZ + WW
        if (iff == AXV && off == VEC && GB == ZZWW)
            return (-(-(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) / (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (256. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (1 - 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (256. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((2 + cost) * elq * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (256. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-1 + cost) * (2 + cost) * elq * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((2 + cost) * elq * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (256. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * elq * (2 * mws + cost * mzs) * (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) + C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) * Cplx(0, 1))) / (64. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-1 + cost) * (2 + cost) * elq * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) + C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((1 + cost) * elq * mzs * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost) * (1 + cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) /
                          (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((3 + cost) * (-3 + cost * cost) * elq * mzs * (-1 + 8 * (sw * sw + 4 * (sw * sw * sw * sw * sw * sw) - 3 * swq)) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));

        if (iff == AXV && off == AXV && GB == ZZWW)
            //ZZ + WW
            return (-(-(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (cost * ((1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) - (1 - 4 * (sw * sw) + 8 * swq) * (1 - 4 * (sw * sw) + 8 * swq)) *
                       (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (256. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (32. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (1 - 8 * (sw * sw) - 32 * (sw * sw * sw * sw * sw * sw) + 32 * (sw * sw * sw * sw * sw * sw * sw * sw) + 24 * swq) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (256. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) *
                       (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (-3 * cost * ((1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) + cost * cost * cost * ((1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) + 2 * ((1 - 4 * (sw * sw) + 8 * swq) * (1 - 4 * (sw * sw) + 8 * swq))) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) +
                        C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (64. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * (1 - 8 * (sw * sw) - 32 * (sw * sw * sw * sw * sw * sw) + 32 * (sw * sw * sw * sw * sw * sw * sw * sw) + 24 * swq) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) +
                        C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) *
                       (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (-3 * cost * ((1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) + cost * cost * cost * ((1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) + 2 * ((1 - 4 * (sw * sw) + 8 * swq) * (1 - 4 * (sw * sw) + 8 * swq))) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(((-1 - cost) * mzs) / 2., 0, 0, 0, 0, mzs) +
                        C0Im(((-1 - cost) * mzs) / 2., 0, 0, 0, 0, mzs) * Cplx(0, 1))) /
                          (64. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * elq * (2 * mws + cost * mzs) *
                       (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) +
                        C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) * Cplx(0, 1))) /
                          (64. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * (2 + cost) * elq * (1 - 8 * (sw * sw) - 32 * (sw * sw * sw * sw * sw * sw) + 32 * (sw * sw * sw * sw * sw * sw * sw * sw) + 24 * swq) *
                       (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) +
                        C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mzs) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * ((-1 + cost) * (-1 + cost) * (1 + cost) - 8 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 64 * (-4 + cost + cost * cost) * (sw * sw * sw * sw * sw * sw) + 64 * (-4 + cost + cost * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (-3 + cost * cost * cost) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost) * (1 + cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) /
                          (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (-9 - 3 * cost + 3 * (cost * cost) + cost * cost * cost - 8 * (3 + cost) * (-3 + cost * cost) * (sw * sw) + 64 * (4 + cost - cost * cost) * (sw * sw * sw * sw * sw * sw) + 64 * (-4 + (-1 + cost) * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (-13 + cost * (-4 + cost * (4 + cost))) * swq) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == AXV && off == AXV && GB == GG)
            //gamma-gamma
            return (-((el * el * el * el * (-4 * (cost * cost) * (Pi * Pi) + B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Refin(-((1 + cost) * mzs) / 2., 0, 0) - 4 * (cost * cost) * (clog1(mzs) * clog1(mzs)) - 2 * cost * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + 4 * (cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + 2 * cost * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + 4 * (cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) - cost * cost * cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1) - 2 * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) + 2 * (cost * cost) * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)))) /
                      (8. * ((-1 + cost) * (-1 + cost)) * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi))));
        //gamma-Z
        if (iff == AXV && off == AXV && GB == GZ)
            return ((elq * ((-1 + cost) * (1 + cost) * (Pi * Pi) * (cost - (1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) +
                            6 * (-1 + cost) * (1 + cost) * (cost - (1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) * B0Refin(mz * mz, 0, mz * mz) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (1 - 4 * (sw * sw) + 8 * (sw * sw * sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) +
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + 2 * (sw * sw)) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (-1 + cost) * (-1 + cost * cost) * (1 - 4 * (sw * sw) + 8 * (sw * sw * sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -6) +
                            (-1 + cost) * (-1 + cost) * (-1 + cost * cost) * (mz * mz) * (1 - 4 * (sw * sw) + 8 * (sw * sw * sw * sw)) * Cplx(0, -3) *
                                (C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, -1)) -
                            12 * ((1 + cost) * (1 + cost)) * (1 - cost * cost) * (mz * mz) * (sw * sw) * (1 - 2 * (sw * sw)) *
                                (C0Re(-((1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + C0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, 1)) +
                            (-1 + cost) * (1 + cost) * (cost - (1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, 6) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + 2 * (sw * sw)) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 24) +
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + 2 * (sw * sw)) * li2((-1 + cost) / (1 + cost)) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (1 - 4 * (sw * sw) + 8 * (sw * sw * sw * sw)) * li2((1 + cost) / (-1 + cost)) +
                            6 * (-1 + cost) * (-1 + cost * cost) * log(-2 / (-1 + cost)) + 6 * (1 + cost) * (-1 + cost * cost) * log(2 / (1 + cost)) +
                            3 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * log((1 - cost) / (1 + cost)) -
                            3 * (-1 + cost) * (1 + cost) * (cost - (1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) * (log(mz * mz) * log(mz * mz)) -
                            3 * (-1 + cost) * (-1 + cost * cost) * (1 - 4 * (sw * sw) + 8 * (sw * sw * sw * sw)) *
                                (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) +
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + 2 * (sw * sw)) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * (-((-1 + cost) * (1 - 4 * (sw * sw) + 8 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.)) + 4 * (1 + cost) * (sw * sw) * (-1 + 2 * (sw * sw)) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            12 * ((1 - cost) * (1 + cost) * (cost - (1 - 4 * (sw * sw)) * (1 - 4 * (sw * sw))) * log(mz * mz) + (-1 + cost) * (-1 + cost * cost) * (1 - 4 * (sw * sw) + 8 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) - 4 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + 2 * (sw * sw)) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (192. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        else
            return (0);
        break;

    case NEU:
        if (iff == VEC && off == VEC && GB == ZZWW)
            //ZZ+WW
            return (-(-(elq * (-1 + cost - 4 * (-1 + cost) * (sw * sw) + 8 * cost * swq) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (256. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (256. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (2 - 8 * (sw * sw) + cost * (-3 + cost * cost) * (1 - 4 * (sw * sw) + 8 * swq)) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (128. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (2 - 8 * (sw * sw) + cost * (-3 + cost * cost) * (1 - 4 * (sw * sw) + 8 * swq)) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (128. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((1 + cost) * elq * mzs * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) / (256. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * ((-1 + cost) * (-1 + cost) * (1 + cost) - 4 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) + 8 * (5 + cost * (-2 + (-2 + cost) * cost)) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (-4 + 4 * cost + ((3 + cost) * (-3 + cost * cost) * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw)))) / ((1 + cost) * (1 + cost) * (sw * sw * sw * sw))) * D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * (cw * cw * cw * cw) * (Pi * Pi))));
        if (iff == VEC && off == AXV && GB == ZZWW)
            return (-(-(elq * (-1 + cost - 4 * (-1 + cost) * (sw * sw) + 8 * cost * swq) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (256. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (256. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (2 - 8 * (sw * sw) + cost * (-3 + cost * cost) * (1 - 4 * (sw * sw) + 8 * swq)) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (128. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (2 - 8 * (sw * sw) + cost * (-3 + cost * cost) * (1 - 4 * (sw * sw) + 8 * swq)) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (128. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((1 + cost) * elq * mzs * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) / (256. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * ((-1 + cost) * (-1 + cost) * (1 + cost) - 4 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) + 8 * (5 + cost * (-2 + (-2 + cost) * cost)) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (-4 + 4 * cost + ((3 + cost) * (-3 + cost * cost) * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw)))) / ((1 + cost) * (1 + cost) * (sw * sw * sw * sw))) * D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * (cw * cw * cw * cw) * (Pi * Pi))));
        if (iff == AXV && off == VEC && GB == ZZWW)
            return (-((elq * (1 - cost + 4 * (-1 + cost) * (sw * sw) + 8 * swq) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (256. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (256. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (2 - 3 * cost + cost * cost * cost - 4 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) + 16 * swq) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (128. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (2 - 3 * cost + cost * cost * cost - 4 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) + 16 * swq) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (128. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((1 + cost) * elq * mzs * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) / (256. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * ((-1 + cost) * (-1 + cost) * (1 + cost) - 4 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) + 8 * (-4 + cost + cost * cost) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (-9 - 3 * cost + 3 * (cost * cost) + cost * cost * cost - 4 * (3 + cost) * (-3 + cost * cost) * (sw * sw) + 8 * (-4 + (-1 + cost) * cost) * swq) * D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == AXV && off == AXV && GB == ZZWW)
            return (-((elq * (1 - cost + 4 * (-1 + cost) * (sw * sw) + 8 * swq) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (256. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (256. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (2 - 3 * cost + cost * cost * cost - 4 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) + 16 * swq) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (128. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (2 - 3 * cost + cost * cost * cost - 4 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) + 16 * swq) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (256. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (128. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * ((1 - 2 * (sw * sw)) * (1 - 2 * (sw * sw))) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (512. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * ((-1 + cost) * (-1 + cost) * (1 + cost) - 4 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) + 8 * (-4 + cost + cost * cost) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (-9 - 3 * cost + 3 * (cost * cost) + cost * cost * cost - 4 * (3 + cost) * (-3 + cost * cost) * (sw * sw) + 8 * (-4 + (-1 + cost) * cost) * swq) * D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (1024. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((1 + cost) * elq * mzs * D0(-((1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) / (256. * (Pi * Pi) * (sw * sw * sw * sw))));
        else
            return (0);

    case CQU:
    case UQU:
        if (iff == VEC && off == VEC && GB == ZZWW)
            return (-((elq * (9 - 9 * cost + 60 * (-1 + cost) * (sw * sw) + 320 * cost * (sw * sw * sw * sw * sw * sw) - 256 * cost * (sw * sw * sw * sw * sw * sw * sw * sw) + 8 * (12 - 25 * cost) * swq) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (13 - 40 * (sw * sw) + 32 * swq) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (576. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (9 + 4 * (sw * sw) * (-15 + 37 * (sw * sw) + 32 * (sw * sw * sw * sw * sw * sw) - 40 * swq)) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 320 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw) + 256 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 8 * (24 + 25 * cost * (-3 + cost * cost)) * swq) * (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (13 - 40 * (sw * sw) + 32 * swq) *
                       (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (1152. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 4 * (sw * sw) * (-15 + 37 * (sw * sw) + 32 * (sw * sw * sw * sw * sw * sw) - 40 * swq)) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 320 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw) + 256 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 8 * (24 + 25 * cost * (-3 + cost * cost)) * swq) * (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (13 - 40 * (sw * sw) + 32 * swq) *
                       (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (1152. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 4 * (sw * sw) * (-15 + 37 * (sw * sw) + 32 * (sw * sw * sw * sw * sw * sw) - 40 * swq)) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((1 + cost) * elq * mzs * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) / (256. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 320 * (5 + cost * (-2 + (-2 + cost) * cost)) * (sw * sw * sw * sw * sw * sw) + 256 * (5 + cost * (-2 + (-2 + cost) * cost)) * (sw * sw * sw * sw * sw * sw * sw * sw) + 8 * (77 + cost * (-38 + cost * (-38 + 25 * cost))) * swq) *
                       D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 60 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 320 * (-5 + cost * (-2 + cost * (2 + cost))) * (sw * sw * sw * sw * sw * sw) + 256 * (-5 + cost * (-2 + cost * (2 + cost))) * (sw * sw * sw * sw * sw * sw * sw * sw) + 8 * (-173 + cost * (-62 + cost * (62 + 25 * cost))) * swq) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == VEC && off == VEC && GB == GG)

            return ((4. / 9. )*(-((el * el * el * el * (2 * cost * (Pi * Pi) + 2 * (cost * cost * cost) * (Pi * Pi) + B0Refin(-((1 - cost) * mzs) / 2., 0, 0) - B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + clog1((1 - cost) / (1 + cost)) - 2 * (cost * cost) * clog1((1 - cost) / (1 + cost)) + cost * cost * cost * cost * clog1((1 - cost) / (1 + cost)) - clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) + 2 * (cost * cost) * clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) - cost * cost * cost * cost * clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) + 2 * cost * (clog1(mzs) * clog1(mzs)) + 2 * (cost * cost * cost) * (clog1(mzs) * clog1(mzs)) + 2 * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * cost * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + 2 * (cost * cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * cost * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * (cost * cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + 2 * (cost * cost) * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + cost * cost * Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, -4) + Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, -2) + cost * cost * cost * cost * Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, -2) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, -1) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1) + 2 * cost * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - 2 * (cost * cost * cost) * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, 2) + cost * cost * cost * cost * Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, 2) + cost * cost * Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, 4))) /
                      (8. * ((-1 + cost) * (-1 + cost)) * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi)))));
        if (iff == VEC && off == VEC && GB == GZ)

            return ((elq * ((-1 + cost) * (1 + cost) * (Pi * Pi) * (-3 + cost * (3 - 20 * (sw * sw) + 32 * (sw * sw * sw * sw))) +
                            6 * (-1 + cost) * (1 + cost) * (-3 + cost * (3 - 20 * (sw * sw) + 32 * (sw * sw * sw * sw))) * B0Refin(mz * mz, 0, mz * mz) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) -
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -12) +
                            (-1 + cost) * (-1 + cost * cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -6) +
                            (-1 + cost) * (-1 + cost) * (-1 + cost * cost) * (mz * mz) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * Cplx(0, -3) *
                                (C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, -1)) +
                            6 * ((1 + cost) * (1 + cost)) * (-1 + cost * cost) * (mz * mz) * (sw * sw) * (-5 + 8 * (sw * sw)) *
                                (C0Re(-((1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + C0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, 1)) +
                            (-1 + cost) * (1 + cost) * (-3 + cost * (3 - 20 * (sw * sw) + 32 * (sw * sw * sw * sw))) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, 6) -
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * li2((-1 + cost) / (1 + cost)) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * li2((1 + cost) / (-1 + cost)) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 20 * (sw * sw) + 32 * (sw * sw * sw * sw)) * log(-2 / (-1 + cost)) +
                            6 * (1 + cost) * (-1 + cost * cost) * (3 - 20 * (sw * sw) + 32 * (sw * sw * sw * sw)) * log(2 / (1 + cost)) +
                            3 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * (-1 + 4 * (sw * sw)) * (-3 + 8 * (sw * sw)) * log((1 - cost) / (1 + cost)) -
                            3 * (-1 + cost) * (1 + cost) * (-3 + cost * (3 - 20 * (sw * sw) + 32 * (sw * sw * sw * sw))) * (log(mz * mz) * log(mz * mz)) -
                            3 * (-1 + cost) * (-1 + cost * cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) *
                                (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) -
                            6 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * (-((-1 + cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.)) + 2 * (1 + cost) * (sw * sw) * (5 - 8 * (sw * sw)) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * ((-1 + cost) * (1 + cost) * (-3 + cost * (3 - 20 * (sw * sw) + 32 * (sw * sw * sw * sw))) * log(mz * mz) - (-1 + cost) * (-1 + cost * cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) - 2 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (864. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        if (iff == VEC && off == AXV && GB == ZZWW)

            return (-((elq * (9 - 9 * cost + 60 * (-1 + cost) * (sw * sw) + 64 * (-2 + 3 * cost) * (sw * sw * sw * sw * sw * sw) + 8 * (16 - 21 * cost) * swq) *
                       (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (5 - 8 * (sw * sw)) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (576. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (9 - 60 * (sw * sw) - 160 * (sw * sw * sw * sw * sw * sw) + 148 * swq) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 64 * (4 - 9 * cost + 3 * (cost * cost * cost)) * (sw * sw * sw * sw * sw * sw) + 8 * (32 + 21 * cost * (-3 + cost * cost)) * swq) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (-5 + 8 * (sw * sw)) * (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) / (1152. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) -
                      ((-1 + cost) * (2 + cost) * elq * (-9 + 60 * (sw * sw) + 160 * (sw * sw * sw * sw * sw * sw) - 148 * swq) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 64 * (4 - 9 * cost + 3 * (cost * cost * cost)) * (sw * sw * sw * sw * sw * sw) + 8 * (32 + 21 * cost * (-3 + cost * cost)) * swq) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (-5 + 8 * (sw * sw)) * (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) / (1152. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) -
                      ((-1 + cost) * (2 + cost) * elq * (-9 + 60 * (sw * sw) + 160 * (sw * sw * sw * sw * sw * sw) - 148 * swq) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((1 + cost) * elq * mzs * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) / (256. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 64 * (7 + (-2 + cost) * cost * (2 + 3 * cost)) * (sw * sw * sw * sw * sw * sw) + 8 * (41 + cost * (-26 + cost * (-26 + 21 * cost))) * swq) *
                       D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 60 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 64 * (-23 + cost * (-8 + cost * (8 + 3 * cost))) * (sw * sw * sw * sw * sw * sw) + 8 * (-169 + cost * (-58 + cost * (58 + 21 * cost))) * swq) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == VEC && off == AXV && GB == GZ)
            return ((elq * (-((-1 + cost) * (1 + cost) * (Pi * Pi) * (3 - 8 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw)))) -
                            6 * (-1 + cost) * (1 + cost) * (3 - 8 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * B0Refin(mz * mz, 0, mz * mz) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 10 * (sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) +
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (-1 + cost) * (1 + cost) * (3 - 8 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, -6) +
                            3 * ((-1 + cost) * (-1 + cost)) * (-1 + cost * cost) * (mz * mz) * (-3 + 10 * (sw * sw)) *
                                (C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, 1)) +
                            6 * ((1 + cost) * (1 + cost)) * (1 - cost * cost) * (mz * mz) * (sw * sw) *
                                (C0Re(-((1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + C0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, 1)) +
                            (-1 + cost) * (-1 + cost * cost) * (-3 + 10 * (sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 6) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 12) +
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * li2((-1 + cost) / (1 + cost)) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 10 * (sw * sw)) * li2((1 + cost) / (-1 + cost)) -
                            18 * (-1 + cost) * (-1 + cost * cost) * (-1 + 4 * (sw * sw)) * log(-2 / (-1 + cost)) -
                            18 * (1 + cost) * (-1 + cost * cost) * (-1 + 4 * (sw * sw)) * log(2 / (1 + cost)) -
                            9 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * (-1 + 4 * (sw * sw)) * log((1 - cost) / (1 + cost)) +
                            3 * (-1 + cost) * (1 + cost) * (3 - 8 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * (log(mz * mz) * log(mz * mz)) +
                            3 * (-1 + cost) * (-1 + cost * cost) * (-3 + 10 * (sw * sw)) * (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * ((-1 + cost) * (-3 + 10 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) + 2 * (1 + cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * (-((-1 + cost) * (1 + cost) * (3 - 8 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * log(mz * mz)) + (-1 + cost) * (-1 + cost * cost) * (-3 + 10 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) + 2 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (864. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        if (iff == AXV && off == VEC && GB == ZZWW)

            return (-((elq * (9 - 9 * cost + 60 * (-1 + cost) * (sw * sw) + 64 * (-3 + 2 * cost) * (sw * sw * sw * sw * sw * sw) + 8 * (21 - 16 * cost) * swq) *
                       (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (-5 + 8 * (sw * sw)) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (576. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (9 - 60 * (sw * sw) - 160 * (sw * sw * sw * sw * sw * sw) + 148 * swq) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 128 * (3 - 3 * cost + cost * cost * cost) * (sw * sw * sw * sw * sw * sw) + 16 * (21 + 8 * cost * (-3 + cost * cost)) * swq) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (-5 + 8 * (sw * sw)) * (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) / (1152. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) -
                      ((-1 + cost) * (2 + cost) * elq * (-9 + 60 * (sw * sw) + 160 * (sw * sw * sw * sw * sw * sw) - 148 * swq) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 128 * (3 - 3 * cost + cost * cost * cost) * (sw * sw * sw * sw * sw * sw) + 16 * (21 + 8 * cost * (-3 + cost * cost)) * swq) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (-5 + 8 * (sw * sw)) * (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) / (1152. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) -
                      ((-1 + cost) * (2 + cost) * elq * (-9 + 60 * (sw * sw) + 160 * (sw * sw * sw * sw * sw * sw) - 148 * swq) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((1 + cost) * elq * mzs * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) / (256. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) + 64 * (2 + cost + cost * cost - 2 * (cost * cost * cost)) * (sw * sw * sw * sw * sw * sw) + 8 * (-4 + cost * (-11 + cost * (-11 + 16 * cost))) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 60 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 64 * (-22 + cost * (-7 + cost * (7 + 2 * cost))) * (sw * sw * sw * sw * sw * sw) + 8 * (-164 + cost * (-53 + cost * (53 + 16 * cost))) * swq) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == AXV && off == VEC && GB == GZ)

            return ((elq * (-((-1 + cost) * (1 + cost) * (Pi * Pi) * (3 - 12 * (sw * sw) + cost * (-3 + 8 * (sw * sw)))) -
                            6 * (-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 8 * (sw * sw))) * B0Refin(mz * mz, 0, mz * mz) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 10 * (sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) -
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -12) +
                            (-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 8 * (sw * sw))) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, -6) +
                            3 * ((-1 + cost) * (-1 + cost)) * (-1 + cost * cost) * (mz * mz) * (-3 + 10 * (sw * sw)) *
                                (C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, 1)) +
                            6 * ((1 + cost) * (1 + cost)) * (-1 + cost * cost) * (mz * mz) * (sw * sw) *
                                (C0Re(-((1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + C0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, 1)) +
                            (-1 + cost) * (-1 + cost * cost) * (-3 + 10 * (sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 6) -
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * li2((-1 + cost) / (1 + cost)) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 10 * (sw * sw)) * li2((1 + cost) / (-1 + cost)) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * log(-2 / (-1 + cost)) -
                            6 * (1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * log(2 / (1 + cost)) -
                            3 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * log((1 - cost) / (1 + cost)) +
                            3 * (-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 8 * (sw * sw))) * (log(mz * mz) * log(mz * mz)) +
                            3 * (-1 + cost) * (-1 + cost * cost) * (-3 + 10 * (sw * sw)) * (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) -
                            6 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * ((-1 + cost) * (-3 + 10 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) - 2 * (1 + cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * (-((-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 8 * (sw * sw))) * log(mz * mz)) + (-1 + cost) * (-1 + cost * cost) * (-3 + 10 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) - 2 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (864. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        if (iff == AXV && off == AXV && GB == ZZWW)

            return (-((elq * (9 - 9 * cost + 60 * (-1 + cost) * (sw * sw) - 320 * (sw * sw * sw * sw * sw * sw) + 256 * (sw * sw * sw * sw * sw * sw * sw * sw) + 8 * (25 - 12 * cost) * swq) *
                       (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (13 - 40 * (sw * sw) + 32 * swq) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (576. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (9 + 4 * (sw * sw) * (-15 + 37 * (sw * sw) + 32 * (sw * sw * sw * sw * sw * sw) - 40 * swq)) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 640 * (sw * sw * sw * sw * sw * sw) + 512 * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (25 + 6 * cost * (-3 + cost * cost)) * swq) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (13 - 40 * (sw * sw) + 32 * swq) *
                       (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (1152. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 4 * (sw * sw) * (-15 + 37 * (sw * sw) + 32 * (sw * sw * sw * sw * sw * sw) - 40 * swq)) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) / (128. * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 640 * (sw * sw * sw * sw * sw * sw) + 512 * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (25 + 6 * cost * (-3 + cost * cost)) * swq) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (13 - 40 * (sw * sw) + 32 * swq) *
                       (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (1152. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 4 * (sw * sw) * (-15 + 37 * (sw * sw) + 32 * (sw * sw * sw * sw * sw * sw) - 40 * swq)) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 60 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 320 * (-4 + cost + cost * cost) * (sw * sw * sw * sw * sw * sw) + 256 * (-4 + cost + cost * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 8 * (-40 + cost + cost * cost + 12 * (cost * cost * cost)) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 60 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 320 * (-4 + (-1 + cost) * cost) * (sw * sw * sw * sw * sw * sw) + 256 * (-4 + (-1 + cost) * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 8 * (-160 + cost * (-49 + cost * (49 + 12 * cost))) * swq) * D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((1 + cost) * elq * mzs * D0(-((1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) / (256. * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == AXV && off == AXV && GB == GG)

            return ((4. /9.)*(-((el * el * el * el * (-4 * (cost * cost) * (Pi * Pi) + B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Refin(-((1 + cost) * mzs) / 2., 0, 0) - 4 * (cost * cost) * (clog1(mzs) * clog1(mzs)) - 2 * cost * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + 4 * (cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + 2 * cost * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + 4 * (cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) - cost * cost * cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1) - 2 * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) + 2 * (cost * cost) * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)))) /
                      (8. * ((-1 + cost) * (-1 + cost)) * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi)))));
        if (iff == AXV && off == AXV && GB == GZ)

            return ((elq * ((-1 + cost) * (1 + cost) * (Pi * Pi) * (-3 + 3 * cost + 20 * (sw * sw) - 32 * (sw * sw * sw * sw)) +
                            6 * (-1 + cost) * (1 + cost) * (-3 + 3 * cost + 20 * (sw * sw) - 32 * (sw * sw * sw * sw)) * B0Refin(mz * mz, 0, mz * mz) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) +
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (-1 + cost) * (-1 + cost * cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -6) +
                            (-1 + cost) * (-1 + cost) * (-1 + cost * cost) * (mz * mz) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * Cplx(0, -3) *
                                (C0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + C0Re(((-1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, -1)) -
                            6 * ((1 + cost) * (1 + cost)) * (1 - cost * cost) * (mz * mz) * (sw * sw) * (5 - 8 * (sw * sw)) *
                                (C0Re(-((1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) + C0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0, 0, 0, mz * mz) * Cplx(0, 1)) +
                            (-1 + cost) * (1 + cost) * (-3 + 3 * cost + 20 * (sw * sw) - 32 * (sw * sw * sw * sw)) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, 6) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 12) +
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * li2((-1 + cost) / (1 + cost)) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * li2((1 + cost) / (-1 + cost)) +
                            18 * (-1 + cost) * (-1 + cost * cost) * log(-2 / (-1 + cost)) + 18 * (1 + cost) * (-1 + cost * cost) * log(2 / (1 + cost)) +
                            9 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * log((1 - cost) / (1 + cost)) -
                            3 * (-1 + cost) * (1 + cost) * (-3 + 3 * cost + 20 * (sw * sw) - 32 * (sw * sw * sw * sw)) * (log(mz * mz) * log(mz * mz)) -
                            3 * (-1 + cost) * (-1 + cost * cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) *
                                (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * (-((-1 + cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.)) + 2 * (1 + cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * ((-1 + cost) * (1 + cost) * (-3 + 3 * cost + 20 * (sw * sw) - 32 * (sw * sw * sw * sw)) * log(mz * mz) - (-1 + cost) * (-1 + cost * cost) * (3 - 10 * (sw * sw) + 16 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) + 2 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-5 + 8 * (sw * sw)) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (864. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        else
            return (0);

    case DQU:
    case SQU:
        if (iff == VEC && off == VEC && GB == ZZWW)

            return (-(-(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) / (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 - 9 * cost + 48 * (-1 + cost) * (sw * sw) + 128 * cost * (sw * sw * sw * sw * sw * sw) - 64 * cost * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (3 - 8 * cost) * swq) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (5 - 8 * (sw * sw) + 4 * swq) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (288. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 128 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw) + 64 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 32 * (3 + 4 * cost * (-3 + cost * cost)) * swq) * (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (5 - 8 * (sw * sw) + 4 * swq) *
                       (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (576. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 128 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw) + 64 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 32 * (3 + 4 * cost * (-3 + cost * cost)) * swq) * (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (5 - 8 * (sw * sw) + 4 * swq) *
                       (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (576. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * elq * (2 * mws + cost * mzs) * (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) + C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) * Cplx(0, 1))) / (64. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 128 * (5 + cost * (-2 + (-2 + cost) * cost)) * (sw * sw * sw * sw * sw * sw) + 64 * (5 + cost * (-2 + (-2 + cost) * cost)) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (28 + cost * (-13 + cost * (-13 + 8 * cost))) * swq) *
                       D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost - cost * cost + cost * cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs) * D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) /
                          (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 48 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 128 * (-5 + cost * (-2 + cost * (2 + cost))) * (sw * sw * sw * sw * sw * sw) + 64 * (-5 + cost * (-2 + cost * (2 + cost))) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (-52 + cost * (-19 + cost * (19 + 8 * cost))) * swq) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == VEC && off == VEC && GB == GG)

            return ((1. / 9. )*(-((el * el * el * el * (2 * cost * (Pi * Pi) + 2 * (cost * cost * cost) * (Pi * Pi) + B0Refin(-((1 - cost) * mzs) / 2., 0, 0) - B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + clog1((1 - cost) / (1 + cost)) - 2 * (cost * cost) * clog1((1 - cost) / (1 + cost)) + cost * cost * cost * cost * clog1((1 - cost) / (1 + cost)) - clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) + 2 * (cost * cost) * clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) - cost * cost * cost * cost * clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) + 2 * cost * (clog1(mzs) * clog1(mzs)) + 2 * (cost * cost * cost) * (clog1(mzs) * clog1(mzs)) + 2 * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * cost * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + 2 * (cost * cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * cost * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * (cost * cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + 2 * (cost * cost) * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + cost * cost * Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, -4) + Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, -2) + cost * cost * cost * cost * Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, -2) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, -1) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1) + 2 * cost * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - 2 * (cost * cost * cost) * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, 2) + cost * cost * cost * cost * Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, 2) + cost * cost * Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, 4))) /
                      (8. * ((-1 + cost) * (-1 + cost)) * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi)))));
        if (iff == VEC && off == VEC && GB == GZ)

            return ((elq * ((-1 + cost) * (1 + cost) * (Pi * Pi) * (-3 + cost * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw))) +
                            6 * (-1 + cost) * (1 + cost) * (-3 + cost * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw))) * B0Refin(mz * mz, 0, mz * mz) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) -
                            48 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -48) +
                            (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -6) +
                            (-1 + cost) * (-1 + cost) * (-1 + cost * cost) * (mz * mz) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * Cplx(0, -3) *
                                (C0Im(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Re(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, -1)) +
                            24 * ((1 + cost) * (1 + cost)) * (-1 + cost * cost) * (mz * mz) * (sw * sw) * (-1 + sw * sw) *
                                (C0Re(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Im(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, 1)) +
                            (-1 + cost) * (1 + cost) * (-3 + cost * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw))) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, 6) -
                            48 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * li2((-1 + cost) / (1 + cost)) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * li2((1 + cost) / (-1 + cost)) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw)) * log(-2 / (-1 + cost)) +
                            6 * (1 + cost) * (-1 + cost * cost) * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw)) * log(2 / (1 + cost)) +
                            3 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * (-3 + 4 * (sw * sw)) * (-1 + 4 * (sw * sw)) * log(-((-1 + cost) / (1 + cost))) -
                            3 * (-1 + cost) * (1 + cost) * (-3 + cost * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw))) * (log(mz * mz) * log(mz * mz)) -
                            3 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) *
                                (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) -
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * (-((-1 + cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.)) - 8 * (1 + cost) * (sw * sw) * (-1 + sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * ((-1 + cost) * (1 + cost) * (-3 + cost * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw))) * log(mz * mz) - (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) - 8 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (1728. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        if (iff == VEC && off == AXV && GB == ZZWW)
            return (-(-(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) / (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 - 9 * cost + 48 * (-1 + cost) * (sw * sw) + 32 * (-1 + 3 * cost) * (sw * sw * sw * sw * sw * sw) + 8 * (7 - 15 * cost) * swq) *
                       (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-1 + sw * sw) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (72. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 - 48 * (sw * sw) - 64 * (sw * sw * sw * sw * sw * sw) + 88 * swq) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 32 * (2 - 9 * cost + 3 * (cost * cost * cost)) * (sw * sw * sw * sw * sw * sw) + 8 * (14 + 15 * cost * (-3 + cost * cost)) * swq) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (-1 + sw * sw) * (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) / (144. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) -
                      ((-1 + cost) * (2 + cost) * elq * (-3 + 8 * (sw * sw)) * (3 - 8 * (sw * sw) + 8 * swq) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 32 * (2 - 9 * cost + 3 * (cost * cost * cost)) * (sw * sw * sw * sw * sw * sw) + 8 * (14 + 15 * cost * (-3 + cost * cost)) * swq) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (-1 + sw * sw) * (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) / (144. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) -
                      ((-1 + cost) * (2 + cost) * elq * (-3 + 8 * (sw * sw)) * (3 - 8 * (sw * sw) + 8 * swq) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * elq * (2 * mws + cost * mzs) * (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) + C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) * Cplx(0, 1))) / (64. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 32 * (11 + cost * (-5 + cost * (-5 + 3 * cost))) * (sw * sw * sw * sw * sw * sw) + 8 * (47 + cost * (-23 + cost * (-23 + 15 * cost))) * swq) *
                       D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost) * (1 + cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) /
                          (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 48 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 32 * (-19 + cost * (-7 + cost * (7 + 3 * cost))) * (sw * sw * sw * sw * sw * sw) + 8 * (-103 + cost * (-37 + cost * (37 + 15 * cost))) * swq) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == VEC && off == AXV && GB == GZ)
            return ((elq * (-((-1 + cost) * (1 + cost) * (Pi * Pi) * (3 - 4 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw)))) -
                            6 * (-1 + cost) * (1 + cost) * (3 - 4 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * B0Refin(mz * mz, 0, mz * mz) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) +
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (-1 + cost) * (1 + cost) * (3 - 4 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, -6) +
                            3 * ((-1 + cost) * (-1 + cost)) * (-1 + cost * cost) * (mz * mz) * (-3 + 8 * (sw * sw)) *
                                (C0Re(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Im(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, 1)) +
                            12 * ((1 + cost) * (1 + cost)) * (1 - cost * cost) * (mz * mz) * (sw * sw) *
                                (C0Re(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Im(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, 1)) +
                            (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 6) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 24) +
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * li2((-1 + cost) / (1 + cost)) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * li2((1 + cost) / (-1 + cost)) -
                            18 * (-1 + cost) * (-1 + cost * cost) * (-1 + 4 * (sw * sw)) * log(-2 / (-1 + cost)) -
                            18 * (1 + cost) * (-1 + cost * cost) * (-1 + 4 * (sw * sw)) * log(2 / (1 + cost)) -
                            9 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * (-1 + 4 * (sw * sw)) * log(-((-1 + cost) / (1 + cost))) +
                            3 * (-1 + cost) * (1 + cost) * (3 - 4 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * (log(mz * mz) * log(mz * mz)) +
                            3 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) +
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * ((-1 + cost) * (-3 + 8 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) + 4 * (1 + cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * (-((-1 + cost) * (1 + cost) * (3 - 4 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * log(mz * mz)) + (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) + 4 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (1728. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        if (iff == AXV && off == VEC && GB == ZZWW)
            return (-(-(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) / (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 - 9 * cost + 48 * (-1 + cost) * (sw * sw) + 32 * (-3 + cost) * (sw * sw * sw * sw * sw * sw) + 8 * (15 - 7 * cost) * swq) *
                       (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (-1 + sw * sw) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (72. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 - 48 * (sw * sw) - 64 * (sw * sw * sw * sw * sw * sw) + 88 * swq) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 32 * (6 - 3 * cost + cost * cost * cost) * (sw * sw * sw * sw * sw * sw) + 8 * (30 + 7 * cost * (-3 + cost * cost)) * swq) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (-1 + sw * sw) * (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) / (144. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) -
                      ((-1 + cost) * (2 + cost) * elq * (-3 + 8 * (sw * sw)) * (3 - 8 * (sw * sw) + 8 * swq) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 32 * (6 - 3 * cost + cost * cost * cost) * (sw * sw * sw * sw * sw * sw) + 8 * (30 + 7 * cost * (-3 + cost * cost)) * swq) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (-1 + sw * sw) * (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) / (144. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) -
                      ((-1 + cost) * (2 + cost) * elq * (-3 + 8 * (sw * sw)) * (3 - 8 * (sw * sw) + 8 * swq) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * elq * (2 * mws + cost * mzs) * (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) + C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) * Cplx(0, 1))) / (64. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 32 * (-7 + cost + cost * cost + cost * cost * cost) * (sw * sw * sw * sw * sw * sw) + 8 * (-25 + cost + cost * cost + 7 * (cost * cost * cost)) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost) * (1 + cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) /
                          (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 48 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 32 * (-17 + cost * (-5 + cost * (5 + cost))) * (sw * sw * sw * sw * sw * sw) + 8 * (-95 + cost * (-29 + cost * (29 + 7 * cost))) * swq) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == AXV && off == VEC && GB == GZ)
            return ((elq * (-((-1 + cost) * (1 + cost) * (Pi * Pi) * (3 - 12 * (sw * sw) + cost * (-3 + 4 * (sw * sw)))) -
                            6 * (-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 4 * (sw * sw))) * B0Refin(mz * mz, 0, mz * mz) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) -
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -24) +
                            (-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 4 * (sw * sw))) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, -6) +
                            3 * ((-1 + cost) * (-1 + cost)) * (-1 + cost * cost) * (mz * mz) * (-3 + 8 * (sw * sw)) *
                                (C0Re(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Im(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, 1)) +
                            12 * ((1 + cost) * (1 + cost)) * (-1 + cost * cost) * (mz * mz) * (sw * sw) *
                                (C0Re(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Im(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, 1)) +
                            (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 6) -
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * li2((-1 + cost) / (1 + cost)) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * li2((1 + cost) / (-1 + cost)) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 4 * (sw * sw)) * log(-2 / (-1 + cost)) -
                            6 * (1 + cost) * (-1 + cost * cost) * (-3 + 4 * (sw * sw)) * log(2 / (1 + cost)) -
                            3 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * (-3 + 4 * (sw * sw)) * log(-((-1 + cost) / (1 + cost))) +
                            3 * (-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 4 * (sw * sw))) * (log(mz * mz) * log(mz * mz)) +
                            3 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) -
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * ((-1 + cost) * (-3 + 8 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) - 4 * (1 + cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * (-((-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 4 * (sw * sw))) * log(mz * mz)) + (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) - 4 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (1728. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        if (iff == AXV && off == AXV && GB == ZZWW)
            return (-(-(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) / (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 - 9 * cost + 48 * (-1 + cost) * (sw * sw) - 128 * (sw * sw * sw * sw * sw * sw) + 64 * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (8 - 3 * cost) * swq) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (5 - 8 * (sw * sw) + 4 * swq) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (288. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, 0, mzs, mws, 0, mws) + C0Im(0, 0, mzs, mws, 0, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 256 * (sw * sw * sw * sw * sw * sw) + 128 * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (16 - 9 * cost + 3 * (cost * cost * cost)) * swq) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (5 - 8 * (sw * sw) + 4 * swq) *
                       (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (576. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (-4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 256 * (sw * sw * sw * sw * sw * sw) + 128 * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (16 - 9 * cost + 3 * (cost * cost * cost)) * swq) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (5 - 8 * (sw * sw) + 4 * swq) *
                       (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (576. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * elq * (2 * mws + cost * mzs) * (C0Re(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) + C0Im(((-1 + cost) * mzs) / 2., 0, 0, 0, 0, mws) * Cplx(0, 1))) / (64. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 128 * (-4 + cost + cost * cost) * (sw * sw * sw * sw * sw * sw) + 64 * (-4 + cost + cost * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (-17 + cost * (2 + cost * (2 + 3 * cost))) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost) * (1 + cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mws, mws)) /
                          (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 48 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 128 * (-4 + (-1 + cost) * cost) * (sw * sw * sw * sw * sw * sw) + 64 * (-4 + (-1 + cost) * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (-47 + cost * (-14 + cost * (14 + 3 * cost))) * swq) * D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == AXV && off == AXV && GB == GG)
            return ((1. /9.)*(-((el * el * el * el * (-4 * (cost * cost) * (Pi * Pi) + B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Refin(-((1 + cost) * mzs) / 2., 0, 0) - 4 * (cost * cost) * (clog1(mzs) * clog1(mzs)) - 2 * cost * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + 4 * (cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + 2 * cost * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + 4 * (cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) - cost * cost * cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1) - 2 * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) + 2 * (cost * cost) * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)))) /
                      (8. * ((-1 + cost) * (-1 + cost)) * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi)))));
        if (iff == AXV && off == AXV && GB == GZ)
            return ((elq * ((-1 + cost) * (1 + cost) * (Pi * Pi) * (-3 + 3 * cost + 16 * (sw * sw) - 16 * (sw * sw * sw * sw)) +
                            6 * (-1 + cost) * (1 + cost) * (-3 + 3 * cost + 16 * (sw * sw) - 16 * (sw * sw * sw * sw)) * B0Refin(mz * mz, 0, mz * mz) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) +
                            48 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -6) +
                            (-1 + cost) * (-1 + cost) * (-1 + cost * cost) * (mz * mz) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * Cplx(0, -3) *
                                (C0Im(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Re(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, -1)) +
                            (1 + cost) * (1 + cost) * (-1 + cost * cost) * (mz * mz) * (sw * sw) * (-1 + sw * sw) * Cplx(0, -24) *
                                (C0Im(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Re(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, -1)) +
                            (-1 + cost) * (1 + cost) * (-3 + 3 * cost + 16 * (sw * sw) - 16 * (sw * sw * sw * sw)) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, 6) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 48) +
                            48 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * li2((-1 + cost) / (1 + cost)) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * li2((1 + cost) / (-1 + cost)) +
                            18 * (-1 + cost) * (-1 + cost * cost) * log(-2 / (-1 + cost)) + 18 * (1 + cost) * (-1 + cost * cost) * log(2 / (1 + cost)) +
                            9 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * log(-((-1 + cost) / (1 + cost))) -
                            3 * (-1 + cost) * (1 + cost) * (-3 + 3 * cost + 16 * (sw * sw) - 16 * (sw * sw * sw * sw)) * (log(mz * mz) * log(mz * mz)) -
                            3 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) *
                                (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) +
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * (-((-1 + cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.)) + 8 * (1 + cost) * (sw * sw) * (-1 + sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * ((-1 + cost) * (1 + cost) * (-3 + 3 * cost + 16 * (sw * sw) - 16 * (sw * sw * sw * sw)) * log(mz * mz) - (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) + 8 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (1728. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        else
            return (0);
    case BQU:
        if (iff == VEC && off == VEC && GB == ZZWW)

            return (-(-(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) / (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 - 9 * cost + 48 * (-1 + cost) * (sw * sw) + 128 * cost * (sw * sw * sw * sw * sw * sw) - 64 * cost * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (3 - 8 * cost) * swq) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (5 - 8 * (sw * sw) + 4 * swq) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (288. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, mt * mt) + B0Im(((-1 + cost) * mzs) / 2., 0, mt * mt) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (2 * (mt * mt) - 4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, 0, mzs, mws, mt * mt, mws) + C0Im(0, 0, mzs, mws, mt * mt, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 128 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw) + 64 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 32 * (3 + 4 * cost * (-3 + cost * cost)) * swq) * (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (5 - 8 * (sw * sw) + 4 * swq) *
                       (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (576. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * elq * (-(mt * mt) + 2 * mws + cost * mzs) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (2 * (mt * mt) - 4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 128 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw) + 64 * cost * (-3 + cost * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 32 * (3 + 4 * cost * (-3 + cost * cost)) * swq) * (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      ((-2 + cost) * (1 + cost) * elq * (5 - 8 * (sw * sw) + 4 * swq) *
                       (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (576. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * elq * (-(mt * mt) + 2 * mws + cost * mzs) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt, 0) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 128 * (5 + cost * (-2 + (-2 + cost) * cost)) * (sw * sw * sw * sw * sw * sw) + 64 * (5 + cost * (-2 + (-2 + cost) * cost)) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (28 + cost * (-13 + cost * (-13 + 8 * cost))) * swq) *
                       D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 48 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 128 * (-5 + cost * (-2 + cost * (2 + cost))) * (sw * sw * sw * sw * sw * sw) + 64 * (-5 + cost * (-2 + cost * (2 + cost))) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (-52 + cost * (-19 + cost * (19 + 8 * cost))) * swq) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (-4 * (mt * mt * mt * mt) + 4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost - cost * cost + cost * cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs - 2 * (mt * mt) * (2 * (-3 + cost) * mws + (2 - cost + cost * cost) * mzs)) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, mt * mt, mws, mws)) /
                          (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == VEC && off == VEC && GB == GG)

            return (-((1. / 9.) * (el * el * el * el * (2 * cost * (Pi * Pi) + 2 * (cost * cost * cost) * (Pi * Pi) + B0Refin(-((1 - cost) * mzs) / 2., 0, 0) - B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + clog1((1 - cost) / (1 + cost)) - 2 * (cost * cost) * clog1((1 - cost) / (1 + cost)) + cost * cost * cost * cost * clog1((1 - cost) / (1 + cost)) - clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) + 2 * (cost * cost) * clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) - cost * cost * cost * cost * clog1((1 - cost) / (1 + cost)) * clog1(((1 - cost) * (1 + cost) * (mz * mz * mz * mz)) / 4.) + 2 * cost * (clog1(mzs) * clog1(mzs)) + 2 * (cost * cost * cost) * (clog1(mzs) * clog1(mzs)) + 2 * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * cost * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + 2 * (cost * cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * cost * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - 2 * (cost * cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + 2 * (cost * cost) * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + cost * cost * Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, -4) + Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, -2) + cost * cost * cost * cost * Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, -2) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, -1) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1) + 2 * cost * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - 2 * (cost * cost * cost) * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, 2) + cost * cost * cost * cost * Pi * clog1(-((1 + cost) * mzs) / 2.) * Cplx(0, 2) + cost * cost * Pi * clog1(-((1 - cost) * mzs) / 2.) * Cplx(0, 4))) /
                      (8. * ((-1 + cost) * (-1 + cost)) * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi))));
        if (iff == VEC && off == VEC && GB == GZ)

            return ((elq * ((-1 + cost) * (1 + cost) * (Pi * Pi) * (-3 + cost * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw))) +
                            6 * (-1 + cost) * (1 + cost) * (-3 + cost * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw))) * B0Refin(mz * mz, 0, mz * mz) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) -
                            48 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -48) +
                            (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -6) +
                            (-1 + cost) * (-1 + cost) * (-1 + cost * cost) * (mz * mz) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * Cplx(0, -3) *
                                (C0Im(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Re(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, -1)) +
                            24 * ((1 + cost) * (1 + cost)) * (-1 + cost * cost) * (mz * mz) * (sw * sw) * (-1 + sw * sw) *
                                (C0Re(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Im(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, 1)) +
                            (-1 + cost) * (1 + cost) * (-3 + cost * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw))) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, 6) -
                            48 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * li2((-1 + cost) / (1 + cost)) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * li2((1 + cost) / (-1 + cost)) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw)) * log(-2 / (-1 + cost)) +
                            6 * (1 + cost) * (-1 + cost * cost) * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw)) * log(2 / (1 + cost)) +
                            3 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * (-3 + 4 * (sw * sw)) * (-1 + 4 * (sw * sw)) * log(-((-1 + cost) / (1 + cost))) -
                            3 * (-1 + cost) * (1 + cost) * (-3 + cost * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw))) * (log(mz * mz) * log(mz * mz)) -
                            3 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) *
                                (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) -
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * (-((-1 + cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.)) - 8 * (1 + cost) * (sw * sw) * (-1 + sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * ((-1 + cost) * (1 + cost) * (-3 + cost * (3 - 16 * (sw * sw) + 16 * (sw * sw * sw * sw))) * log(mz * mz) - (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) - 8 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (1728. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        if (iff == VEC && off == AXV && GB == ZZWW)
            return (-(
                -(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) / (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (9 - 9 * cost + 48 * (-1 + cost) * (sw * sw) + 32 * (-1 + 3 * cost) * (sw * sw * sw * sw * sw * sw) + 8 * (7 - 15 * cost) * swq) *
                 (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                    (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                (elq * (-1 + sw * sw) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                    (72. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                (elq * (9 - 48 * (sw * sw) - 64 * (sw * sw * sw * sw * sw * sw) + 88 * swq) *
                 (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                    (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, mt * mt) + B0Im(((-1 + cost) * mzs) / 2., 0, mt * mt) * Cplx(0, 1))) /
                    (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                (elq * (2 * (mt * mt) - 4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, 0, mzs, mws, mt * mt, mws) + C0Im(0, 0, mzs, mws, mt * mt, mws) * Cplx(0, 1))) /
                    (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 32 * (2 - 9 * cost + 3 * (cost * cost * cost)) * (sw * sw * sw * sw * sw * sw) + 8 * (14 + 15 * cost * (-3 + cost * cost)) * swq) *
                 (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                    (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                ((-2 + cost) * (1 + cost) * elq * (-1 + sw * sw) * (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) / (144. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                ((-1 + cost) * elq * (-(mt * mt) + 2 * mws + cost * mzs) *
                 (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt) * Cplx(0, 1))) /
                    (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                ((-1 + cost) * (2 + cost) * elq * (-3 + 8 * (sw * sw)) * (3 - 8 * (sw * sw) + 8 * swq) *
                 (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                    (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                (elq * (2 * (mt * mt) - 4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                    (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 32 * (2 - 9 * cost + 3 * (cost * cost * cost)) * (sw * sw * sw * sw * sw * sw) + 8 * (14 + 15 * cost * (-3 + cost * cost)) * swq) *
                 (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                    (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                ((-2 + cost) * (1 + cost) * elq * (-1 + sw * sw) * (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) / (144. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                ((-1 + cost) * elq * (-(mt * mt) + 2 * mws + cost * mzs) *
                 (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt, 0) * Cplx(0, 1))) /
                    (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                ((-1 + cost) * (2 + cost) * elq * (-3 + 8 * (sw * sw)) * (3 - 8 * (sw * sw) + 8 * swq) *
                 (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                    (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 32 * (11 + cost * (-5 + cost * (-5 + 3 * cost))) * (sw * sw * sw * sw * sw * sw) + 8 * (47 + cost * (-23 + cost * (-23 + 15 * cost))) * swq) *
                 D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                    (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 48 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 32 * (-19 + cost * (-7 + cost * (7 + 3 * cost))) * (sw * sw * sw * sw * sw * sw) + 8 * (-103 + cost * (-37 + cost * (37 + 15 * cost))) * swq) *
                 D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                    (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost) * (1 + cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs - 2 * (mt * mt) * (2 * (mt * mt) + 2 * (-3 + cost) * mws + (2 + (-1 + cost) * cost) * mzs)) *
                 D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, mt * mt, mws, mws)) /
                    (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == VEC && off == AXV && GB == GZ)

            return ((elq * (-((-1 + cost) * (1 + cost) * (Pi * Pi) * (3 - 4 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw)))) -
                            6 * (-1 + cost) * (1 + cost) * (3 - 4 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * B0Refin(mz * mz, 0, mz * mz) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) +
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (-1 + cost) * (1 + cost) * (3 - 4 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, -6) +
                            3 * ((-1 + cost) * (-1 + cost)) * (-1 + cost * cost) * (mz * mz) * (-3 + 8 * (sw * sw)) *
                                (C0Re(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Im(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, 1)) +
                            12 * ((1 + cost) * (1 + cost)) * (1 - cost * cost) * (mz * mz) * (sw * sw) *
                                (C0Re(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Im(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, 1)) +
                            (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 6) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 24) +
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * li2((-1 + cost) / (1 + cost)) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * li2((1 + cost) / (-1 + cost)) -
                            18 * (-1 + cost) * (-1 + cost * cost) * (-1 + 4 * (sw * sw)) * log(-2 / (-1 + cost)) -
                            18 * (1 + cost) * (-1 + cost * cost) * (-1 + 4 * (sw * sw)) * log(2 / (1 + cost)) -
                            9 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * (-1 + 4 * (sw * sw)) * log(-((-1 + cost) / (1 + cost))) +
                            3 * (-1 + cost) * (1 + cost) * (3 - 4 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * (log(mz * mz) * log(mz * mz)) +
                            3 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) +
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * ((-1 + cost) * (-3 + 8 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) + 4 * (1 + cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * (-((-1 + cost) * (1 + cost) * (3 - 4 * (sw * sw) + 3 * cost * (-1 + 4 * (sw * sw))) * log(mz * mz)) + (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) + 4 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (1728. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        if (iff == AXV && off == VEC && GB == ZZWW)

            return (-(
                -(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) / (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (9 - 9 * cost + 48 * (-1 + cost) * (sw * sw) + 32 * (-3 + cost) * (sw * sw * sw * sw * sw * sw) + 8 * (15 - 7 * cost) * swq) *
                 (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                    (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (-1 + sw * sw) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                    (72. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                (elq * (9 - 48 * (sw * sw) - 64 * (sw * sw * sw * sw * sw * sw) + 88 * swq) *
                 (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                    (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, mt * mt) + B0Im(((-1 + cost) * mzs) / 2., 0, mt * mt) * Cplx(0, 1))) /
                    (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                (elq * (2 * (mt * mt) - 4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, 0, mzs, mws, mt * mt, mws) + C0Im(0, 0, mzs, mws, mt * mt, mws) * Cplx(0, 1))) /
                    (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 32 * (6 - 3 * cost + cost * cost * cost) * (sw * sw * sw * sw * sw * sw) + 8 * (30 + 7 * cost * (-3 + cost * cost)) * swq) *
                 (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                    (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                ((-2 + cost) * (1 + cost) * elq * (-1 + sw * sw) * (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) / (144. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                ((-1 + cost) * elq * (-(mt * mt) + 2 * mws + cost * mzs) *
                 (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt) * Cplx(0, 1))) /
                    (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                ((-1 + cost) * (2 + cost) * elq * (-3 + 8 * (sw * sw)) * (3 - 8 * (sw * sw) + 8 * swq) *
                 (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                    (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                (elq * (2 * (mt * mt) - 4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                    (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 32 * (6 - 3 * cost + cost * cost * cost) * (sw * sw * sw * sw * sw * sw) + 8 * (30 + 7 * cost * (-3 + cost * cost)) * swq) *
                 (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                    (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                ((-2 + cost) * (1 + cost) * elq * (-1 + sw * sw) * (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) / (144. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                ((-1 + cost) * elq * (-(mt * mt) + 2 * mws + cost * mzs) *
                 (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt, 0) * Cplx(0, 1))) /
                    (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                ((-1 + cost) * (2 + cost) * elq * (-3 + 8 * (sw * sw)) * (3 - 8 * (sw * sw) + 8 * swq) *
                 (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                    (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 32 * (-7 + cost + cost * cost + cost * cost * cost) * (sw * sw * sw * sw * sw * sw) + 8 * (-25 + cost + cost * cost + 7 * (cost * cost * cost)) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                    (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 48 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 32 * (-17 + cost * (-5 + cost * (5 + cost))) * (sw * sw * sw * sw * sw * sw) + 8 * (-95 + cost * (-29 + cost * (29 + 7 * cost))) * swq) *
                 D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                    (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost) * (1 + cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs - 2 * (mt * mt) * (2 * (mt * mt) + 2 * (-3 + cost) * mws + (2 + (-1 + cost) * cost) * mzs)) *
                 D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, mt * mt, mws, mws)) /
                    (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == AXV && off == VEC && GB == GZ)
            return ((elq * (-((-1 + cost) * (1 + cost) * (Pi * Pi) * (3 - 12 * (sw * sw) + cost * (-3 + 4 * (sw * sw)))) -
                            6 * (-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 4 * (sw * sw))) * B0Refin(mz * mz, 0, mz * mz) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) -
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -24) +
                            (-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 4 * (sw * sw))) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, -6) +
                            3 * ((-1 + cost) * (-1 + cost)) * (-1 + cost * cost) * (mz * mz) * (-3 + 8 * (sw * sw)) *
                                (C0Re(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Im(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, 1)) +
                            12 * ((1 + cost) * (1 + cost)) * (-1 + cost * cost) * (mz * mz) * (sw * sw) *
                                (C0Re(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Im(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, 1)) +
                            (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 6) -
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * li2((-1 + cost) / (1 + cost)) +
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * li2((1 + cost) / (-1 + cost)) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (-3 + 4 * (sw * sw)) * log(-2 / (-1 + cost)) -
                            6 * (1 + cost) * (-1 + cost * cost) * (-3 + 4 * (sw * sw)) * log(2 / (1 + cost)) -
                            3 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * (-3 + 4 * (sw * sw)) * log(-((-1 + cost) / (1 + cost))) +
                            3 * (-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 4 * (sw * sw))) * (log(mz * mz) * log(mz * mz)) +
                            3 * (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) -
                            12 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * ((-1 + cost) * (-3 + 8 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) - 4 * (1 + cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * (-((-1 + cost) * (1 + cost) * (3 - 12 * (sw * sw) + cost * (-3 + 4 * (sw * sw))) * log(mz * mz)) + (-1 + cost) * (-1 + cost * cost) * (-3 + 8 * (sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) - 4 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (1728. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        if (iff == AXV && off == AXV && GB == ZZWW)
            return (-(-(elq * (B0Refin(mzs, mws, mws) + B0Im(mzs, mws, mws) * Cplx(0, 1))) / (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 - 9 * cost + 48 * (-1 + cost) * (sw * sw) - 128 * (sw * sw * sw * sw * sw * sw) + 64 * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (8 - 3 * cost) * swq) * (B0Refin(mzs, mzs, mzs) + B0Im(mzs, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * (-1 + cost * cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (5 - 8 * (sw * sw) + 4 * swq) * (B0Refin(((-1 - cost) * mzs) / 2., 0, 0) + B0Im(((-1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (288. * (-1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi)) +
                      (elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (B0Refin(((-1 + cost) * mzs) / 2., 0, 0) + B0Im(((-1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1))) /
                          (2304. * (1 + cost) * (cw * cw * cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (B0Refin(((-1 + cost) * mzs) / 2., 0, mt * mt) + B0Im(((-1 + cost) * mzs) / 2., 0, mt * mt) * Cplx(0, 1))) /
                          (64. * (1 + cost) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (2 * (mt * mt) - 4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, 0, mzs, mws, mt * mt, mws) + C0Im(0, 0, mzs, mws, mt * mt, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 256 * (sw * sw * sw * sw * sw * sw) + 128 * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (16 - 9 * cost + 3 * (cost * cost * cost)) * swq) *
                       (C0Re(0, 0, mzs, mzs, 0, mzs) + C0Im(0, 0, mzs, mzs, 0, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (5 - 8 * (sw * sw) + 4 * swq) *
                       (C0Re(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 - cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (576. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * elq * (-(mt * mt) + 2 * mws + cost * mzs) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (C0Re(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) + C0Im(0, 0, ((-1 + cost) * mzs) / 2., 0, mzs, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) -
                      (elq * (2 * (mt * mt) - 4 * mws + (1 + cost * cost) * mzs) * (C0Re(0, mzs, 0, 0, mws, mws) + C0Im(0, mzs, 0, 0, mws, mws) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (9 * ((-1 + cost) * (-1 + cost)) * (2 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (2 + cost) * (sw * sw) - 256 * (sw * sw * sw * sw * sw * sw) + 128 * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (16 - 9 * cost + 3 * (cost * cost * cost)) * swq) *
                       (C0Re(0, mzs, 0, 0, mzs, mzs) + C0Im(0, mzs, 0, 0, mzs, mzs) * Cplx(0, 1))) /
                          (2304. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-2 + cost) * (1 + cost) * elq * (5 - 8 * (sw * sw) + 4 * swq) *
                       (C0Re(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 - cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (576. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi)) +
                      ((-1 + cost) * elq * (-(mt * mt) + 2 * mws + cost * mzs) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mws, mt * mt, 0) * Cplx(0, 1))) /
                          (128. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      ((-1 + cost) * (2 + cost) * elq * (9 + 8 * (sw * sw) * (-6 + 11 * (sw * sw) + 4 * (sw * sw * sw * sw * sw * sw) - 8 * swq)) *
                       (C0Re(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) + C0Im(0, ((-1 + cost) * mzs) / 2., 0, mzs, 0, 0) * Cplx(0, 1))) /
                          (4608. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * ((-1 + cost) * (-1 + cost)) * (1 + cost) - 48 * ((-1 + cost) * (-1 + cost)) * (1 + cost) * (sw * sw) - 128 * (-4 + cost + cost * cost) * (sw * sw * sw * sw * sw * sw) + 64 * (-4 + cost + cost * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (-17 + cost * (2 + cost * (2 + 3 * cost))) * swq) * D0(((-1 - cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((-1 + cost) * (-1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * mzs * (9 * (3 + cost) * (-3 + cost * cost) - 48 * (3 + cost) * (-3 + cost * cost) * (sw * sw) - 128 * (-4 + (-1 + cost) * cost) * (sw * sw * sw * sw * sw * sw) + 64 * (-4 + (-1 + cost) * cost) * (sw * sw * sw * sw * sw * sw * sw * sw) + 16 * (-47 + cost * (-14 + cost * (14 + 3 * cost))) * swq) * D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, 0, mzs, mzs)) /
                          (9216. * ((1 + cost) * (1 + cost)) * (cw * cw * cw * cw) * (Pi * Pi) * (sw * sw * sw * sw)) +
                      (elq * (4 * (-3 + cost) * (mw * mw * mw * mw) + (-1 + cost) * (1 + cost * cost) * (mz * mz * mz * mz) + 4 * ((-1 + cost) * (-1 + cost)) * mws * mzs - 2 * (mt * mt) * (2 * (mt * mt) + 2 * (-3 + cost) * mws + (2 + (-1 + cost) * cost) * mzs)) *
                       D0(((-1 + cost) * mzs) / 2., 0, mzs, 0, 0, 0, 0, mt * mt, mws, mws)) /
                          (256. * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi) * (sw * sw * sw * sw))));
        if (iff == AXV && off == AXV && GB == GG)

            return ((1. /9.)*(-((el * el * el * el * (-4 * (cost * cost) * (Pi * Pi) + B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Refin(-((1 + cost) * mzs) / 2., 0, 0) - 4 * (cost * cost) * (clog1(mzs) * clog1(mzs)) - 2 * cost * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + 4 * (cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) - 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 - cost) * mzs) / 2.) + cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + cost * cost * cost * (clog1(-((1 - cost) * mzs) / 2.) * clog1(-((1 - cost) * mzs) / 2.)) + 2 * cost * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + 4 * (cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) + 2 * (cost * cost * cost) * clog1(mzs) * clog1(-((1 + cost) * mzs) / 2.) - cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) - 2 * (cost * cost) * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) - cost * cost * cost * (clog1(-((1 + cost) * mzs) / 2.) * clog1(-((1 + cost) * mzs) / 2.)) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1) - 2 * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) + 2 * (cost * cost) * (B0Refin(mzs, 0, 0) + B0Im(mzs, 0, 0) * Cplx(0, 1)) - cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * cost * cost * (B0Refin(-((1 - cost) * mzs) / 2., 0, 0) + B0Im(-((1 - cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) + cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)) - cost * cost * cost * (B0Refin(-((1 + cost) * mzs) / 2., 0, 0) + B0Im(-((1 + cost) * mzs) / 2., 0, 0) * Cplx(0, 1)))) /
                      (8. * ((-1 + cost) * (-1 + cost)) * ((1 + cost) * (1 + cost)) * (mz * mz) * (Pi * Pi)))));
        if (iff == AXV && off == AXV && GB == GZ)

            return ((elq * ((-1 + cost) * (1 + cost) * (Pi * Pi) * (-3 + 3 * cost + 16 * (sw * sw) - 16 * (sw * sw * sw * sw)) +
                            6 * (-1 + cost) * (1 + cost) * (-3 + 3 * cost + 16 * (sw * sw) - 16 * (sw * sw * sw * sw)) * B0Refin(mz * mz, 0, mz * mz) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * B0Refin(((-1 + cost) * (mz * mz)) / 2., 0, 0) +
                            48 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * B0Refin(-((1 + cost) * (mz * mz)) / 2., 0, 0) +
                            (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * B0Im(((-1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, -6) +
                            (-1 + cost) * (-1 + cost) * (-1 + cost * cost) * (mz * mz) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * Cplx(0, -3) *
                                (C0Im(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Re(0, 0, ((-1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, -1)) +
                            (1 + cost) * (1 + cost) * (-1 + cost * cost) * (mz * mz) * (sw * sw) * (-1 + sw * sw) * Cplx(0, -24) *
                                (C0Im(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) + C0Re(0, 0, -((1 + cost) * (mz * mz)) / 2., 0, mz * mz, 0) * Cplx(0, -1)) +
                            (-1 + cost) * (1 + cost) * (-3 + 3 * cost + 16 * (sw * sw) - 16 * (sw * sw * sw * sw)) * B0Im(mz * mz, 0, mz * mz) * Cplx(0, 6) +
                            (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * B0Im(-((1 + cost) * (mz * mz)) / 2., 0, 0) * Cplx(0, 48) +
                            48 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * li2((-1 + cost) / (1 + cost)) -
                            6 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * li2((1 + cost) / (-1 + cost)) +
                            18 * (-1 + cost) * (-1 + cost * cost) * log(-2 / (-1 + cost)) + 18 * (1 + cost) * (-1 + cost * cost) * log(2 / (1 + cost)) +
                            9 * (-1 + cost) * (1 + cost) * (-1 + cost * cost) * log(-((-1 + cost) / (1 + cost))) -
                            3 * (-1 + cost) * (1 + cost) * (-3 + 3 * cost + 16 * (sw * sw) - 16 * (sw * sw * sw * sw)) * (log(mz * mz) * log(mz * mz)) -
                            3 * (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) *
                                (log(-((-1 + cost) * (mz * mz)) / 2.) * log(-((-1 + cost) * (mz * mz)) / 2.)) +
                            24 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * (log(((1 + cost) * (mz * mz)) / 2.) * log(((1 + cost) * (mz * mz)) / 2.)) +
                            6 * (1 - cost * cost) * log(mz * mz) * (-((-1 + cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.)) + 8 * (1 + cost) * (sw * sw) * (-1 + sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) -
                            12 * ((-1 + cost) * (1 + cost) * (-3 + 3 * cost + 16 * (sw * sw) - 16 * (sw * sw * sw * sw)) * log(mz * mz) - (-1 + cost) * (-1 + cost * cost) * (3 - 8 * (sw * sw) + 8 * (sw * sw * sw * sw)) * log(-((-1 + cost) * (mz * mz)) / 2.) + 8 * (1 + cost) * (-1 + cost * cost) * (sw * sw) * (-1 + sw * sw) * log(((1 + cost) * (mz * mz)) / 2.)) * log(1 - s / (mzs - I * gamz * mz)))) /
                    (1728. * ((-1 + cost * cost) * (-1 + cost * cost)) * (cw * cw) * (mz * mz) * (Pi * Pi) * (sw * sw)));
        else
            return (0);

        break;
    }
}

// Leading pole contribution of gamma-Z boxes (after IR subtraction)
Cplx bgzR(const int it, const int ot, const inval &input)
{
    const inval *ival;
    ival = &input;
    double cost = 0.5;
    double el = sqrt(4 * Pi * ival->get(al));
    double z = (1 - cost) / (1 + cost);
    return (
        -(el * el / (4 * Pi * Pi)) * Qf[it] * Qf[ot] * (li2(-z) - li2((-1.) / z) + 0.5 * log(z) * (log(4 / (1 - cost * cost)) + 1)));
}

//axuiliary functions

Cplx deX2(const int it, const int ot, const inval &input)

{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)), cost = 0.5,
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           gamz = ival->get(GamZ),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz,
           elq = el * el * el * el,
           swq = sw * sw * sw * sw,
           s = 88. * 88.;
    Cplx s0 = mzs - I * mz * gamz;
    return (-imag(SigZ1p(input)) * imag(SigZ1p(input)) + 2. * bgzR(it, ot, input));
}
Cplx rAAI(const int it, const int ot, const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)), cost = 0.5,
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           gamz = ival->get(GamZ),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz,
           elq = el * el * el * el,
           swq = sw * sw * sw * sw,
           s = 88. * 88.;
    Cplx s0 = mzs - I * mz * gamz;

    return (imag(Z1(it, AXV, input)) / z0(it, AXV, input) + imag(Z1(ot, AXV, input)) / z0(ot, AXV, input) - imag(SigZ1p(input)));
}
Cplx Ixf(const int form, const int ftyp, const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)), cost = 0.5,
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           gamz = ival->get(GamZ),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz,
           elq = el * el * el * el,
           swq = sw * sw * sw * sw;
           
    Cplx s0 = mzs - I * mz * gamz;

    if (form == VEC)
        return ((az0(ftyp, input) * imag(Z1(ftyp, VEC, input)) - vz0(ftyp, input) * imag(Z1(ftyp, AXV, input))) / (az0(ftyp, input) * az0(ftyp, input)));
    else
        return (0.);
}
Cplx xIij(const int it, const int ot, const int iff, const int off, const inval &input)
{
    const inval *ival;
    ival = &input;
    double el = sqrt(4 * Pi * ival->get(al)), cost = 0.5,
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           gamz = ival->get(GamZ),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz,
           elq = el * el * el * el,
           swq = sw * sw * sw * sw;
    Cplx s0 = mzs - I * mz * gamz;

    return (imag(Z1p(it, iff, input)) / z0(it, iff, input) + imag(Z1p(ot, off, input)) / z0(ot, off, input) - 0.5 * imag(SigZ1p2(input)));
}
/*Cplx Bij( const int ot, const int if2, const int of2)
{   Cplx bsum;
    for(int x = 0, x < 3, x++)
    bsum += B1(ot, if2, of2, x, input)
    return( bsum);
}*/
//it is hard to define QWf since it passes an object from SW_SMNLO to the argument, I have not figured out how to do so.
/*Cplx QWf(const int form, const int ftyp, const psobs &SWf, const inval &input) // define &SWf as the obj
{
    const inval *ival;    // define ival are pointers towards the objs of class inval
    const psobs *swef;    
    ival = &input;        // pass the initial address of those objs declared in functions to pointers
    swef = &SWf;
    double el = sqrt(4 * Pi * ival->get(al)), cost = 0.5,
           mz = ival->get(MZ),
           mw = ival->get(MW),
           mt = ival->get(MT),
           mh = ival->get(MH),
           delal = ival->get(Delal),
           gamz = ival->get(GamZ),
           cw = mw / mz,
           sw = sqrt(1 - cw * cw),
           mws = mw * mw,
           mhs = mh * mh,
           mts = mt * mt,
           mzs = mz * mz,
           elq = el * el * el * el,
           swq = sw * sw * sw * sw,
           s = 88. * 88.;
    Cplx s0 = mzs - I * mz * gamz;
    Cplx   sweff = swef->result();
    switch ( form )
    {
    case VEC:
    
    if (ftyp == NEU)
    return(1.);

    else
    return(
        1 - 4. *fabs(Qf[ftyp]) * sweff);
    
    break;
    case AXV:
    return (1.);
    break;
    }
}*/
/*Cplx Uijkl(const int it, const int ot, const int if1, const int of1, const int  if2, const int of2)
{
    return(QWf(if1, it, SWf, input) * QWf(if2, it, SWf, input) * Ixf(of1, ot, input) * Ixf(of2, ot, input)
    + QWf(of1, ot, SWf, input) * QWf(of2, ot, SWf, input) * Ixf(if1, it, input) * Ixf(if2, it, input)
    +(QWf(if1, it, SWf, input) * Ixf(if2, it ,input) - QWf(if2, it, SWf, input) * Ixf(if1, it, input))
    * (QWf(of2, ot, SWf, input) * Ixf(of1, ot, input)-QWf(of1, ot, SWf, input) * Ixf(of2, ot, input)));

}
Cplx Xijkl( const int it, const int ot, const int if1, const int of1, const int if2, const int of2)
{
    return(z0(it, if1, input) * z0(ot, of1, input) * (z0(it, if2, input) * real(Z1p(ot, if2, input)) 
    + z0(ot, of2, input)* real(Z1p(it, if2, input)) - 0.5 * z0(it, if2, input) * z0(ot ,of2, input) real(SigZ1p2(input)) 
    )
    +z0(it, if2, input) * z0(ot, of2, input) * (z0(it, if1, input) * real(Z1p(ot, of1, input)) 
    + z0(ot, of1, input)* real(Z1p(it, if1, input)) - 0.5 * z0(it, if1, input) * z0(ot ,of1, input) real(SigZ1p2(input)) 
    )
    );
}

Cplx Yijkl( const int it, const int ot, const int if1, const int of1, const int if2, const int of2)
{
    return ( z0(it, if1, input) * z0(ot,of1,input) * g0(it, if2, input) * G1(ot, of2, input)
+ z0(it, if1, input) * z0(ot, of1, input) * G1(it, if2, input) * g0(ot, of2, input)
+ z0(it, if1, input) *conj( Z1(ot, of1, input)) * g0(it, if2, input) * g0(ot, of2, input)
+ conj(Z1(it, if2, input)) * z0(ot, of1, input) * g0(it, if2, input) * g0(ot ,of2, input)*(1.-conj(SigZ1p(input))-SigG1(input)/mzs)
+ mzs * z0(it, if1, input) * z0(ot, of1, input) * Bij(ot, if2, of2)
+ z0(it, if2, input) * z0(ot,of2,input) * g0(it, if1, input) * G1(ot, of1, input)
+ z0(it, if2, input) * z0(ot, of2, input) * G1(it, if1, input) * g0(ot, of1, input)
+ z0(it, if2, input) *conj( Z1(ot, of2, input)) * g0(it, if1, input) * g0(ot, of1, input)
+ conj(Z1(it, if1, input)) * z0(ot, of2, input) * g0(it, if1, input) * g0(ot ,of1, input)*(1.-conj(SigZ1p(input))-SigG1(input)/mzs)
+ mzs * z0(it, if2, input) * z0(ot, of2, input) * Bij(ot, if1, of1)
);
}

Cplx Vijkl(const int it, const int ot, const int if1, const int of1, const int if2, const int of2)
{
    return( g0(it, if1, input) * g0(ot, of1, input)*(0.5 *(g0(it, if2, input) * g0(ot, of2, input)-z0(ot, if2, input) * z0(ot, of2, input)))
    + g0(it, if2, input) * g0(ot, of2, input)*(0.5 *(g0(it, if1, input) * g0(ot, of1, input)-z0(ot, if1, input) * z0(ot, of1, input))));
}*/
/*
int main()
{
    inval myinput(StanMod);
    myinput.set(al, 0.303 * 0.303 / (4. * Pi));
    myinput.set(MZ, 91.0);
    myinput.set(GamZ, 2.5);
    myinput.set(MW, 80.0);
    myinput.set(MT, 173.0);
    myinput.set(MH, 125.0);
    myinput.set(Delal, 0.059);
    double mz = myinput.get(MZ);
    double s = 88. * 88.;

    cout << "Bee is" << B1(LEP, VEC, VEC, GG, myinput, s)+ B1(LEP, VEC, VEC, ZZWW, myinput, s)+B1(LEP, VEC, VEC, GZ, myinput, s)<< endl;

    cout << " Ixf difference " << "lep" << Ixf(VEC, LEP, myinput) <<endl;
    cout << " Ixf difference " << "el" << Ixf(VEC, LEP, myinput) <<endl;
}
*/