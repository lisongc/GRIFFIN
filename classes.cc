#include "ff.h"

// the base classes compute results at NLO:

Cplx coeffR::result(void) const
{
       inval myinput(StanMod);
       double
           mz = ival->get(MZ),
           el = sqrt(4 * Pi * ival->get(al)),
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
           cost = 0.5,
           s = 88. * 88.;
       Cplx s0 = mzs - I * mz * gamz;
       
       myinput.set(al, 0.303 * 0.303 / (4. * Pi));
       myinput.set(MZ, 91.0);
       myinput.set(GamZ, 2.5);
       myinput.set(MW, 80.0);
       myinput.set(MT, 173.0);
       myinput.set(MH, 125.0);
       myinput.set(Delal, 0.059);

       Cplx front = 4. * I3f[it] * I3f[ot] * sqrt(FAi->result() * FAo->result());
       Cplx qwe = 1 - 4. * fabs(Qf[it]) * SWi->result();
       Cplx qwf = 1 - 4. * fabs(Qf[ot]) * SWo->result();
       Cplx rnX = 1. + I * rAAI(it, ot, myinput) - 0.5 * sqr(rAAI(it, ot, myinput)) + 0.5 * deX2(it, ot, myinput);

       if (iff == VEC && off == VEC)
       {
              if (ot == NEU)
                     return (front * (qwe * 1. * rnX +
                                      (I - rAAI(it, ot, myinput)) * (qwe * Ixf(VEC, ot, myinput) + 1. * Ixf(VEC, it, myinput)) - Ixf(VEC, it, myinput) * Ixf(VEC, ot, myinput)) +
                             mz * gamz * vz0(it, myinput) * vz0(ot, myinput) * xIij(it, ot, iff, off, myinput));
              else
                     return (front * (qwe * qwf * rnX + (I - rAAI(it, ot, myinput)) * (qwe * Ixf(VEC, ot, myinput) + qwf * Ixf(VEC, it, myinput)) -
                                      Ixf(VEC, it, myinput) * Ixf(VEC, ot, myinput)) +
                             mz * gamz * vz0(it, myinput) * vz0(ot, myinput) * xIij(it, ot, iff, off, myinput));
       }
       if (iff == VEC && off == AXV)

              return (front * (qwe * rnX +
                               (I - rAAI(it, ot, myinput)) * Ixf(VEC, it, myinput))) +
                     mz * gamz * vz0(it, myinput) * az0(ot, myinput) * xIij(it, ot, VEC, AXV, myinput);

       if (iff == AXV && off == VEC)
       {
              if (ot == NEU)
                     return (front * (1. * rnX +
                                      (I - rAAI(it, ot, myinput)) * Ixf(VEC, ot, myinput)) +
                             mz * gamz * az0(it, myinput) * vz0(ot, myinput) * xIij(it, ot, AXV, VEC, myinput));
              else
                     return (front * (qwf * rnX +
                                      (I - rAAI(it, ot, myinput)) * Ixf(VEC, ot, myinput)) +
                             mz * gamz * az0(it, myinput) * vz0(ot, myinput) * xIij(it, ot, AXV, VEC, myinput));
       }
       if (iff == AXV && off == AXV)
              return (front * rnX +
                      mz * gamz * az0(it, myinput) * az0(ot, myinput) * xIij(it, ot, iff, off, myinput));
       else
              return 0;
}

Cplx coeffS::result(void) const
{

       inval myinput(StanMod);

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
              s = mzs;
       Cplx s0 = mzs - I * mz * gamz;

       myinput.set(al, 0.303 * 0.303 / (4. * Pi));
       myinput.set(MZ, 91.0);
       myinput.set(GamZ, 2.5);
       myinput.set(MW, 80.0);
       myinput.set(MT, 173.0);
       myinput.set(MH, 125.0);
       myinput.set(Delal, 0.059);

       return ((1. + I * gamz / mz) * g0(it, iff, myinput) * g0(ot, off, myinput) / (mz * mz) // Sij(0), why complex mass here??? cplx mass issue at tree-level...

               + z0(it, iff, myinput) * Z1p(ot, off, myinput) + Z1p(it, iff, myinput) * z0(ot, off, myinput)

               - 0.5 * z0(it, iff, myinput) * z0(ot, off, myinput) * SigZ1p2(myinput) + (g0(it, iff, myinput) * G1(ot, off, myinput)

                                                                                         + G1(it, iff, myinput) * g0(ot, off, myinput)) /
                                                                                            (mz * mz) +
               (-SigG1(myinput)) / (mzs * mzs) * g0(it, iff, myinput) * g0(ot, off, myinput)

               + B1(ot, iff, off, GG, myinput, mzs) + B1(ot, iff, off, ZZWW, myinput, mzs) + B1(ot, iff, off, GZ, myinput, mzs));
}

Cplx coeffSp::result(void) const
{
       inval myinput(StanMod);

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

       myinput.set(al, 0.303 * 0.303 / (4. * Pi));
       myinput.set(MZ, 91.0);
       myinput.set(GamZ, 2.5);
       myinput.set(MW, 80.0);
       myinput.set(MT, 173.0);
       myinput.set(MH, 125.0);
       myinput.set(Delal, 0.059);

       return (-g0(it, iff, myinput) * g0(ot, off, myinput) / (mz * mz * mz * mz));
}

Cplx matelsq::result(void) const
{
       inval myinput(StanMod);

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

       myinput.set(al, 0.303 * 0.303 / (4. * Pi));
       myinput.set(MZ, 91.0);
       myinput.set(GamZ, 2.5);
       myinput.set(MW, 80.0);
       myinput.set(MT, 173.0);
       myinput.set(MH, 125.0);
       myinput.set(Delal, 0.059);
       /*double gie0 = g0(it, if1, myinput), gke0 = g0(it, if2, myinput),
              gjf0 = g0(ot, of1, myinput), glf0 = g0(ot, of2, myinput),
              zie0 = z0(it, if1, myinput), zke0 = z0(it, if2, myinput),
              zjf0 = z0(ot, of1, myinput), zlf0 = z0(ot, of2, myinput);
       double Yijkl = zie0 * zjf0 * gke0 * glf0 + zke0 * zlf0 * gie0 * gjf0,
              Vijkl = gie0 * gjf0 * (gke0 * glf0 / 2 - zke0 * zlf0) + gke0 * glf0 * (gie0 * gjf0 / 2 - zie0 * zjf0);
       return ((FAi->result() * FAo->result() * ((if1 == VEC) ? (1. - 4 * fabs(Qf[it]) * SWi->result()) : 1) * ((if2 == VEC) ? (1. - 4 * fabs(Qf[it]) * SWi->result()) : 1) * ((of1 == VEC) ? (1. - 4 * fabs(Qf[ot]) * SWo->result()) : 1) * ((of2 == VEC) ? (1. - 4 * fabs(Qf[ot]) * SWo->result()) : 1) + gamz * gamz / (mz * mz) * gie0 * gjf0 * gke0 * glf0 + (s / (mz * mz) - 1) * Yijkl + sqr(s / (mz * mz) - 1) * Vijkl) / (sqr(s - mz * mz) + sqr(mz * gamz)));*/
       Cplx Bij; /*= B1(ot, if1, of1, GG, myinput, s) + B1(ot, if1, of1, GZ, myinput, s) + B1(ot, if1, of1, ZZWW, myinput, s);*/
       for (int x = 0; x < 4; x++)
              Bij += B1(ot, if1, of1, x, myinput, mzs);

       Cplx Bkl;
       for (int n = 0; n < 4; n++)
              Bkl += B1(ot, if2, of2, n, myinput, mzs);

       Cplx qit1 = ((if1 == VEC) ? (1 - 4. * fabs(Qf[it]) * SWi->result()) : 1.);
       Cplx qit2 = ((if2 == VEC) ? (1 - 4. * fabs(Qf[it]) * SWi->result()) : 1.);
       Cplx qot1;
       if (ot == NEU)
              qot1 = 1;
       else
              qot1 = ((of1 == VEC) ? (1 - 4. * fabs(Qf[ot]) * SWo->result()) : 1.);
       Cplx qot2;
       if (ot == NEU)
              qot2 = 1;
       else
              qot2 = ((of2 == VEC) ? (1 - 4. * fabs(Qf[ot]) * SWo->result()) : 1.);

       Cplx Uijkl = qit1 * qit2 * Ixf(of1, ot, myinput) * Ixf(of2, ot, myinput) + qot1 * qot2 * Ixf(if1, it, myinput) * Ixf(if2, it, myinput) + (qit1 * Ixf(if2, it, myinput) - qit2 * Ixf(if1, it, myinput)) * (qot2 * Ixf(of1, ot, myinput) - qot1 * Ixf(of2, ot, myinput));

       Cplx Xijkl =
           z0(it, if1, myinput) * z0(ot, of1, myinput) * (z0(it, if2, myinput) * real(Z1p(ot, of2, myinput)) + z0(ot, of2, myinput) * real(Z1p(it, if2, myinput)) - 0.5 * z0(it, if2, myinput) * z0(ot, of2, myinput) * real(SigZ1p2(myinput))) + z0(it, if2, myinput) * z0(ot, of2, myinput) * (z0(it, if1, myinput) * real(Z1p(ot, of1, myinput)) + z0(ot, of1, myinput) * real(Z1p(it, if1, myinput)) - 0.5 * z0(it, if1, myinput) * z0(ot, of1, myinput) * real(SigZ1p2(myinput)));

       Cplx Yijkl = z0(it, if1, myinput) * z0(ot, of1, myinput) * g0(it, if2, myinput) * G1(ot, of2, myinput) + z0(it, if1, myinput) * z0(ot, of1, myinput) * G1(it, if2, myinput) * g0(ot, of2, myinput) + z0(it, if1, myinput) * conj(Z1(ot, of1, myinput)) * g0(it, if2, myinput) * g0(ot, of2, myinput) + conj(Z1(it, if1, myinput)) * z0(ot, of1, myinput) * g0(it, if2, myinput) * g0(ot, of2, myinput) + z0(it, if1, myinput) * z0(ot, of1, myinput) * g0(it, if2, myinput) * g0(ot, of2, myinput) * (1. - conj(SigZ1p(myinput)) - SigG1(myinput) / mzs) + mzs * z0(it, if1, myinput) * z0(ot, of1, myinput) * Bkl

                    + z0(it, if2, myinput) * z0(ot, of2, myinput) * g0(it, if1, myinput) * G1(ot, of1, myinput) + z0(it, if2, myinput) * z0(ot, of2, myinput) * G1(it, if1, myinput) * g0(ot, of1, myinput) + z0(it, if2, myinput) * conj(Z1(ot, of2, myinput)) * g0(it, if1, myinput) * g0(ot, of1, myinput) + conj(Z1(it, if2, myinput)) * z0(ot, of2, myinput) * g0(it, if1, myinput) * g0(ot, of1, myinput) + z0(it, if2, myinput) * z0(ot, of2, myinput) * g0(it, if1, myinput) * g0(ot, of1, myinput) * (1. - conj(SigZ1p(myinput)) - SigG1(myinput) / mzs) + mzs * z0(it, if2, myinput) * z0(ot, of2, myinput) * Bij;

       Cplx Vijkl = g0(it, if1, myinput) * g0(ot, of1, myinput) * (0.5 * (g0(it, if2, myinput) * g0(ot, of2, myinput)) - z0(it, if2, myinput) * z0(ot, of2, myinput)) + g0(it, if2, myinput) * g0(ot, of2, myinput) * (0.5 * (g0(it, if1, myinput) * g0(ot, of1, myinput)) - z0(it, if1, myinput) * z0(ot, of1, myinput));

       if (ot == NEU)
              return ((FAi->result() * FAo->result()) * (qit1 * 1 * qit2 * 1 * (1 + deX2(it, ot, myinput)) + Uijkl) / pow(fabs(s - s0), 2) - ((gamz / mz) * imag(Yijkl) - pow(gamz / mz, 2) * g0(it, if1, myinput) * g0(ot, of1, myinput) * g0(it, if2, myinput) * g0(ot, of2, myinput)) / pow(fabs(s - s0), 2) + ((s - mzs) / pow(fabs(s - s0), 2)) * (Xijkl + real(Yijkl) / mzs) + ((pow((s - mzs), 2) / pow(fabs(s - s0), 2))) * Vijkl * (1 / (mzs * mzs))

              );
       else
              return (
                  (FAi->result() * FAo->result()) * (qit1 * qot1 * qit2 * qot2 * (1 + deX2(it, ot, myinput)) + Uijkl) / pow(fabs(s - s0), 2) - ((gamz / mz) * imag(Yijkl) - pow(gamz / mz, 2) * g0(it, if1, myinput) * g0(ot, of1, myinput) * g0(it, if2, myinput) * g0(ot, of2, myinput)) / pow(fabs(s - s0), 2) + ((s - mzs) / pow(fabs(s - s0), 2)) * (Xijkl + real(Yijkl) / mzs) + ((pow((s - mzs), 2) / pow(fabs(s - s0), 2))) * Vijkl * (1 / (mzs * mzs))

              );
       /*
    return( qit1 );
    */
}

/*******************SM prediction of axial form factors and effective weak-mixing angles **********************************/

Cplx FA_SMLO::result(void) const
{
       inval myinput(StanMod);
       myinput.set(al, 0.303 * 0.303 / (4. * Pi));
       myinput.set(MZ, 91.0);
       myinput.set(GamZ, 2.5);
       myinput.set(MW, 80.0);
       myinput.set(MT, 173.0);
       myinput.set(MH, 125.0);
       myinput.set(Delal, 0.059);

       double el = sqrt(4 * Pi * ival->get(al)),
              mw = ival->get(MW),
              mz = ival->get(MZ);
       return (sqr(az0(ftyp, myinput)));
}

Cplx SW_SMLO::result(void) const
{
       inval myinput(StanMod);
       myinput.set(al, 0.303 * 0.303 / (4. * Pi));
       myinput.set(MZ, 91.0);
       myinput.set(GamZ, 2.5);
       myinput.set(MW, 80.0);
       myinput.set(MT, 173.0);
       myinput.set(MH, 125.0);
       myinput.set(Delal, 0.059);

       double el = sqrt(4 * Pi * ival->get(al)),
              mw = ival->get(MW),
              mz = ival->get(MZ);
       return ((1 - vz0(ftyp, myinput) / az0(ftyp, myinput)) / (4 * fabs(Qf[ftyp])));
}

Cplx FA_SMNLO ::result(void) const
{
       inval myinput(StanMod);
       myinput.set(al, 0.303 * 0.303 / (4. * Pi));
       myinput.set(MZ, 91.0);
       myinput.set(GamZ, 2.5);
       myinput.set(MW, 80.0);
       myinput.set(MT, 173.0);
       myinput.set(MH, 125.0);
       myinput.set(Delal, 0.059);

       double el = sqrt(4 * Pi * ival->get(al)),
              mz = ival->get(MZ),
              mw = ival->get(MW),
              mt = ival->get(MT),
              mh = ival->get(MH),
              gamz = ival->get(GamZ),
              delal = ival->get(Delal),
              cw = mw / mz,
              sw = sqrt(1 - cw * cw),
              mws = mw * mw,
              mhs = mh * mh,
              mts = mt * mt,
              mzs = mz * mz,
              az0sq = real(az0(ftp, myinput) * az0(ftp, myinput));

       return (az0sq * (1 - real(SigZ1p(myinput))) + 2 * az0(ftp, myinput) * real(Z1(ftp, AXV, myinput)) - 0.5 * az0sq * mz * gamz * imag(SigZ1p2(myinput)));
}
Cplx SW_SMNLO ::result(void) const
{
       inval myinput(StanMod);
       myinput.set(al, 0.303 * 0.303 / (4. * Pi));
       myinput.set(MZ, 91.0);
       myinput.set(GamZ, 2.5);
       myinput.set(MW, 80.0);
       myinput.set(MT, 173.0);
       myinput.set(MH, 125.0);
       myinput.set(Delal, 0.059);
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
              mzs = mz * mz,

              az0sq = real(az0(ftp, myinput) * az0(ftp, myinput)),
              rva0 = vz0(ftp, myinput) / az0(ftp, myinput);

       return ((real(Z1(ftp, AXV, myinput)) * vz0(ftp, myinput) - real(Z1(ftp, VEC, myinput)) * az0(ftp, myinput)) / (4. * fabs(Qf[ftp]) * az0sq) + (1. - rva0) / (4. * fabs(Qf[ftp])));
}
/*
int main()
{
  inval myinput(StanMod);
  myinput.set(al, 1 / 137.037);
  myinput.set(MZ, 91.0);
  myinput.set(GamZ, 2.5);
  myinput.set(MW, 80.0);
  myinput.set(MT, 173.0);
  myinput.set(MH, 125.0);
  myinput.set(Delal, 0.059);

  FA_SMNLO nloFAe(UQU, myinput), nloFAb(BQU, myinput), nloFAs(SQU, myinput);
  cout << " the axial form-factor FA^e is " << nloFAe.result() << endl;
  cout << " the axial form-factor FA^b is " << nloFAb.result() << endl;
  cout << "the axial form-factor FA^s is " << nloFAs.result() << endl;
}*/