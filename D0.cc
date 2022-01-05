/****************************************************************
 *   D0.h, D0.cc
 *   scalar one-loop box integral in terms of dilogs
 *   Lisong Chen Ayres Freitas
 *   Last edition oct 7 2021
 *   
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * **************************************************************/
#include "D0.h"
#include "C0.h"
#include "li.h"

#define Errstream cerr
#define ZERO_LIMIT 1e-14
#define ZEROS_LIMIT 1e-8
#define db0 double(0)

static const Cplx IEPS(0, 1e-15);
static const double zms = 8100.0; // this the reference zero mass defined in
                                  // A. Denner et al. NPB 367(1991) upon wich the results will not depend.

//define sign function (from C0.cc)

double sign0(double var);

//define sign function (from C0.cc)
Cplx Eta(Cplx c1, Cplx c2);

/*
Cplx P(Cplx y0, Cplx y1, Cplx y2, Cplx y3)
{
    return (k01 * y0 * y1 + k02 * y0 * y2 + k03 * y0 * y3 +
            k12 * y1 * y2 + k13 * y1 * y3 + k23 * y2 * y3 + y0 * y0 +
            y1 * y1 + y2 * y2 + y3 * y3);
}
Cplx Q(Cplx y0, Cplx y1, Cplx y2, Cplx y3)
{
    return ((1 / r02 - r02) * y0 + (k12 - k01 * r02) * y1 + y2 - r02 * r02 * y2 + (k23 - k03 * r02) * y3

    );
}

Cplx Qb(Cplx y0, Cplx y1, Cplx y2, Cplx y3)
{

    return ((k01 - k03 * r13) * y0 + (2 - r13 * (1 / r13 + r13)) * y1 + (k12 - k23 * r13) * y2 +
            (1 / r13 - r13) * y3);
}
*/

//define the etatilde function used in D0 scalar integral
Cplx etatil(Cplx c1, Cplx c2)
{
    double im1 = imag(c1), im2 = imag(c2), im12 = imag(c1 * c2),
           re1 = real(c1), re2 = real(c2), re12 = real(c1 * c2);

    if (im2 > ZERO_LIMIT)
        return Eta(c1, c2);
    if (re2 < db0 && im2 < ZERO_LIMIT)
    {
        if (im1 > db0 && im2 > db0)
            return (-TwoPiI);
        if ((im1 > db0 && im2 < db0) || (im1 < db0 && im2 > db0) || (im1 == db0))
            return (Cplx(0));
        if (im1 < db0 && im2 < db0)
            return (TwoPiI);
    }
    if (re2 > 0 && im2 < ZERO_LIMIT)
        return Cplx(0);
}

//the solutions of rij, and rijtilde.
inline Cplx solr(double k)
{
    return ((k + sign0(k) * sqrt(Cplx(k * k - 4.0))) / 2.0);
}

Cplx solir(Cplx k)
{
    return ((k - IEPS + sign0(real(k)) * sqrt((k - IEPS) * (k - IEPS) - Cplx(4., 0.))) / 2.0);
}
//define function in a compact form ( li2 + eta * log)
Cplx del(Cplx a, Cplx b)
{
    return (
        li2(1. + a * b) + Eta(-b, a) * log(1. + a * b));
}

//D0(s1,s2,s3,s4;s12,s23;m0s,0,0,m3s)
Cplx D0(double ps10, double ps20, double ps30, double ps40, double s12,
        double s23, double msq0, double msq1, double msq2, double msq3)
{
    double k01, k02, k03, k12, k13, k23, temp;
    Cplx r01, r02, r03, r12, r13, r23;
    Cplx ik01, ik02, ik03, ik12, ik13, ik23;
    double kij[6], gammajk[3][3], ms[4], ps[6];
    Cplx rij[6], rt[6], x[3], xkj[3][4], s[4], ikij[6];
    Cplx a, b, c, d, res0, res1;
    int n, j, k;

    if ((msq0 == 0. && msq1 == 0.) && (msq2 != 0. && msq3 != 0.)) //D0(s1,s2,s3,s4;s12,s23;m0s,0,0,m3s)
    {
        temp = ps10; // swapping variables according to NPB 361 1991
        ps10 = ps20;
        ps20 = temp;

        temp = ps30;
        ps30 = ps40;
        ps40 = temp;

        temp = msq3;
        msq3 = msq0;
        msq0 = temp;

        temp = msq2;
        msq2 = msq3;
        msq3 = temp;

        k01 = (msq0 - ps10) / sqrt(zms * msq0);
        k02 = (msq0 - s12) / sqrt(zms * msq0);
        k03 = (msq3 + msq0 - ps40) / sqrt(msq3 * msq0);
        k12 = (-ps20) / zms;
        k13 = (msq3 - s23) / sqrt(zms * msq3);
        k23 = (msq3 - ps30) / sqrt(zms * msq3);

        kij[0] = k01;
        kij[1] = k02;
        kij[2] = k03;
        kij[3] = k12;
        kij[4] = k13;
        kij[5] = k23;

        for (int n = 0; n < 6; n++)
        {
            rij[n] = solr(kij[n]);
            // if (imag(rij[n]) == db0)
            // {

            rt[n] = solir(kij[n]);
            // }
            // else
            //     rt[n] = Cplx(0.);
            ikij[n] = Cplx(kij[n]) - IEPS;
        }

        ik01 = ikij[0];
        ik02 = ikij[1];
        ik03 = ikij[2];
        ik12 = ikij[3];
        ik13 = ikij[4];
        ik23 = ikij[5];

        r01 = rij[0];
        r02 = rij[1];
        r03 = rij[2];
        r12 = rij[3];
        r13 = rij[4];
        r23 = rij[5];

        s[0] = rt[0];
        s[1] = rt[3];
        s[2] = rt[5];
        s[3] = rt[2];
        a = k13 * k23 - k12;
        b = k02 * k13 + k01 * k23 - k03 * k12;
        c = k01 * k02 - k12;
        d = k12;

        Cplx disc = sqrt(b * b - 4 * a * (c + d * IEPS));
        x[0] = 0;
        x[1] = (-b - disc) / (2 * a);
        x[2] = (-b + disc) / (2 * a);

        // if (abs(x[1]) < abs(x[2]))
        //     x[2] = (c + d * IEPS) / (a * x[1]);
        // else
        //    x[1] = (c + d * IEPS) / (a * x[2]);
        xkj[0][0] = Cplx(0.); //convention
        xkj[1][0] = x[1] / r13;
        xkj[2][0] = x[2] / r13;

        xkj[0][1] = Cplx(0.);
        xkj[1][1] = x[1] * r02 / r13;
        xkj[2][1] = x[2] * r02 / r13;

        xkj[0][2] = Cplx(0);
        xkj[1][2] = (x[1] * r02);
        xkj[2][2] = (x[2] * r02);

        xkj[0][3] = Cplx(0);
        xkj[1][3] = x[1];
        xkj[2][3] = x[2];

        for (k = 1; k < 3; k++)
        {
            res1 += pow(-1., k) * (li2(1. + s[3] * xkj[k][3]) + Eta(-xkj[k][3], s[3]) * log(1. + s[3] * xkj[k][3]) + li2(1. + xkj[k][3] / s[3]) + Eta(-xkj[k][3], 1 / s[3]) * log(1. + xkj[k][3] / s[3]) - li2(1. + (ik23) / (ik02)*x[k]) - Eta(-x[k], (ik23) / (ik02)) * log(1 + ik23 / ik02 * x[k]) - li2(1. + (ik13) / (ik01)*x[k]) - Eta(-x[k], ik13 / ik01) * log(1. + (ik13 / ik01) * x[k]) + log(-x[k]) * (log(ik01) + log(ik02) - log(ik12)));
        }
        res1 /= zms * sqrt(msq0 * msq3) * a * (x[2] - x[1]);

        for (k = 1; k < 3; k++)
        {
            res0 += pow(-1., k) * (del(s[3], xkj[k][3]) + del(1 / s[3], xkj[k][3]) - del(ik23 / ik02, x[k]) - del(ik13 / ik01, x[k]) + log(-x[k]) * (log(ik01) + log(ik02) - log(ik12)));
        }

        res0 /= zms * sqrt(msq0 * msq3) * a * (x[1] - x[2]);

        return (res0);
    }
    else if (msq0 == 0. && msq1 != 0. && msq2 != 0. && msq3 != 0.)
    {
        // transpositioning arguments to Denner's case. This needs a better implementation later on.
        temp = ps30;
        ps30 = ps10;
        ps10 = temp;

        temp = ps40;
        ps40 = ps20;
        ps20 = temp;

        temp = msq2;
        msq2 = msq0;
        msq0 = temp;

        temp = msq3;
        msq3 = msq1;
        msq1 = temp;

        k01 = (msq0 + msq1 - ps10) / sqrt(msq0 * msq1);
        k02 = (msq0 - s12) / sqrt(zms * msq0);
        k03 = (msq0 + msq3 - ps40) / sqrt(msq0 * msq3);
        k12 = (msq1 - ps20) / sqrt(zms * msq1);
        k13 = (msq1 + msq3 - s23) / sqrt(msq1 * msq3);
        k23 = (msq3 - ps30) / sqrt(zms * msq3);

        kij[0] = k01;
        kij[1] = k02;
        kij[2] = k03;
        kij[3] = k12;
        kij[4] = k13;
        kij[5] = k23;

        for (int n = 0; n < 6; n++)
        {
            rij[n] = solr(kij[n]);
            // if (imag(rij[n]) == db0)
            // {

            rt[n] = solir(kij[n]);
            // }
            // else
            //     rt[n] = Cplx(0.);
            ikij[n] = Cplx(kij[n]) - IEPS;
        }

        ik01 = ikij[0];
        ik02 = ikij[1];
        ik03 = ikij[2];
        ik12 = ikij[3];
        ik13 = ikij[4];
        ik23 = ikij[5];

        r01 = rij[0];
        r02 = rij[1];
        r03 = rij[2];
        r12 = rij[3];
        r13 = rij[4];
        r23 = rij[5];

        s[0] = rt[0];
        s[1] = rt[3];
        s[2] = rt[5];
        s[3] = rt[2];

        a = k23 / r13 - k12;
        b = k02 * (1 / r13 - r13) + k01 * k23 - k03 * k12;
        c = k01 * k02 - k02 * k03 * r13 + r13 * k23 - k12;
        d = k12 - r13 * k23;

        Cplx disc = sqrt(b * b - 4 * a * (c + d * IEPS));
        x[0] = 0;
        x[1] = (-b - disc) / (2 * a);
        x[2] = (-b + disc) / (2 * a);

        xkj[0][0] = Cplx(0.); //convention
        xkj[1][0] = x[1] / r13;
        xkj[2][0] = x[2] / r13;

        xkj[0][1] = Cplx(0.);
        xkj[1][1] = x[1] * r02 / r13;
        xkj[2][1] = x[2] * r02 / r13;

        xkj[0][2] = Cplx(0);
        xkj[1][2] = (x[1] * r02);
        xkj[2][2] = (x[2] * r02);

        xkj[0][3] = Cplx(0);
        xkj[1][3] = x[1];
        xkj[2][3] = x[2];

        for (k = 1; k < 3; k++)
        {
            res0 += pow(-1., k) * (del(s[3], xkj[k][3]) + del(1 / s[3], xkj[k][3]) - del(s[0], xkj[k][0]) - del(1 / s[0], xkj[k][0]) - del(ik23 / ik02, xkj[k][3]) + del(ik12 / ik02, xkj[k][0]) - etatil(-x[k], 1 / rt[4]) * (log(((k01 - k03 * r13) * 1. + (1 / r13 - r13) * x[k])) + log(ik02)));
        }
        res0 /= sqrt(zms) * sqrt(msq0) * sqrt(msq1) * sqrt(msq3) * a * (x[1] - x[2]);

        return res0;
    }
    else
    {
        Errstream << " the rest part is yet to be done" << endl;
        return ((0));
    }
}
