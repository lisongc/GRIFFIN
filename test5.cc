/* testfile for comoputing matrix element far-away from the Z-pole*/

#include <iostream>
using namespace std;

#include "EWPOZ.h"
#include "EWPOZ2.h"
#include "xscnnlo.h"
#include "SMval.h"
#include "li.h"

int main()
{
    SMval myinput; // convert masses from PDG values to complex pole scheme
    myinput.set(al, 1 / 137.03599976);
    myinput.set(als, 0.118);
    myinput.set(MZ, 91.1876);
    myinput.set(MW, 80.4044);
    myinput.set(GamZ, 2.4952);
    myinput.set(GamW, 2.115);
    myinput.set(MH, 125.1);
    myinput.set(MT, 172.5);
    myinput.set(MB, 0);
    myinput.set(Delal, 0.059);
    // process ee -> ll
    int ini = LEP, fin = DQU, iff = VEC, off = VEC;

    double cme0 = 5, // initial center-of-mass energy for the patch far away from the peak (left)
                     //        cme1=80.1876, //..... around the peak

        cost = 0.5, // scattering angle
        pe = 0.,    // polarization degree of the incoming electron beam
        Nc = 1;     // the color factor of the outgoing fermion
    double mz = myinput.get(MZc), gz = myinput.get(GZc), mw = myinput.get(MWc);

    // compute SM form factors at LO and NLO
    FA_SMLO FAe(ini, myinput), FAf(fin, myinput);
    SW_SMLO SWe(ini, myinput), SWf(fin, myinput);
    FA_SMNLO FAe1(ini, myinput), FAf1(fin, myinput);
    SW_SMNLO SWe1(ini, myinput), SWf1(fin, myinput);
    FA_SMNNLO FAe2(ini, myinput), FAf2(fin, myinput);
    SW_SMNNLO SWe2(ini, myinput), SWf2(fin, myinput);

    cout << "mz in cplx-mass scheme = " << mz << endl;
    cout << "mw in cplx-mass scheme = " << mw << endl;
    cout << "F_A^i(NLO) = " << FAe1.result() << endl;
    cout << "F_A^f(NLO) = " << FAf1.result() << endl;
    cout << "sineff^i(NLO) = " << SWe1.result() << endl;
    cout << " sineff^f(NLO) = " << SWf1.result() << endl;

    cout << "=====================================================" << endl;
    Cplx amp0[2][2], amp1[2][2], matenpLO, matepLO, matenpNLO, matepNLO;
    string formlist[2] = {"V", "A"};

    mat_SMNNLO Mij(ini, fin, VEC, VEC, FAe1, FAf1, SWe1, SWf1, cme0 * cme0, cost, myinput);
    matel aij(ini, fin, VEC, VEC, FAe, FAf, SWe, SWf, cme0 * cme0, cost, myinput);

    for (int j = 1; j < 40; j++)
    {
        cme0 += 5;
        Mij.setkinvar(cme0 * cme0, cost);

        Mij.setform(1, 1);

        aij.setkinvar(cme0 * cme0, cost);

        aij.setform(0, 0);
        matenpNLO = Mij.mateloffp1();
        matenpLO = aij.mateloffp1();
        matepNLO = Mij.mateloffp2();
        matepLO = aij.mateloffp2();

        cout << "sqrt(s) = " << cme0 << " "
             << "non-exp NLO Matrix element AV =" << matenpNLO << "  "
             << "exp Matrix element AV = " << matepNLO << "  "
             << " the difference = " << Mij.result2() << endl;
       // cout << "cme is " << cme0 << endl;
        //   cout << " the clog1(1-s/(mzs)) = " << clog1( 1-cme0 * cme0 / (mz * mz)) << endl;
        //  cout << " the clog2(1-s/(mzs)) = " << clog2( 1-cme0 * cme0 / (mz * mz)) << endl;
        //cout << " the clog1((-1+cost)*s/2) = " << clog1((-1 + 0.5) * cme0 * cme0 / 2) << endl;
        //cout << " the clog1(-(-1+cost)*s/2) = " << clog1(-(-1 + 0.5) * cme0 * cme0 / 2) << endl;
        /* cout << "sqrt(s) = " << cme0 << " "
             << "non-exp Matrix element VV =" << matenpLO << "  "
             << "exp Matrix element VV = " << matepLO << "  "
             << " the difference = " << aij.result2() << endl; */
    }

    /*for (int i = 1; i < 4; i++)
    {
        cme0 += 0.5;

        for (int m = VEC; m <= AXV; m++)
        {
            for (int n = VEC; n <= AXV; n++)
            {
                Mij.setkinvar(cme0 * cme0, cost);

                Mij.setform(m, n);

                amp0[m][n] = Mij.result2();

                cout <<"sqrt(s) = "<< cme0 << " "<< formlist[m] << formlist[n] << ":" << amp0[m][n] <<endl;
            }

        }
       // cout << cme0 << " " << 1e+6 * real(((Nc * cme0 * cme0) / (32 * Pi)) * ((1 + cost * cost) * (sqr(abs(amp0[0][0])) + sqr(abs(amp0[1][1])) + sqr(abs(amp0[0][1])) + sqr(abs(amp0[1][0]))) + 4 * cost * (amp0[0][0] * conj(amp0[1][1]) + amp0[0][1] * conj(amp0[1][0])) - 2 * pe * (1 + cost * cost) * (amp0[0][0] * conj(amp0[1][0]) + amp0[0][1] * conj(amp0[1][1])) - 4 * pe * cost * (amp0[0][0] * conj(amp0[0][1]) + amp0[1][0] * conj(amp0[1][1])))) << endl;

    }
*/

    return 0;
}
