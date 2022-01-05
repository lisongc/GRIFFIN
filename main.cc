#include <iostream>
using namespace std;

//#include "oneloop.h"
//#include "ff.h"
#include "classes.h"
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
    Cplx s0 = pow(myinput.get(MZ), 2) - I * myinput.get(GamZ) * myinput.get(MZ);
    double s = 88. * 88.;
    double cost = 0.5;

    // compute vertex form factors:
    FA_SMLO FAe(LEP, myinput), FAb(BQU, myinput);
    SW_SMLO SWe(LEP, myinput), SWb(BQU, myinput);

    cout << " leading order axial-vector form factor and sweff " << endl;
    cout << "F_A^e = " << FAe.result() << endl;
    cout << "F_A^b = " << FAb.result() << endl;
    cout << "sineff^e = " << SWe.result() << endl;
    cout << "sineff^b = " << SWb.result() << endl;
    cout << endl;

    double cme = 88.; // center-of-mass energy

    // compute matrix element for ee->bb using SM form factors:
    coeffR R(LEP, LEP, VEC, VEC, FAe, FAb, SWe, SWb, myinput);
    coeffS S(LEP, LEP, VEC, VEC, myinput);
    coeffSp Sp(UQU, BQU, VEC, VEC, myinput);
    
    double mz = myinput.get(MZ), gz = myinput.get(GamZ);
    Cplx sminuss0 = cme * cme - Cplx(mz * mz, -mz * gz);
    
    Cplx matel1 = R.result() / sminuss0 + S.result() + Sp.result() * sminuss0;
    cout << "SM matrix element = " << matel1 << endl;
    cout << " ( \" ) squared = " << sqr(abs(matel1)) << endl;
    cout << " ( \" ) squared, term by term and truncated = "
         << sqr(abs(R.result() / sminuss0)) +
                2 * real(R.result() / sminuss0 * conj(S.result() + Sp.result() * sminuss0)) + sqr(abs(S.result()))
         << endl;
    cout << endl;

    // directly compute squared matrix element for ee->bb using SM form factors:
    matelsq Msq(LEP, LEP, VEC, VEC, VEC, VEC, FAe, FAe, SWe, SWe, s, cost,
                myinput);
    cout << "Direct computation of SM matrix element squared = "
         << Msq.result() << endl;
    cout << endl;

    //compute NLO coefficient Sij
    // ee -> ll
    coeffS Svvl(LEP, LEP, VEC, VEC, myinput), Sval(LEP, LEP, VEC, AXV, myinput),
        Savl(LEP, LEP, AXV, VEC, myinput), Saal(LEP, LEP, AXV, AXV, myinput),
        // ee -> neu
        Svvn(LEP, NEU, VEC, VEC, myinput), Svan(LEP, NEU, VEC, AXV, myinput),
        Savn(LEP, NEU, AXV, VEC, myinput), Saan(LEP, NEU, AXV, AXV, myinput),
        // ee -> uqu
        Svvu(LEP, UQU, VEC, VEC, myinput), Svau(LEP, UQU, VEC, AXV, myinput),
        Savu(LEP, UQU, AXV, VEC, myinput), Saau(LEP, UQU, AXV, AXV, myinput),
        // ee -> dqu
        Svvd(LEP, DQU, VEC, VEC, myinput), Svad(LEP, DQU, VEC, AXV, myinput),
        Savd(LEP, DQU, AXV, VEC, myinput), Saad(LEP, DQU, AXV, AXV, myinput),
        // ee -> bqu
        Svvb(LEP, BQU, VEC, VEC, myinput), Svab(LEP, BQU, VEC, AXV, myinput),
        Savb(LEP, BQU, AXV, VEC, myinput), Saab(LEP, BQU, AXV, AXV, myinput);

    cout << "S for ll is " << Svvl.result() << "  " << Sval.result() << " " << Savl.result() << " " << Saal.result() << endl;
    cout << "S for nn is " << Svvn.result() << "  " << Svan.result() << " " << Savn.result() << " " << Saan.result() << endl;
    cout << "S for uu is " << Svvu.result() << "  " << Svau.result() << " " << Savu.result() << " " << Saau.result() << endl;
    cout << "S for dd is " << Svvd.result() << "  " << Svad.result() << " " << Savd.result() << " " << Saad.result() << endl;
    cout << "S for bb is " << Svvb.result() << "  " << Svab.result() << " " << Savb.result() << " " << Saab.result() << endl;

    //compute Sp at NLO
    // ee -> ll
    coeffSp Spvvl(LEP, LEP, VEC, VEC, myinput), Spval(LEP, LEP, VEC, AXV, myinput),
        Spavl(LEP, LEP, AXV, VEC, myinput), Spaal(LEP, LEP, AXV, AXV, myinput),
        // ee -> neu
        Spvvn(LEP, NEU, VEC, VEC, myinput), Spvan(LEP, NEU, VEC, AXV, myinput),
        Spavn(LEP, NEU, AXV, VEC, myinput), Spaan(LEP, NEU, AXV, AXV, myinput),
        // ee -> uqu
        Spvvu(LEP, UQU, VEC, VEC, myinput), Spvau(LEP, UQU, VEC, AXV, myinput),
        Spavu(LEP, UQU, AXV, VEC, myinput), Spaau(LEP, UQU, AXV, AXV, myinput),
        // ee -> dqu
        Spvvd(LEP, DQU, VEC, VEC, myinput), Spvad(LEP, DQU, VEC, AXV, myinput),
        Spavd(LEP, DQU, AXV, VEC, myinput), Spaad(LEP, DQU, AXV, AXV, myinput),
        // ee -> bqu
        Spvvb(LEP, BQU, VEC, VEC, myinput), Spvab(LEP, BQU, VEC, AXV, myinput),
        Spavb(LEP, BQU, AXV, VEC, myinput), Spaab(LEP, BQU, AXV, AXV, myinput);

    cout << "Sp for ll is " << Spvvl.result() << "  " << Spval.result() << " " << Spavl.result() << " " << Spaal.result() << endl;
    cout << "Sp for nn is " << Spvvn.result() << "  " << Spvan.result() << " " << Spavn.result() << " " << Spaan.result() << endl;
    cout << "Sp for uu is " << Spvvu.result() << "  " << Spvau.result() << " " << Spavu.result() << " " << Spaau.result() << endl;
    cout << "Sp for dd is " << Spvvd.result() << "  " << Spvad.result() << " " << Spavd.result() << " " << Spaad.result() << endl;
    cout << "Sp for bb is " << Spvvb.result() << "  " << Spvab.result() << " " << Spavb.result() << " " << Spaab.result() << endl;
    //compute the FA^f at NLO order

    FA_SMNLO nloFAe(LEP, myinput), nloFAn(NEU, myinput), nloFAu(UQU, myinput),
        nloFAd(DQU, myinput), nloFAb(BQU, myinput), nloFAs(SQU, myinput);
    cout << " the axial form-factor FA^e is " << nloFAe.result() << endl;
    cout << " the axial form-factor FA^nu is " << nloFAn.result() << endl;
    cout << " the axial form-factor FA^u is " << nloFAu.result() << endl;
    cout << " the axial form-factor FA^d is " << nloFAd.result() << endl;
    cout << " the axial form-factor FA^b is " << nloFAb.result() << endl;
    cout << "the axial form-factor FA^s is " << nloFAs.result() << endl;

    //compute the effective weak-mixing sweff at NLO
    SW_SMNLO nloswefflep(LEP, myinput), nlosweffb(BQU, myinput), nlosweffs(SQU, myinput), nlosweffn(NEU, myinput),
        nlosweffu(UQU, myinput), nlosweffd(DQU, myinput), nlosweffc(CQU, myinput);
    cout << "the sweff leptonic is " << nloswefflep.result() << endl;
    cout << "the sweff bottom is " << nlosweffb.result() << endl;
    cout << "the sweff strange is " << nlosweffs.result() << endl;
    cout << "the sweff neutrino is " << nlosweffn.result() << endl;
    cout << "the sweff upquark is " << nlosweffu.result() << endl;
    cout << "the sweff downquark is " << nlosweffd.result() << endl;
    cout << "the sweff charmquark is " << nlosweffc.result() << endl;

    //compute NLO coefficient Ri (partially NNLO (no NNLO psobs objs involved))

    //ee -> ll
    coeffR Rvvl(LEP, LEP, VEC, VEC, nloFAe, nloFAe, nloswefflep, nloswefflep, myinput),
        Rval(LEP, LEP, VEC, AXV, nloFAe, nloFAe, nloswefflep, nloswefflep, myinput),
        Ravl(LEP, LEP, AXV, VEC, nloFAe, nloFAe, nloswefflep, nloswefflep, myinput),
        Raal(LEP, LEP, AXV, AXV, nloFAe, nloFAe, nloswefflep, nloswefflep, myinput),
        //ee -> nn
        Rvvn(LEP, NEU, VEC, VEC, nloFAe, nloFAn, nloswefflep, nlosweffn, myinput),
        Rvan(LEP, NEU, VEC, AXV, nloFAe, nloFAn, nloswefflep, nlosweffn, myinput),
        Ravn(LEP, NEU, AXV, VEC, nloFAe, nloFAn, nloswefflep, nlosweffn, myinput),
        Raan(LEP, NEU, AXV, AXV, nloFAe, nloFAn, nloswefflep, nlosweffn, myinput),
        //ee -> uu
        Rvvu(LEP, UQU, VEC, VEC, nloFAe, nloFAu, nloswefflep, nlosweffu, myinput),
        Rvau(LEP, UQU, VEC, AXV, nloFAe, nloFAu, nloswefflep, nlosweffu, myinput),
        Ravu(LEP, UQU, AXV, VEC, nloFAe, nloFAu, nloswefflep, nlosweffu, myinput),
        Raau(LEP, UQU, AXV, AXV, nloFAe, nloFAu, nloswefflep, nlosweffu, myinput),
        //ee -> dd
        Rvvd(LEP, DQU, VEC, VEC, nloFAe, nloFAd, nloswefflep, nlosweffd, myinput),
        Rvad(LEP, DQU, VEC, AXV, nloFAe, nloFAd, nloswefflep, nlosweffd, myinput),
        Ravd(LEP, DQU, AXV, VEC, nloFAe, nloFAd, nloswefflep, nlosweffd, myinput),
        Raad(LEP, DQU, AXV, AXV, nloFAe, nloFAd, nloswefflep, nlosweffd, myinput),
        //ee -> bb
        Rvvb(LEP, BQU, VEC, VEC, nloFAe, nloFAb, nloswefflep, nlosweffb, myinput),
        Rvab(LEP, BQU, VEC, AXV, nloFAe, nloFAb, nloswefflep, nlosweffb, myinput),
        Ravb(LEP, BQU, AXV, VEC, nloFAe, nloFAb, nloswefflep, nlosweffb, myinput),
        Raab(LEP, BQU, AXV, AXV, nloFAe, nloFAb, nloswefflep, nlosweffb, myinput);

    cout << "R for ll is " << Rvvl.result() << "  " << Rval.result() << " " << Ravl.result() << " " << Raal.result() << endl;
    cout << "R for nn is " << Rvvn.result() << "  " << Rvan.result() << " " << Ravn.result() << " " << Raan.result() << endl;
    cout << "R for uu is " << Rvvu.result() << "  " << Rvau.result() << " " << Ravu.result() << " " << Raau.result() << endl;
    cout << "R for dd is " << Rvvd.result() << "  " << Rvad.result() << " " << Ravd.result() << " " << Raad.result() << endl;
    cout << "R for bb is " << Rvvb.result() << "  " << Rvab.result() << " " << Ravb.result() << " " << Raab.result() << endl;

    // SM amplitude in terms of R, S, and S'

    Cplx amplvv = Rvvl.result() / sminuss0 + Svvl.result() + Spvvl.result() * sminuss0,
         amplva = Rval.result() / sminuss0 + Sval.result() + Spval.result() * sminuss0,
         amplav = Ravl.result() / sminuss0 + Savl.result() + Spavl.result() * sminuss0,
         amplaa = Raal.result() / sminuss0 + Saal.result() + Spaal.result() * sminuss0,
         ampnvv = Rvvn.result() / sminuss0 + Svvn.result() + Spvvn.result() * sminuss0,
         ampnva = Rvan.result() / sminuss0 + Svan.result() + Spvan.result() * sminuss0,
         ampnav = Ravn.result() / sminuss0 + Savn.result() + Spavn.result() * sminuss0,
         ampnaa = Raan.result() / sminuss0 + Saan.result() + Spaan.result() * sminuss0,
         ampuvv = Rvvu.result() / sminuss0 + Svvu.result() + Spvvu.result() * sminuss0,
         ampuva = Rvau.result() / sminuss0 + Svau.result() + Spvau.result() * sminuss0,
         ampuav = Ravu.result() / sminuss0 + Savu.result() + Spavu.result() * sminuss0,
         ampuaa = Raau.result() / sminuss0 + Saau.result() + Spaau.result() * sminuss0,
         ampdvv = Rvvd.result() / sminuss0 + Svvd.result() + Spvvd.result() * sminuss0,
         ampdva = Rvad.result() / sminuss0 + Svad.result() + Spvad.result() * sminuss0,
         ampdav = Ravd.result() / sminuss0 + Savd.result() + Spavd.result() * sminuss0,
         ampdaa = Raad.result() / sminuss0 + Saad.result() + Spaad.result() * sminuss0,
         ampbvv = Rvvb.result() / sminuss0 + Svvb.result() + Spvvb.result() * sminuss0,
         ampbva = Rvab.result() / sminuss0 + Svab.result() + Spvab.result() * sminuss0,
         ampbav = Ravb.result() / sminuss0 + Savb.result() + Spavb.result() * sminuss0,
         ampbaa = Raab.result() / sminuss0 + Saab.result() + Spaab.result() * sminuss0;

    Cplx matesqlvv = sqr(abs(Rvvl.result() / sminuss0)) +
                     2 * real(Rvvl.result() / sminuss0 * conj(Svvl.result() + Spvvl.result() * sminuss0)) + sqr(abs(Svvl.result())),
         matesqlva = sqr(abs(Rval.result() / sminuss0)) +
                     2 * real(Rval.result() / sminuss0 * conj(Sval.result() + Spval.result() * sminuss0)) + sqr(abs(Sval.result())),
         matesqlav = sqr(abs(Ravl.result() / sminuss0)) +
                     2 * real(Ravl.result() / sminuss0 * conj(Savl.result() + Spavl.result() * sminuss0)) + sqr(abs(Savl.result())),
         matesqlaa = sqr(abs(Raal.result() / sminuss0)) +
                     2 * real(Raal.result() / sminuss0 * conj(Saal.result() + Spaal.result() * sminuss0)) + sqr(abs(Saal.result())),
         matesqnvv = sqr(abs(Rvvn.result() / sminuss0)) +
                     2 * real(Rvvn.result() / sminuss0 * conj(Svvn.result() + Spvvn.result() * sminuss0)) + sqr(abs(Svvn.result())),
         matesqnva = sqr(abs(Rvan.result() / sminuss0)) +
                     2 * real(Rvan.result() / sminuss0 * conj(Svan.result() + Spvan.result() * sminuss0)) + sqr(abs(Svan.result())),
         matesqnav = sqr(abs(Ravn.result() / sminuss0)) +
                     2 * real(Ravn.result() / sminuss0 * conj(Savn.result() + Spavn.result() * sminuss0)) + sqr(abs(Savn.result())),
         matesqnaa = sqr(abs(Raan.result() / sminuss0)) +
                     2 * real(Raan.result() / sminuss0 * conj(Saan.result() + Spaan.result() * sminuss0)) + sqr(abs(Saan.result())),
         matesquvv = sqr(abs(Rvvu.result() / sminuss0)) +
                     2 * real(Rvvu.result() / sminuss0 * conj(Svvu.result() + Spvvu.result() * sminuss0)) + sqr(abs(Svvu.result())),
         matesquva = sqr(abs(Rvau.result() / sminuss0)) +
                     2 * real(Rvau.result() / sminuss0 * conj(Svau.result() + Spvau.result() * sminuss0)) + sqr(abs(Svau.result())),
         matesquav = sqr(abs(Ravu.result() / sminuss0)) +
                     2 * real(Ravu.result() / sminuss0 * conj(Savu.result() + Spavu.result() * sminuss0)) + sqr(abs(Savu.result())),
         matesquaa = sqr(abs(Raau.result() / sminuss0)) +
                     2 * real(Raau.result() / sminuss0 * conj(Saau.result() + Spaau.result() * sminuss0)) + sqr(abs(Saau.result())),
         matesqdvv = sqr(abs(Rvvd.result() / sminuss0)) +
                     2 * real(Rvvd.result() / sminuss0 * conj(Svvd.result() + Spvvd.result() * sminuss0)) + sqr(abs(Svvd.result())),
         matesqdva = sqr(abs(Rvad.result() / sminuss0)) +
                     2 * real(Rvad.result() / sminuss0 * conj(Svad.result() + Spvad.result() * sminuss0)) + sqr(abs(Svad.result())),
         matesqdav = sqr(abs(Ravd.result() / sminuss0)) +
                     2 * real(Ravd.result() / sminuss0 * conj(Savd.result() + Spavd.result() * sminuss0)) + sqr(abs(Savd.result())),
         matesqdaa = sqr(abs(Raad.result() / sminuss0)) +
                     2 * real(Raad.result() / sminuss0 * conj(Saad.result() + Spaad.result() * sminuss0)) + sqr(abs(Saad.result())),
         matesqbvv = sqr(abs(Rvvb.result() / sminuss0)) +
                     2 * real(Rvvb.result() / sminuss0 * conj(Svvb.result() + Spvvb.result() * sminuss0)) + sqr(abs(Svvb.result())),
         matesqbva = sqr(abs(Rvab.result() / sminuss0)) +
                     2 * real(Rvab.result() / sminuss0 * conj(Svab.result() + Spvab.result() * sminuss0)) + sqr(abs(Svab.result())),
         matesqbav = sqr(abs(Ravb.result() / sminuss0)) +
                     2 * real(Ravb.result() / sminuss0 * conj(Savb.result() + Spavb.result() * sminuss0)) + sqr(abs(Savb.result())),
         matesqbaa = sqr(abs(Raab.result() / sminuss0)) +
                     2 * real(Raab.result() / sminuss0 * conj(Saab.result() + Spaab.result() * sminuss0)) + sqr(abs(Saab.result()));
    cout << " ============ NLO SM amplitude in terms of R, S, and S================" << endl;
    cout << " ============ e+e- -> l+l-================" << endl;
    cout << "amplvv = " << amplvv << " square of amp directly " << sqr(abs(amplvv))
         << " truncated..." << matesqlvv << endl;
    cout << "amplva = " << amplva << " square of amp directly " << sqr(abs(amplva))
         << " truncated..." << matesqlva << endl;
    cout << "amplav = " << amplav << " square of amp directly " << sqr(abs(amplav))
         << " truncated..." << matesqlav << endl;
    cout << "amplaa = " << amplaa << " square of amp directly " << sqr(abs(amplaa))
         << " truncated..." << matesqlaa << endl;
    cout << "============== e+e- -> neu neu===================" << endl;

    cout << "ampnvv = " << ampnvv << " square of amp directly " << sqr(abs(ampnvv)) << " truncated " << matesqnvv << endl;
    cout << "ampnva = " << ampnva << " square of amp directly " << sqr(abs(ampnva)) << " truncated " << matesqnva << endl;
    cout << "ampnav = " << ampnav << " square of amp directly " << sqr(abs(ampnav)) << " truncated " << matesqnav << endl;
    cout << "ampnaa = " << ampnaa << " square of amp directly " << sqr(abs(ampnaa)) << " truncated " << matesqnaa << endl;
    cout << "============== e+e- -> u u===================" << endl;
    cout << "ampuvv = " << ampuvv << " square of amp directly " << sqr(abs(ampuvv)) << " truncated " << matesquvv << endl;
    cout << "ampuva = " << ampuva << " square of amp directly " << sqr(abs(ampuva)) << " truncated " << matesquva << endl;
    cout << "ampuav = " << ampuav << " square of amp directly " << sqr(abs(ampuav)) << " truncated " << matesquav << endl;
    cout << "ampuaa = " << ampuaa << " square of amp directly " << sqr(abs(ampuaa)) << " truncated " << matesquaa << endl;
    cout << "============== e+e- -> d d===================" << endl;
    cout << "ampdvv = " << ampdvv << " square of amp directly " << sqr(abs(ampdvv)) << " truncated " << matesqdvv << endl;
    cout << "ampdva = " << ampdva << " square of amp directly " << sqr(abs(ampdva)) << " truncated " << matesqdva << endl;
    cout << "ampdav = " << ampdav << " square of amp directly " << sqr(abs(ampdav)) << " truncated " << matesqdav << endl;
    cout << "ampdaa = " << ampdaa << " square of amp directly " << sqr(abs(ampdaa)) << " truncated " << matesqdaa << endl;
    cout << "============== e+e- -> b b===================" << endl;
    cout << "ampbvv = " << ampbvv << " square of amp directly " << sqr(abs(ampbvv)) << " truncated " << matesqbvv << endl;
    cout << "ampbva = " << ampbva << " square of amp directly " << sqr(abs(ampbva)) << " truncated " << matesqbva << endl;
    cout << "ampbav = " << ampbav << " square of amp directly " << sqr(abs(ampbav)) << " truncated " << matesqbav << endl;
    cout << "ampbaa = " << ampbaa << " square of amp directly " << sqr(abs(ampbaa)) << " truncated " << matesqbaa << endl;

    //direct computation of matrix element square at NLO

    //direct computation of matrix element square at NLO
    //ee -> ll
    matelsq Msqelvvvv(LEP, LEP, VEC, VEC, VEC, VEC, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);
    matelsq Msqelvvva(LEP, LEP, VEC, VEC, VEC, AXV, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);
    matelsq Msqelvvav(LEP, LEP, VEC, VEC, AXV, VEC, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);
    matelsq Msqelvvaa(LEP, LEP, VEC, VEC, AXV, AXV, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);
    matelsq Msqelvava(LEP, LEP, VEC, AXV, VEC, AXV, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);
    matelsq Msqelvaav(LEP, LEP, VEC, AXV, AXV, VEC, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);
    matelsq Msqelvaaa(LEP, LEP, VEC, AXV, AXV, AXV, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);
    matelsq Msqelavav(LEP, LEP, AXV, VEC, AXV, VEC, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);
    matelsq Msqelavaa(LEP, LEP, AXV, VEC, AXV, AXV, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);
    matelsq Msqelaaaa(LEP, LEP, AXV, AXV, AXV, AXV, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);

    //ee -> nn

    matelsq Msqenvvvv(LEP, NEU, VEC, VEC, VEC, VEC, nloFAe, nloFAn, nloswefflep, nlosweffn, s, cost, myinput);
    matelsq Msqenvvva(LEP, NEU, VEC, VEC, VEC, AXV, nloFAe, nloFAn, nloswefflep, nlosweffn, s, cost, myinput);
    matelsq Msqenvvav(LEP, NEU, VEC, VEC, AXV, VEC, nloFAe, nloFAn, nloswefflep, nlosweffn, s, cost, myinput);
    matelsq Msqenvvaa(LEP, NEU, VEC, VEC, AXV, AXV, nloFAe, nloFAn, nloswefflep, nlosweffn, s, cost, myinput);
    matelsq Msqenvava(LEP, NEU, VEC, AXV, VEC, AXV, nloFAe, nloFAn, nloswefflep, nlosweffn, s, cost, myinput);
    matelsq Msqenvaaa(LEP, NEU, VEC, AXV, AXV, AXV, nloFAe, nloFAn, nloswefflep, nlosweffn, s, cost, myinput);
    matelsq Msqenavav(LEP, NEU, AXV, VEC, AXV, VEC, nloFAe, nloFAn, nloswefflep, nlosweffn, s, cost, myinput);
    matelsq Msqenavaa(LEP, NEU, AXV, VEC, AXV, AXV, nloFAe, nloFAn, nloswefflep, nlosweffn, s, cost, myinput);
    matelsq Msqenaaaa(LEP, NEU, AXV, AXV, AXV, AXV, nloFAe, nloFAn, nloswefflep, nlosweffn, s, cost, myinput);

    //ee -> uu

    matelsq Msqeuvvvv(LEP, UQU, VEC, VEC, VEC, VEC, nloFAe, nloFAu, nloswefflep, nlosweffu, s, cost, myinput);
    matelsq Msqeuvvva(LEP, UQU, VEC, VEC, VEC, AXV, nloFAe, nloFAu, nloswefflep, nlosweffu, s, cost, myinput);
    matelsq Msqeuvvav(LEP, UQU, VEC, VEC, AXV, VEC, nloFAe, nloFAu, nloswefflep, nlosweffu, s, cost, myinput);
    matelsq Msqeuvvaa(LEP, UQU, VEC, VEC, AXV, AXV, nloFAe, nloFAu, nloswefflep, nlosweffu, s, cost, myinput);
    matelsq Msqeuvava(LEP, UQU, VEC, AXV, VEC, AXV, nloFAe, nloFAu, nloswefflep, nlosweffu, s, cost, myinput);
    matelsq Msqeuvaaa(LEP, UQU, VEC, AXV, AXV, AXV, nloFAe, nloFAu, nloswefflep, nlosweffu, s, cost, myinput);
    matelsq Msqeuavav(LEP, UQU, AXV, VEC, AXV, VEC, nloFAe, nloFAu, nloswefflep, nlosweffu, s, cost, myinput);
    matelsq Msqeuavaa(LEP, UQU, AXV, VEC, AXV, AXV, nloFAe, nloFAu, nloswefflep, nlosweffu, s, cost, myinput);
    matelsq Msqeuaaaa(LEP, UQU, AXV, AXV, AXV, AXV, nloFAe, nloFAu, nloswefflep, nlosweffu, s, cost, myinput);

    //ee -> dd

    matelsq Msqedvvvv(LEP, DQU, VEC, VEC, VEC, VEC, nloFAe, nloFAd, nloswefflep, nlosweffd, s, cost, myinput);
    matelsq Msqedvvva(LEP, DQU, VEC, VEC, VEC, AXV, nloFAe, nloFAd, nloswefflep, nlosweffd, s, cost, myinput);
    matelsq Msqedvvav(LEP, DQU, VEC, VEC, AXV, VEC, nloFAe, nloFAd, nloswefflep, nlosweffd, s, cost, myinput);
    matelsq Msqedvvaa(LEP, DQU, VEC, VEC, AXV, AXV, nloFAe, nloFAd, nloswefflep, nlosweffd, s, cost, myinput);
    matelsq Msqedvava(LEP, DQU, VEC, AXV, VEC, AXV, nloFAe, nloFAd, nloswefflep, nlosweffd, s, cost, myinput);
    matelsq Msqedvaav(LEP, DQU, VEC, AXV, AXV, VEC, nloFAe, nloFAd, nloswefflep, nlosweffd, s, cost, myinput);
    matelsq Msqedvaaa(LEP, DQU, VEC, AXV, AXV, AXV, nloFAe, nloFAd, nloswefflep, nlosweffd, s, cost, myinput);
    matelsq Msqedavav(LEP, DQU, AXV, VEC, AXV, VEC, nloFAe, nloFAd, nloswefflep, nlosweffd, s, cost, myinput);
    matelsq Msqedavaa(LEP, DQU, AXV, VEC, AXV, AXV, nloFAe, nloFAd, nloswefflep, nlosweffd, s, cost, myinput);
    matelsq Msqedaaaa(LEP, DQU, AXV, AXV, AXV, AXV, nloFAe, nloFAd, nloswefflep, nlosweffd, s, cost, myinput);

    //ee -> bb

    matelsq Msqebvvvv(LEP, BQU, VEC, VEC, VEC, VEC, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);
    matelsq Msqebvvva(LEP, BQU, VEC, VEC, VEC, AXV, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);
    matelsq Msqebvvav(LEP, BQU, VEC, VEC, AXV, VEC, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);
    matelsq Msqebvvaa(LEP, BQU, VEC, VEC, AXV, AXV, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);
    matelsq Msqebvava(LEP, BQU, VEC, AXV, VEC, AXV, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);
    matelsq Msqebvaav(LEP, BQU, VEC, AXV, AXV, VEC, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);
    matelsq Msqebvaaa(LEP, BQU, VEC, AXV, AXV, AXV, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);
    matelsq Msqebavav(LEP, BQU, AXV, VEC, AXV, VEC, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);
    matelsq Msqebavaa(LEP, BQU, AXV, VEC, AXV, AXV, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);
    matelsq Msqebaaaa(LEP, BQU, AXV, AXV, AXV, AXV, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);

    //ee -> bb LO
    matelsq Msqebvvvvlo(LEP, BQU, VEC, VEC, VEC, VEC, nloFAe, nloFAb, nloswefflep, nlosweffb, s, cost, myinput);
    //ee -> ll LO
    matelsq Msqelvvvvlo(LEP, LEP, VEC, VEC, VEC, VEC, nloFAe, nloFAe, nloswefflep, nloswefflep, s, cost, myinput);
    cout << "======================================== " << endl;
    cout << " Re { Mlvv Mlvv* } = " << Msqelvvvv.result() << endl;
    cout << " Re { Mlva Mlva* } = " << Msqelvava.result() << endl;
    cout << " Re { Mlav Mlav* } = " << Msqelavav.result() << endl;
    cout << " Re { Mlaa Mlaa* } = " << Msqelaaaa.result() << endl;
    cout << "=======================================" << endl;
    cout << " Re { Mnvv Mnvv* } = " << Msqenvvvv.result() << endl;
    cout << " Re { Mnva Mnva* } = " << Msqenvava.result() << endl;
    cout << " Re { Mnav Mnav* } = " << Msqenavav.result() << endl;
    cout << " Re { Mnaa Mnaa* } = " << Msqenaaaa.result() << endl;
    cout << "=======================================" << endl;
    cout << " Re { Muvv Muvv* } = " << Msqeuvvvv.result() << endl;
    cout << " Re { Muva Muva* } = " << Msqeuvava.result() << endl;
    cout << " Re { Muav Muav* } = " << Msqeuavav.result() << endl;
    cout << " Re { Muaa Muaa* } = " << Msqeuaaaa.result() << endl;
    cout << "=======================================" << endl;
    cout << " Re { Mdvv Mdvv* } = " << Msqedvvvv.result() << endl;
    cout << " Re { Mdva Mdva* } = " << Msqedvava.result() << endl;
    cout << " Re { Mdav Mdav* } = " << Msqedavav.result() << endl;
    cout << " Re { Mdaa Mdaa* } = " << Msqedaaaa.result() << endl;
    cout << "=======================================" << endl;
    cout << " Re { Mbvv Mbvv* } = " << Msqebvvvv.result() << endl;
    cout << " Re { Mbva Mbva* } = " << Msqebvava.result() << endl;
    cout << " Re { Mbav Mbav* } = " << Msqebavav.result() << endl;
    cout << " Re { Mbaa Mbaa* } = " << Msqebaaaa.result() << endl;

    cout << " the leading order matrix element vv square of ee -> bb " << Msqebvvvvlo.result() << endl;
    cout << " the leading order matrix element vv square of ee -> ll " << Msqelvvvvlo.result() << endl;
}
