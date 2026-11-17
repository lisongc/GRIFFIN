// Direct Z-vertex corrections at the mass supplied by the caller.
// Finite Gmu input resolution is tested separately in test_native_smeft.
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "griffin/sm/SMvalG.h"
#include "griffin/smeft/ff0_smeft.h"

using namespace griffin;
using smeft::Generation;
using smeft::Input;

namespace {
void require(bool condition, const std::string& message)
{
  if(!condition) throw std::runtime_error(message);
}
void close(double actual, double expected, const std::string& message)
{
  require(std::isfinite(actual) && std::isfinite(expected)
      && std::fabs(actual-expected) <= 1e-15+2e-12*std::fabs(expected),message);
}
SMval makeSM()
{
  SMval input;
  input.set(InputPar::ZMassExperiment,91.1876); input.set(InputPar::WMassExperiment,80.377);
  input.set(InputPar::ZWidthExperiment,2.4952); input.set(InputPar::WWidthExperiment,2.085);
  input.set(InputPar::alpha,1.0/137.03599976); input.set(InputPar::alphaS,0.1179);
  input.set(InputPar::HiggsMass,125.1); input.set(InputPar::topMass,172.5); input.set(InputPar::bottomMass,2.87);
  input.set(InputPar::charmMass,0.652); input.set(InputPar::tauMass,1.777);
  input.set(InputPar::strangeMass,0.0); input.set(InputPar::upMass,0.0); input.set(InputPar::downMass,0.0);
  input.set(InputPar::muonMass,0.0); input.set(InputPar::electronMass,0.0);
  input.set(InputPar::DeltaAlpha,0.059); input.set(InputPar::fermiConstant,1.166379e-5);
  return input;
}
enum Family { charged, neutrino, up, down };
struct Flavor {
  Fermion fermion;
  Generation generation;
  Family family;
  double t3, charge;
};
const Flavor flavors[] = {
  {Fermion::electron,Generation::first,charged,-0.5,-1.0},
  {Fermion::muon,Generation::second,charged,-0.5,-1.0},
  {Fermion::tau,Generation::third,charged,-0.5,-1.0},
  {Fermion::nuElectron,Generation::first,neutrino,0.5,0.0},
  {Fermion::nuMuon,Generation::second,neutrino,0.5,0.0},
  {Fermion::nuTau,Generation::third,neutrino,0.5,0.0},
  {Fermion::u,Generation::first,up,0.5,2.0/3.0},
  {Fermion::c,Generation::second,up,0.5,2.0/3.0},
  {Fermion::d,Generation::first,down,-0.5,-1.0/3.0},
  {Fermion::s,Generation::second,down,-0.5,-1.0/3.0},
  {Fermion::b,Generation::third,down,-0.5,-1.0/3.0}
};
Input empty(double lambda=1000.0)
{
  Input input; input.setLambdaGeV(lambda); return input;
}
void setCurrents(Input& input, const Flavor& f, double singlet,
                 double triplet, double right)
{
  const Generation g=f.generation;
  if(f.family==charged || f.family==neutrino)
  {
    input.setCphil1(g,g,singlet); input.setCphil3(g,g,triplet);
    input.setCphie(g,g,right);
  }
  else
  {
    input.setCphiq1(g,g,singlet); input.setCphiq3(g,g,triplet);
    if(f.family==up) input.setCphiu(g,g,right);
    else input.setCphid(g,g,right);
  }
}
void check(const Flavor& f, const inval& sm, const Input& eft,
           double axial, double vector, const std::string& label)
{
  close(smeft::directZAxialShift(f.fermion,sm,eft),axial,label+" axial");
  close(smeft::directZVectorShift(f.fermion,sm,eft),vector,label+" vector");
}
void testDirectFormulas(const inval& sm)
{
  const double mw=sm.get(InputPar::WMassComplex);
  const double mz=sm.get(InputPar::ZMassComplex);
  const double e=std::sqrt(4.0*Pi*sm.get(InputPar::alpha));
  const double cw=mw/mz, sw=std::sqrt(1.0-cw*cw);
  const double unit=-mw*mw*sw/(e*cw*1000.0*1000.0);
  for(const Flavor& f:flavors)
  {
    check(f,sm,empty(),0.0,0.0,"zero coefficients");
    // Independent single-operator left/right relations for every supported
    // family and generation; no right-handed neutrino is manufactured.
    Input left=empty(); setCurrents(left,f,.017,0.0,0.0);
    check(f,sm,left,unit*.017,unit*.017,"left singlet");
    Input triplet=empty(); setCurrents(triplet,f,0.0,-.013,0.0);
    const double tripletValue=unit*(-2.0*f.t3)*(-.013);
    check(f,sm,triplet,tripletValue,tripletValue,"left triplet");
    Input right=empty(); setCurrents(right,f,0.0,0.0,.011);
    const double rightValue=f.family==neutrino?0.0:unit*.011;
    check(f,sm,right,-rightValue,rightValue,"right current");

    // Direct coupling formulas at this mass, not a mass response.
    Input wb=empty(); wb.setCphiWB(.012);
    const double wbValue=2.0*f.t3*.012*mw*mw/(e*1.0e6);
    check(f,sm,wb,wbValue,wbValue,"universal phiWB");
    Input d=empty(); d.setCphiD(-.018);
    const double dAxial=f.t3*(-.018)*mw*mw*(2.0*cw*cw-1.0)
        /(2.0*e*sw*cw*1.0e6);
    const double dVector=dAxial+(-.018)*mw*mw*f.charge*sw/(e*cw*1.0e6);
    check(f,sm,d,dAxial,dVector,"universal phiD");

    Input mixed=empty(); mixed.setCphiWB(.012); mixed.setCphiD(-.018);
    setCurrents(mixed,f,.017,-.013,.011);
    const double a=wbValue+dAxial+unit*.017+tripletValue-rightValue;
    const double v=wbValue+dVector+unit*.017+tripletValue+rightValue;
    check(f,sm,mixed,a,v,"sum of direct corrections");
    Input opposite=empty(); opposite.setCphiWB(-.012); opposite.setCphiD(.018);
    setCurrents(opposite,f,-.017,.013,-.011);
    check(f,sm,opposite,-a,-v,"direct sign reversal");
    Input scale=mixed; scale.setLambdaGeV(2000.0);
    check(f,sm,scale,a/4.0,v/4.0,"inverse-scale normalization");
  }
}
void testIsolation(const inval& sm)
{
  const double mw=sm.get(InputPar::WMassComplex);
  const double mz=sm.get(InputPar::ZMassComplex);
  const double gw=sm.get(InputPar::WWidthComplex);
  Input contacts=empty();
  contacts.setCll(Generation::first,Generation::second,
                  Generation::second,Generation::first,.15);
  contacts.setCll(Generation::second,Generation::first,
                  Generation::first,Generation::second,-.07);
  for(const Flavor& f:flavors)
    check(f,sm,contacts,0.0,0.0,"contact-only input has no direct Z correction");
  Input electron=empty();
  electron.setCphil1(Generation::first,Generation::first,.1);
  require(smeft::directZAxialShift(Fermion::electron,sm,electron)!=0.0,
          "electron current isolation test is degenerate");
  for(const Flavor& f:flavors)
    if(f.fermion!=Fermion::electron && f.fermion!=Fermion::nuElectron)
      check(f,sm,electron,0.0,0.0,"unrelated flavor changed");
  close(sm.get(InputPar::WMassComplex),mw,"direct function changed supplied MW");
  close(sm.get(InputPar::ZMassComplex),mz,"direct function changed supplied MZ");
  close(sm.get(InputPar::WWidthComplex),gw,"direct function changed supplied width");
}
} // namespace

int main()
{
  try
  {
    const SMval direct=makeSM();
    const SMvalGmu derived(direct);
    // The same fixed-mass formulas hold regardless of how the caller obtained
    // that mass; the direct functions have no input-scheme selector.
    testDirectFormulas(direct); testIsolation(direct);
    testDirectFormulas(derived); testIsolation(derived);
    inval changed(direct);
    changed.set(InputPar::WMassComplex,79.9);
    testDirectFormulas(changed); testIsolation(changed);
    std::cout << "Direct fixed-mass Z shifts: formulas, flavors, scaling, isolation PASS\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "Direct Z-shift test failed: " << error.what() << '\n';
    return 1;
  }
}
