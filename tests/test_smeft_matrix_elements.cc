/*
  File: test_smeft_matrix_elements.cc
  Purpose: Test standalone contact amplitudes, flavor/generation matching,
           coefficient aliases, normalization and invalid-process rejection.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-16
*/

// Standalone contacts: normalization, flavor matching and existing aliases.
// Full finite amplitudes and contact counting belong to test_native_smeft.
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include "griffin/smeft/xsc_smeft.h"

using namespace griffin;
using smeft::Generation;
using smeft::Input;

namespace {
const Generation one=Generation::first, two=Generation::second, three=Generation::third;
void require(bool condition, const std::string& message)
{
  if(!condition) throw std::runtime_error(message);
}
void close(Cplx value, Cplx expected, const std::string& message)
{
  require(std::isfinite(real(value)) && std::isfinite(imag(value))
      && std::abs(value-expected)<=1e-20+2e-12*std::abs(expected),message);
}
template<class Exception,class Function>
void throws(Function f, const std::string& message)
{
  bool caught=false;
  try { f(); } catch(const Exception&) { caught=true; }
  require(caught,message);
}
Input empty()
{
  Input input; input.setLambdaGeV(1000.0); return input;
}
double projected(double ll,double lr,double rl,double rr,Current ci,Current cf,
                 double inverseLambdaSquared)
{
  // Left/right signs in the amplitude's vector/axial convention.
  const double initialRight=ci==Current::axial?-1.0:1.0;
  const double finalRight=cf==Current::axial?-1.0:1.0;
  return .25*inverseLambdaSquared*(ll+finalRight*lr+initialRight*rl
                                  +initialRight*finalRight*rr);
}
void check(const Input& input,Fermion final,double ll,double lr,double rl,double rr)
{
  for(Current ci:{Current::vector,Current::axial})
    for(Current cf:{Current::vector,Current::axial})
    {
      const double expected=projected(ll,lr,rl,rr,ci,cf,input.inverseLambdaSquaredGeV());
      const Cplx actual=smeft::fourFermionContactAmplitude(Fermion::electron,final,ci,cf,input);
      close(actual,expected,"contact chiral projection/flavor normalization mismatch");
      Input scaled=input; scaled.setLambdaGeV(2000.0);
      close(smeft::fourFermionContactAmplitude(Fermion::electron,final,ci,cf,scaled),expected/4.0,
            "contact inverse-scale normalization changed");
    }
}
void testQuarks()
{
  const Fermion down[]={Fermion::d,Fermion::s,Fermion::b};
  const Generation gen[]={one,two,three};
  for(int i=0;i<3;++i)
  {
    Input input=empty();
    input.setClq1(one,one,gen[i],gen[i],.8);
    input.setClq3(one,one,gen[i],gen[i],-.2);
    input.setCld(one,one,gen[i],gen[i],-.4);
    input.setCqe(gen[i],gen[i],one,one,.3);
    input.setCed(one,one,gen[i],gen[i],-.1);
    check(input,down[i],.6,-.4,.3,-.1);
    for(int j=0;j<3;++j) if(j!=i) check(input,down[j],0,0,0,0);
  }
  const Fermion up[]={Fermion::u,Fermion::c};
  for(int i=0;i<2;++i)
  {
    Input input=empty();
    input.setClq1(one,one,gen[i],gen[i],.5);
    input.setClq3(one,one,gen[i],gen[i],.2);
    input.setClu(one,one,gen[i],gen[i],.1);
    input.setCeq(one,one,gen[i],gen[i],-.4);
    input.setCeu(one,one,gen[i],gen[i],.6);
    close(input.Cqe(gen[i],gen[i],one,one),-.4,"Ceq/Cqe ordering alias lost");
    check(input,up[i],.3,.1,-.4,.6);
    check(input,up[1-i],0,0,0,0);
  }
}
void testChargedLeptons()
{
  const Fermion final[]={Fermion::muon,Fermion::tau};
  const Generation gen[]={two,three};
  for(int i=0;i<2;++i)
  {
    Input input=empty();
    input.setCll(one,one,gen[i],gen[i],.2);
    input.setCll(one,gen[i],gen[i],one,.3);
    input.setCle(one,one,gen[i],gen[i],-.1);
    input.setCle(gen[i],gen[i],one,one,.4);
    input.setCee(one,gen[i],gen[i],one,.2); // Fierz alias of 11ii
    close(input.Cee(one,one,gen[i],gen[i]),.2,"Cee Fierz alias lost");
    check(input,final[i],.5,-.1,.4,.2);
    check(input,final[1-i],0,0,0,0);
    Input canonical=input;
    canonical.setCee(one,one,gen[i],gen[i],.2);
    check(canonical,final[i],.5,-.1,.4,.2); // no duplicate alias contribution
  }
}
void testNeutrinos()
{
  const Fermion final[]={Fermion::nuElectron,Fermion::nuMuon,Fermion::nuTau};
  const Generation gen[]={one,two,three};
  for(int i=0;i<3;++i)
  {
    Input input=empty();
    input.setCll(one,one,gen[i],gen[i],.6);
    input.setCle(gen[i],gen[i],one,one,-.2);
    if(i!=0) input.setCll(one,gen[i],gen[i],one,9.0);
    // Only nu_e has coincident contractions and therefore twice Cll1111.
    // Charged-lepton exchange contractions do not enter nu_mu/nu_tau.
    check(input,final[i],i==0?1.2:.6,0.0,-.2,0.0);
  }
}

void testEquivalentContactLabels()
{
  const Fermion final[]={Fermion::muon,Fermion::tau};
  const Generation gen[]={two,three};
  for(int i=0;i<2;++i)
  {
    const Generation rr[4][4]={{one,one,gen[i],gen[i]},
        {gen[i],gen[i],one,one},{one,gen[i],gen[i],one},
        {gen[i],one,one,gen[i]}};
    Input repeated=empty();
    for(const auto& label:rr)
    {
      Input single=empty();
      single.setCee(label[0],label[1],label[2],label[3],.2);
      check(single,final[i],0.0,0.0,0.0,.2);
      repeated.setCee(label[0],label[1],label[2],label[3],.2);
      check(repeated,final[i],0.0,0.0,0.0,.2);
    }
    repeated.setCee(rr[3][0],rr[3][1],rr[3][2],rr[3][3],-.3);
    check(repeated,final[i],0.0,0.0,0.0,-.3);
    repeated.setCee(rr[1][0],rr[1][1],rr[1][2],rr[1][3],0.0);
    check(repeated,final[i],0.0,0.0,0.0,0.0);

    const Generation ll[2][2][4]={
      {{one,one,gen[i],gen[i]},{gen[i],gen[i],one,one}},
      {{one,gen[i],gen[i],one},{gen[i],one,one,gen[i]}}
    };
    for(int contraction=0;contraction<2;++contraction)
    {
      Input pair=empty();
      for(const auto& label:ll[contraction])
      {
        Input single=empty();
        single.setCll(label[0],label[1],label[2],label[3],.2);
        check(single,final[i],.2,0.0,0.0,0.0);
        pair.setCll(label[0],label[1],label[2],label[3],.2);
        check(pair,final[i],.2,0.0,0.0,0.0);
      }
      const auto& other=ll[1-contraction][1];
      pair.setCll(other[0],other[1],other[2],other[3],.3);
      check(pair,final[i],.5,0.0,0.0,0.0);
    }
  }
}
void testZeroAndValidation()
{
  Input input=empty();
  for(Fermion f:{Fermion::d,Fermion::u,Fermion::s,Fermion::c,Fermion::b,
                 Fermion::muon,Fermion::tau,Fermion::nuElectron,
                 Fermion::nuMuon,Fermion::nuTau}) check(input,f,0,0,0,0);
  throws<std::invalid_argument>([&]{smeft::fourFermionContactAmplitude(Fermion::muon,
      Fermion::c,Current::vector,Current::vector,input);},"non-electron initial state accepted");
  throws<std::domain_error>([&]{smeft::fourFermionContactAmplitude(Fermion::electron,
      Fermion::electron,Current::vector,Current::vector,input);},"Bhabha contacts accepted");
  for(Current bad:{Current::scalar,Current::pseudoscalar})
  {
    throws<std::invalid_argument>([&]{smeft::fourFermionContactAmplitude(Fermion::electron,
        Fermion::tau,bad,Current::vector,input);},"non-vector/axial initial current accepted");
    throws<std::invalid_argument>([&]{smeft::fourFermionContactAmplitude(Fermion::electron,
        Fermion::tau,Current::vector,bad,input);},"non-vector/axial final current accepted");
  }
  for(int invalid : {-1,0,6,10,17,20,999})
  {
    throws<std::invalid_argument>([&]{smeft::fourFermionContactAmplitude(Fermion::electron,
        static_cast<Fermion>(invalid),Current::vector,Current::vector,input);},
        "typed contact API accepted invalid final fermion");
  }
}
} // namespace

int main()
{
  try
  {
    testQuarks(); testChargedLeptons(); testNeutrinos();
    testEquivalentContactLabels(); testZeroAndValidation();
    std::cout << "Standalone contacts: flavors, normalization, aliases, validation PASS\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "Contact test failed: " << error.what() << '\n'; return 1;
  }
}
