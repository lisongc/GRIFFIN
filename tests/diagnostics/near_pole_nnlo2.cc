// Diagnostic of the separate upstream mat_SMNNLO2 path. This reports
// near-pole behavior without asserting that a known defect must persist.
// Inputs match upstream testmatel2; the default SMEFTCalculator uses mat_SMaas.
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/xscnnlo2.h"
#include "griffin/sm/SMval.h"
using namespace griffin;

SMval inputs() {
  SMval p;
  p.set(InputPar::alpha,1/137.03599976); p.set(InputPar::ZMassExperiment,91.1876); p.set(InputPar::WMassExperiment,80.377);
  p.set(InputPar::ZWidthExperiment,2.4952); p.set(InputPar::WWidthExperiment,2.085); p.set(InputPar::HiggsMass,125.1);
  p.set(InputPar::topMass,172.5); p.set(InputPar::bottomMass,2.87); p.set(InputPar::tauMass,1.777);
  p.set(InputPar::charmMass,1.27); p.set(InputPar::strangeMass,.342); p.set(InputPar::muonMass,.1057);
  p.set(InputPar::upMass,.246); p.set(InputPar::downMass,.246); p.set(InputPar::electronMass,.000511);
  p.set(InputPar::DeltaAlpha,.059); p.set(InputPar::alphaS,.1179); return p;
}
double xsec(const std::array<Cplx,4>& a,double s,double c,double nc) {
 return nc*s/(32*Pi)*.38937966e6*real((1+c*c)*(a[0]*conj(a[0])+a[1]*conj(a[1])+a[2]*conj(a[2])+a[3]*conj(a[3]))+4*c*(a[0]*conj(a[3])+a[2]*conj(a[1])));
}
int main(int argc, char**) {
 auto p=inputs();
 if(argc>1) for(InputPar mass : {InputPar::electronMass,InputPar::muonMass,InputPar::upMass,InputPar::downMass,InputPar::strangeMass}) p.set(mass,0.0);
 const double mz=p.get(InputPar::ZMassComplex), mz2=mz*mz, c=.5;
 std::cout<<std::setprecision(16);
 std::cout<<"MZpole "<<mz<<" MZpole2 "<<mz2<<" difference="<<(mz2-mz*mz)<<" fma_difference="<<std::fma(-mz,mz,mz2)<<'\n';
 std::cout<<"flavor ds aas_AA_re aas_AA_im nnlo2_AA_re nnlo2_AA_im aas_background_AA_re aas_background_AA_im nnlo2_background_AA_re nnlo2_background_AA_im aas_sigma_nb nnlo2_sigma_nb\n";
 for (Fermion fin : {Fermion::muon,Fermion::d}) {
  FA_SMNNLO fai(Fermion::electron,p), faf(fin,p); SW_SMNNLO swi(Fermion::electron,p), swf(fin,p);
  mat_SMNNLO2 m(Fermion::electron,fin,Current::vector,Current::vector,fai,faf,swi,swf,mz2+100,c,p);
  for (double ds : {-100.,-10.01,-10.,-9.99,-1.,-.1,-.01,0.,.01,.1,1.,9.99,10.,10.01,100.}) {
   const double s=mz2+ds; m.setkinvar(s,c);
   std::array<Cplx,4> a,b;
   const Current ifs[]={Current::vector,Current::axial,Current::vector,Current::axial}, ofs[]={Current::vector,Current::vector,Current::axial,Current::axial};
   for(int k=0;k<4;k++) {m.setform(ifs[k],ofs[k]); a[k]=m.mat_SMaas::result(); b[k]=m.result();}
   const auto ba=m.mat_SMaas::resoffZ(),bb=m.resoffZ();
   std::cout<<index(fin)<<' '<<ds<<' '<<real(a[3])<<' '<<imag(a[3])<<' '<<real(b[3])<<' '<<imag(b[3])<<' '<<real(ba)<<' '<<imag(ba)<<' '<<real(bb)<<' '<<imag(bb)<<' '<<xsec(a,s,c,fin==Fermion::d?3:1)<<' '<<xsec(b,s,c,fin==Fermion::d?3:1)<<'\n';
  }
  std::cout<<"terms: flavor ds form 2f_re 2f_im fb_re fb_im vpol_re vpol_im\n";
  for(double ds : {-1.,-.1,-.01,0.,.01,.1,1.}) {
   m.setkinvar(mz2+ds,c);
   for(Current form : {Current::vector,Current::axial}) {
    m.setform(form,form);
    const auto twof=m.resoffZ2f(), fb=m.resoffZfb(), vp=m.resvpol();
    std::cout<<"terms "<<index(fin)<<' '<<ds<<' '<<index(form)<<' '<<real(twof)<<' '<<imag(twof)<<' '<<real(fb)<<' '<<imag(fb)<<' '<<real(vp)<<' '<<imag(vp)<<'\n';
   }
  }
 }
}
