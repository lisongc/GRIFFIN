// Diagnostic only: evaluate the unchanged COMPLETE backgrounds, not a model
// for their singular terms. No subtractions or repairs alter the predictions.
#include <griffin/SMEFT.h>
#include <griffin/sm/xscaas.h>
#include <griffin/sm/ff0.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace griffin;

int main()
{
  try {
    SMval input;
    input.set(InputPar::alpha,1.0/137.03599976);
    input.set(InputPar::alphaS,0.1179);
    input.set(InputPar::ZMassExperiment,91.1876);
    input.set(InputPar::ZWidthExperiment,2.4952);
    input.set(InputPar::WMassExperiment,80.377);
    input.set(InputPar::WWidthExperiment,2.085);
    input.set(InputPar::HiggsMass,125.1);
    input.set(InputPar::topMass,172.5);
    input.set(InputPar::bottomMass,2.87);
    input.set(InputPar::charmMass,0.652);
    input.set(InputPar::tauMass,1.777);
    for(InputPar key : {InputPar::electronMass,InputPar::muonMass,
        InputPar::upMass,InputPar::downMass,InputPar::strangeMass})
      input.set(key,0.0);
    input.set(InputPar::DeltaAlpha,0.059);
    input.set(InputPar::fermiConstant,1.166379e-5);
    const double mz=input.get(InputPar::ZMassComplex), Z=mz*mz;
    const double R0=az0(index(Fermion::electron),input)*az0(index(Fermion::muon),input);
    FA_SMNNLO fae(Fermion::electron,input), fam(Fermion::muon,input);
    SW_SMNNLO swe(Fermion::electron,input), swm(Fermion::muon,input);
    std::cout<<std::setprecision(12)
      <<"Known-limitation diagnostic, NOT a physical near-pole prediction\n"
      <<"SM e+e- -> mu+mu-, AA, cos(theta)=0.5; complete unchanged backgrounds\n"
      <<"Pole MZ [GeV]="<<mz<<", pole MW [GeV]="
      <<input.get(InputPar::WMassComplex)<<", R0="<<R0<<'\n'
      <<"h_GeV2 Im_Q_inherited Im_Q_full_background Im_Q_full_amplitude "
      <<"Re_constant_estimate_GeV2 Im_constant_estimate_GeV2\n";
    for(double requested : {0.03,0.01,0.003}) {
      const double sp=Z+requested, sm=Z-requested, h=(sp-sm)/2.0;
      const mat_SMNNLO bp(Fermion::electron,Fermion::muon,Current::axial,
          Current::axial,fae,fam,swe,swm,sp,0.5,input);
      const mat_SMNNLO bm(Fermion::electron,Fermion::muon,Current::axial,
          Current::axial,fae,fam,swe,swm,sm,0.5,input);
      const mat_SMaas fp(Fermion::electron,Fermion::muon,Current::axial,
          Current::axial,fae,fam,swe,swm,sp,0.5,input);
      const mat_SMaas fm(Fermion::electron,Fermion::muon,Current::axial,
          Current::axial,fae,fam,swe,swm,sm,0.5,input);
      const Cplx inherited=-h/(2*R0)*(bp.resoffZ()-bm.resoffZ());
      const Cplx background=-h/(2*R0)*(fp.resoffZ()-fm.resoffZ());
      const Cplx total=-h/(2*R0)*(fp.result()-fm.result());
      const Cplx constant=-h*h/(2*R0)*(fp.resoffZ()+fm.resoffZ());
      for(Cplx value : {inherited,background,total,constant})
        if(!std::isfinite(real(value)) || !std::isfinite(imag(value)))
          throw std::runtime_error("nonfinite diagnostic output");
      std::cout<<h<<' '<<imag(inherited)<<' '<<imag(background)<<' '
        <<imag(total)<<' '<<real(constant)<<' '<<imag(constant)<<'\n';
    }
    std::cout<<"Q(h)=-h*[B(Z+h)-B(Z-h)]/(2*R0); a regular B gives Q -> 0.\n"
      <<"The even combination estimates the constant numerator residual.\n"
      <<"Finite-difference stability diagnoses this implementation, not its intended matching.\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<error.what()<<'\n';
    return 1;
  }
}
