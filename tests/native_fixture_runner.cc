// Default: archived fixture inputs through the PUBLIC merged API.
// Explicit --reference-* modes are test-only low-level scattering assemblies,
// not installed options on SMEFTCalculator and not derivative projections.
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "griffin/SMEFTCalculator.h"
#include "reference_scattering.h"
#include "griffin/sm/SMvalG.h"

using namespace griffin;
using namespace griffin::prediction;
using smeft::Generation;

namespace {
enum class Evaluation { Public, ReferenceSMNNLO, ReferenceSMaas };

template<class T> void initialize(T& input)
{
  input.set(InputPar::alpha,1.0/137.03599976);
  input.set(InputPar::ZMassExperiment,91.1876);
  input.set(InputPar::WMassExperiment,80.377);
  input.set(InputPar::ZWidthExperiment,2.4952);
  input.set(InputPar::WWidthExperiment,2.085);
  input.set(InputPar::HiggsMass,125.1);
  input.set(InputPar::topMass,172.5);
  input.set(InputPar::bottomMass,2.87);
  input.set(InputPar::charmMass,0.652);
  input.set(InputPar::strangeMass,0.0);
  input.set(InputPar::upMass,0.0);
  input.set(InputPar::downMass,0.0);
  input.set(InputPar::tauMass,1.777);
  input.set(InputPar::muonMass,0.0);
  input.set(InputPar::electronMass,0.0);
  input.set(InputPar::DeltaAlpha,0.059);
  input.set(InputPar::alphaS,0.1179);
  input.set(InputPar::fermiConstant,1.166379e-5);
}

std::vector<std::string> split(const std::string& line)
{
  std::stringstream stream(line);
  std::string field;
  std::vector<std::string> fields;
  while(std::getline(stream,field,',')) fields.push_back(field);
  return fields;
}

void emitArray(const SMEFTCalculator::Amplitudes& values)
{
  std::cout << '[';
  for(unsigned k=0;k<values.size();++k)
  {
    if(k) std::cout << ',';
    std::cout << '[' << real(values[k]) << ',' << imag(values[k]) << ']';
  }
  std::cout << ']';
}

void evaluate(const std::vector<std::string>& row, const inval& input,
              EWInputScheme scheme, Evaluation evaluation)
{
  if(row.size()!=18) throw std::invalid_argument("expected 18 CSV columns");
  if(row[0].find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-.")
      !=std::string::npos) throw std::invalid_argument("unsafe fixture id");
  const Generation one=Generation::first, two=Generation::second;
  smeft::Input eft;
  eft.setLambdaGeV(std::stod(row[4]));
  eft.setCphiD(std::stod(row[5]));
  eft.setCphiWB(std::stod(row[6]));
  eft.setCphil1(one,one,std::stod(row[7]));
  eft.setCphil1(two,two,std::stod(row[8]));
  eft.setCphil3(one,one,std::stod(row[9]));
  eft.setCphil3(two,two,std::stod(row[10]));
  eft.setCphie(one,one,std::stod(row[11]));
  eft.setCphie(two,two,std::stod(row[12]));
  eft.setCll(one,one,two,two,std::stod(row[13]));
  eft.setCll(one,two,two,one,std::stod(row[14]));
  eft.setCee(one,one,two,two,std::stod(row[15]));
  eft.setCle(one,one,two,two,std::stod(row[16]));
  eft.setCle(two,two,one,one,std::stod(row[17]));
  const ElectroweakSetup context(scheme,input);
  const SMEFTCalculator prediction(context,eft);
  const double energy=std::stod(row[2]);
  const ScatteringPoint point(Fermion::electron,Fermion::muon,
      energy*energy,std::stod(row[3]),0.38937966e9);
  const inval& resolved=prediction.resolvedInput();
  SMEFTCalculator::Amplitudes fullAmplitude,contactAmplitude;
  double fullSigma,referenceSigma;
  if(evaluation==Evaluation::Public)
  {
    const auto sigma=prediction.differentialCrossSection(point);
    fullSigma=sigma.total(); referenceSigma=sigma.reference();
    fullAmplitude=prediction.amplitudes(point);
    contactAmplitude=prediction.contactAmplitudes(point);
  }
  else
  {
    const auto content=evaluation==Evaluation::ReferenceSMNNLO
        ? regression_reference::ScatteringContent::ArchivedSMNNLO
        : regression_reference::ScatteringContent::ExtendedSMaas;
    smeft::Input zero=eft;
    zero.clearWilsonCoefficients();
    fullAmplitude=regression_reference::amplitudes(point,resolved,eft,content);
    const auto referenceAmplitude=regression_reference::amplitudes(
        point,prediction.referenceInput(),zero,content);
    fullSigma=regression_reference::crossSection(point,fullAmplitude);
    referenceSigma=regression_reference::crossSection(point,referenceAmplitude);
    contactAmplitude=regression_reference::contacts(point,eft);
  }
  std::cout << "{\"id\":\"" << row[0] << "\",\"sigma_pb\":" << fullSigma
            << ",\"sm_baseline_pb\":" << referenceSigma
            << ",\"delta_sigma_pb\":" << fullSigma-referenceSigma
            << ",\"MW_pole_GeV\":" << resolved.get(InputPar::WMassComplex)
            << ",\"MZ_pole_GeV\":" << resolved.get(InputPar::ZMassComplex)
            << ",\"MW_running_GeV\":" << resolved.get(InputPar::WMassExperiment)
            << ",\"GW_pole_GeV\":" << resolved.get(InputPar::WWidthComplex)
            << ",\"GW_running_GeV\":" << resolved.get(InputPar::WWidthExperiment);
  const PseudoObservable observables[]={PseudoObservable::FA,PseudoObservable::SW,
                                      PseudoObservable::FV};
  const char* names[]={"FA","SW","FV"};
  for(unsigned k=0;k<3;++k)
    for(unsigned f=0;f<2;++f)
    {
      const Fermion flavor=f==0?Fermion::electron:Fermion::muon;
      const Cplx value=evaluation==Evaluation::Public
          ? prediction.pseudoObservable(observables[k],flavor).total()
          : regression_reference::observable(observables[k],flavor,resolved,eft);
      std::cout << ",\"" << names[k] << (f==0?"_e":"_mu") << "\":" << real(value);
    }
  std::cout << ",\"amplitudes\":";
  emitArray(fullAmplitude);
  std::cout << ",\"contact_amplitudes\":";
  emitArray(contactAmplitude);
  std::cout << "}\n";
}
} // namespace

int main(int argc,char** argv)
{
  try {
    Evaluation evaluation=Evaluation::Public;
    if(argc==2 && std::string(argv[1])=="--reference-smnnlo")
      evaluation=Evaluation::ReferenceSMNNLO;
    else if(argc==2 && std::string(argv[1])=="--reference-aas")
      evaluation=Evaluation::ReferenceSMaas;
    else if(argc!=1)
      throw std::invalid_argument("usage: native_fixture_runner [--reference-smnnlo|--reference-aas]");
    SMval alpha; initialize(alpha);
    SMvalGmu gmu; initialize(gmu);
    std::cout << std::setprecision(17);
    std::string line;
    const std::string header="id,scheme,energy_GeV,cos_theta,Lambda_GeV,CphiD,CphiWB,Cphil1_11,Cphil1_22,Cphil3_11,Cphil3_22,Cphie_11,Cphie_22,Cll1122,Cll1221,Cee1122,Cle1122,Cle2211";
    if(!std::getline(std::cin,line) || line!=header)
      throw std::invalid_argument("unexpected fixture CSV header");
    while(std::getline(std::cin,line)) {
      if(line.empty()) continue;
      const auto row=split(line);
      if(row.size()!=18) throw std::invalid_argument("wrong fixture column count");
      if(row[1]=="alpha") evaluate(row,alpha,EWInputScheme::AlphaMWMZ,evaluation);
      else if(row[1]=="gmu") evaluate(row,gmu,EWInputScheme::GmuAlphaMZ,evaluation);
      else throw std::invalid_argument("unknown fixture input scheme");
    }
    return 0;
  } catch(const std::exception& error) {
    std::cerr << "native fixture runner: " << error.what() << '\n';
    return 1;
  }
}
