// Public zero-Wilson SMEFT predictions and resolved inputs for comparison
// with a separately built, frozen upstream v1.1.1 library.
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "griffin/smeft/SMEFTCalculator.h"
#include "griffin/sm/SMvalG.h"

using namespace griffin;
using namespace griffin::prediction;

namespace {
// Preserve the frozen fixture's column order; mass conventions are explicit.
const InputPar rawParameters[] = {
  InputPar::alpha,InputPar::ZMassExperiment,InputPar::WMassExperiment,
  InputPar::ZWidthExperiment,InputPar::WWidthExperiment,InputPar::HiggsMass,
  InputPar::topMass,InputPar::bottomMass,InputPar::charmMass,
  InputPar::strangeMass,InputPar::upMass,InputPar::downMass,InputPar::tauMass,
  InputPar::muonMass,InputPar::electronMass,InputPar::DeltaAlpha,
  InputPar::alphaS,InputPar::fermiConstant
};
const InputPar resolvedParameters[] = {
  InputPar::WMassComplex,InputPar::ZMassComplex,InputPar::HiggsMass,
  InputPar::electronMass,InputPar::muonMass,InputPar::tauMass,InputPar::downMass,
  InputPar::strangeMass,InputPar::bottomMass,InputPar::upMass,InputPar::charmMass,
  InputPar::topMass,InputPar::alpha,InputPar::alphaS,InputPar::DeltaAlpha,
  InputPar::fermiConstant,InputPar::WWidthComplex,InputPar::ZWidthComplex,
  InputPar::WMassExperiment,InputPar::ZMassExperiment,
  InputPar::WWidthExperiment,InputPar::ZWidthExperiment
};

std::vector<std::string> split(const std::string& line)
{
  std::stringstream stream(line);
  std::string value;
  std::vector<std::string> result;
  while(std::getline(stream,value,',')) result.push_back(value);
  return result;
}

void number(double value)
{
  if(!std::isfinite(value)) throw std::runtime_error("nonfinite prediction");
  std::cout << value;
}

void emitInputs(const inval& input)
{
  std::cout << '[';
  for(unsigned k=0;k<sizeof(resolvedParameters)/sizeof(resolvedParameters[0]);++k)
  {
    if(k) std::cout << ',';
    number(input.get(resolvedParameters[k]));
  }
  std::cout << ']';
}

void emitAmplitudes(const SMEFTCalculator::Amplitudes& amplitudes)
{
  std::cout << '[';
  for(unsigned k=0;k<amplitudes.size();++k)
  {
    if(k) std::cout << ',';
    std::cout << '['; number(real(amplitudes[k]));
    std::cout << ','; number(imag(amplitudes[k])); std::cout << ']';
  }
  std::cout << ']';
}

template<class Input> void evaluate(const std::vector<std::string>& row,
                                    EWInputScheme scheme)
{
  Input input;
  for(unsigned k=0;k<sizeof(rawParameters)/sizeof(rawParameters[0]);++k)
    input.set(rawParameters[k],std::stod(row[6+k]));
  smeft::Input zero;
  zero.setLambdaGeV(1000.0);
  const SMEFTCalculator calculator(ElectroweakSetup(scheme,input),zero);
  const double coordinate=std::stod(row[4]);
  double s;
  if(row[3]=="energy") s=coordinate*coordinate;
  else if(row[3]=="delta_s")
  {
    const double mz=calculator.smInputValues().get(InputPar::ZMassComplex);
    s=mz*mz+coordinate;
  }
  else throw std::invalid_argument("unknown kinematic mode");
  const Fermion final=static_cast<Fermion>(std::stoi(row[2]));
  const ScatteringPoint point(Fermion::electron,final,s,std::stod(row[5]),0.38937966e9);
  const auto total=calculator.amplitudes(point);
  const auto sm=calculator.smAmplitudes(point);
  const auto sigma=calculator.differentialCrossSection(point);
  std::cout << "{\"id\":\"" << row[0] << "\",\"s_GeV2\":"; number(s);
  std::cout << ",\"resolved_inputs\":"; emitInputs(calculator.resolvedInput());
  std::cout << ",\"sm_inputs\":"; emitInputs(calculator.smInputValues());
  std::cout << ",\"total_amplitudes\":"; emitAmplitudes(total);
  std::cout << ",\"sm_amplitudes\":"; emitAmplitudes(sm);
  std::cout << ",\"total_sigma_pb\":"; number(sigma.total());
  std::cout << ",\"sm_sigma_pb\":"; number(sigma.smResult());
  std::cout << ",\"model_correction_pb\":"; number(sigma.modelCorrection());
  std::cout << "}\n";
}
} // namespace

int main()
{
  try {
    std::cout << std::setprecision(17);
    std::string line;
    const std::string header="id,scheme,final,mode,coordinate,cos_theta,alpha,MZ_running,MW_running,GammaZ_running,GammaW_running,MH,MT,MB,MC,MS,MU,MD,ML,MM,ME,DeltaAlpha,alphaS,Gmu";
    if(!std::getline(std::cin,line) || line!=header)
      throw std::invalid_argument("unexpected CSV header");
    while(std::getline(std::cin,line))
    {
      if(line.empty()) continue;
      const auto row=split(line);
      if(row.size()!=24) throw std::invalid_argument("expected 24 CSV fields");
      if(row[0].empty() || row[0].find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos)
        throw std::invalid_argument("invalid case id");
      if(row[1]=="alpha") evaluate<SMval>(row,EWInputScheme::AlphaMWMZ);
      else if(row[1]=="gmu") evaluate<SMvalGmu>(row,EWInputScheme::GmuAlphaMZ);
      else throw std::invalid_argument("unknown input scheme");
    }
    return 0;
  } catch(const std::exception& error) {
    std::cerr << "upstream regression runner: " << error.what() << '\n';
    return 1;
  }
}
