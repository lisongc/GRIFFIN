/*
  Exhaustive storage/API regression test for griffin::smeft::Input.

  The test assigns a distinct finite value to every component documented in
  SMEFT_IMPLEMENTATION.md:
    2 scalar values
    7 * 3^2 = 63 two-index values
    10 * 3^4 = 810 four-index values
  for a total of 875 Wilson-coefficient values. It also checks validation,
  default-zero behavior, copying, and clearing.
*/

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "smeft/SMEFTInput.h"

using griffin::smeft::Generation;
using griffin::smeft::Input;

namespace {

typedef void (Input::*MatrixSetter)(Generation, Generation, double);
typedef double (Input::*MatrixGetter)(Generation, Generation) const;
typedef void (Input::*TensorSetter)(Generation, Generation, Generation,
                                    Generation, double);
typedef double (Input::*TensorGetter)(Generation, Generation, Generation,
                                      Generation) const;

// Iterate in the same one-based order used by the public API and JSON formulas.
const Generation generations[3] = {
  Generation::first,
  Generation::second,
  Generation::third
};

int generationNumber(Generation generation)
{
  return static_cast<int>(generation);
}

void require(bool condition, const std::string& message)
{
  if(!condition)
    throw std::runtime_error(message);
}

double nextValue(double& magnitude, std::size_t count)
{
  // Alternating signs and unique magnitudes expose indexing/aliasing mistakes.
  const double value = (count % 2 == 0) ? magnitude : -magnitude;
  magnitude += 0.001;
  return value;
}

void assignMatrix(Input& input, std::ostream& out, const char* name,
                  MatrixSetter setter, MatrixGetter getter,
                  double& magnitude, std::size_t& count)
{
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      const double value = nextValue(magnitude,count);
      (input.*setter)(generations[i],generations[j],value);
      const double stored = (input.*getter)(generations[i],generations[j]);
      require(stored == value, std::string("round-trip failure for ") + name);
      out << name << '[' << generationNumber(generations[i]) << ','
          << generationNumber(generations[j]) << "] = " << stored << '\n';
      ++count;
    }
  }
}

void assignTensor(Input& input, std::ostream& out, const char* name,
                  TensorSetter setter, TensorGetter getter,
                  double& magnitude, std::size_t& count)
{
  for(int p = 0; p < 3; ++p)
  {
    for(int r = 0; r < 3; ++r)
    {
      for(int s = 0; s < 3; ++s)
      {
        for(int t = 0; t < 3; ++t)
        {
          const double value = nextValue(magnitude,count);
          (input.*setter)(generations[p],generations[r],generations[s],
                          generations[t],value);
          const double stored = (input.*getter)(
              generations[p],generations[r],generations[s],generations[t]);
          require(stored == value,
                  std::string("round-trip failure for ") + name);
          out << name << '[' << generationNumber(generations[p]) << ','
              << generationNumber(generations[r]) << ','
              << generationNumber(generations[s]) << ','
              << generationNumber(generations[t]) << "] = " << stored << '\n';
          ++count;
        }
      }
    }
  }
}

template <typename Exception, typename Function>
void requireThrows(Function function, const std::string& message)
{
  bool threw = false;
  try
  {
    function();
  }
  catch(const Exception&)
  {
    threw = true;
  }
  require(threw,message);
}

} // namespace

int main()
{
  try
  {
    Input input;

    require(input.CphiWB() == 0.0, "CphiWB must default to zero");
    require(input.CphiD() == 0.0, "CphiD must default to zero");
    require(input.Cphil3(Generation::first,Generation::first) == 0.0,
            "matrix coefficients must default to zero");
    require(input.Cll(Generation::first,Generation::second,
                      Generation::second,Generation::first) == 0.0,
            "tensor coefficients must default to zero");

    requireThrows<std::logic_error>([&input]() { input.lambdaGeV(); },
                                    "unset Lambda must throw");
    requireThrows<std::invalid_argument>([&input]() { input.setLambdaGeV(0.0); },
                                         "zero Lambda must throw");
    requireThrows<std::invalid_argument>([&input]() { input.setLambdaGeV(-1.0); },
                                         "negative Lambda must throw");
    requireThrows<std::invalid_argument>([&input]() {
      input.setLambdaGeV(std::numeric_limits<double>::infinity());
    }, "infinite Lambda must throw");
    requireThrows<std::invalid_argument>([&input]() {
      input.setCphiWB(std::numeric_limits<double>::quiet_NaN());
    }, "non-finite Wilson coefficient must throw");
    requireThrows<std::out_of_range>([&input]() {
      input.Cphil1(static_cast<Generation>(0),Generation::first);
    }, "invalid generation must throw");

    input.setLambdaGeV(1000.0);

    std::cout << "SMEFT input parameter test\n";
    std::cout << "==========================\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "LambdaGeV = " << input.lambdaGeV() << '\n';
    std::cout << "InverseLambdaSquaredGeV = "
              << input.inverseLambdaSquaredGeV() << '\n';

    double magnitude = 0.001;
    std::size_t count = 0;

    double value = nextValue(magnitude,count);
    input.setCphiWB(value);
    require(input.CphiWB() == value, "CphiWB round-trip failure");
    std::cout << "CphiWB = " << input.CphiWB() << '\n';
    ++count;

    value = nextValue(magnitude,count);
    input.setCphiD(value);
    require(input.CphiD() == value, "CphiD round-trip failure");
    std::cout << "CphiD = " << input.CphiD() << '\n';
    ++count;

    assignMatrix(input,std::cout,"Cphil1",&Input::setCphil1,&Input::Cphil1,
                 magnitude,count);
    assignMatrix(input,std::cout,"Cphil3",&Input::setCphil3,&Input::Cphil3,
                 magnitude,count);
    assignMatrix(input,std::cout,"Cphie",&Input::setCphie,&Input::Cphie,
                 magnitude,count);
    assignMatrix(input,std::cout,"Cphiq1",&Input::setCphiq1,&Input::Cphiq1,
                 magnitude,count);
    assignMatrix(input,std::cout,"Cphiq3",&Input::setCphiq3,&Input::Cphiq3,
                 magnitude,count);
    assignMatrix(input,std::cout,"Cphiu",&Input::setCphiu,&Input::Cphiu,
                 magnitude,count);
    assignMatrix(input,std::cout,"Cphid",&Input::setCphid,&Input::Cphid,
                 magnitude,count);

    assignTensor(input,std::cout,"Cll",&Input::setCll,&Input::Cll,
                 magnitude,count);
    assignTensor(input,std::cout,"Clq1",&Input::setClq1,&Input::Clq1,
                 magnitude,count);
    assignTensor(input,std::cout,"Clq3",&Input::setClq3,&Input::Clq3,
                 magnitude,count);
    assignTensor(input,std::cout,"Clu",&Input::setClu,&Input::Clu,
                 magnitude,count);
    assignTensor(input,std::cout,"Cld",&Input::setCld,&Input::Cld,
                 magnitude,count);
    assignTensor(input,std::cout,"Cqe",&Input::setCqe,&Input::Cqe,
                 magnitude,count);
    assignTensor(input,std::cout,"Ceu",&Input::setCeu,&Input::Ceu,
                 magnitude,count);
    assignTensor(input,std::cout,"Ced",&Input::setCed,&Input::Ced,
                 magnitude,count);
    assignTensor(input,std::cout,"Cle",&Input::setCle,&Input::Cle,
                 magnitude,count);
    assignTensor(input,std::cout,"Cee",&Input::setCee,&Input::Cee,
                 magnitude,count);

    require(count == 875, "unexpected number of Wilson coefficient values");

    // The historical Ceq ordering reverses the two currents and must address
    // the same storage as canonical Warsaw/WCxf Cqe.
    const Generation one = Generation::first;
    const Generation two = Generation::second;
    const Generation three = Generation::third;
    input.setCeq(one,two,three,one,1.2345);
    require(input.Cqe(three,one,one,two) == 1.2345,
            "legacy Ceq setter did not update canonical Cqe storage");
    input.setCqe(two,three,one,one,-2.3456);
    require(input.Ceq(one,one,two,three) == -2.3456,
            "legacy Ceq getter did not read canonical Cqe storage");

    input.setCee(one,two,two,one,3.4567);
    require(input.Cee(one,one,two,two) == 3.4567,
            "Cee Fierz alias did not update its canonical slot");
    input.setCee(one,one,three,three,-4.5678);
    require(input.Cee(one,three,three,one) == -4.5678,
            "Cee Fierz alias did not read its canonical slot");

    Input copy = input;
    require(copy.lambdaGeV() == input.lambdaGeV(), "copy lost Lambda");
    require(copy.CphiWB() == input.CphiWB(), "copy lost scalar coefficient");
    require(copy.Cphil3(Generation::second,Generation::third)
            == input.Cphil3(Generation::second,Generation::third),
            "copy lost matrix coefficient");
    require(copy.Cll(Generation::first,Generation::second,
                     Generation::third,Generation::first)
            == input.Cll(Generation::first,Generation::second,
                         Generation::third,Generation::first),
            "copy lost tensor coefficient");

    copy.clearWilsonCoefficients();
    require(copy.CphiWB() == 0.0, "clear failed for scalar coefficient");
    require(copy.Cphil3(Generation::second,Generation::third) == 0.0,
            "clear failed for matrix coefficient");
    require(copy.Cll(Generation::first,Generation::second,
                     Generation::third,Generation::first) == 0.0,
            "clear failed for tensor coefficient");
    require(input.CphiWB() != 0.0, "clearing copy modified original input");

    std::cout << "DefinedScalarValues = 2\n";
    std::cout << "DefinedTwoIndexValues = 63\n";
    std::cout << "DefinedFourIndexValues = 810\n";
    std::cout << "DefinedWilsonCoefficientValues = " << count << '\n';
    std::cout << "CeqCompatibilityAlias = PASS\n";
    std::cout << "CeeFierzAlias = PASS\n";
    std::cout << "ValidationChecks = PASS\n";
    std::cout << "STATUS = PASS\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "STATUS = FAIL\n";
    std::cerr << error.what() << '\n';
    return 1;
  }
}
