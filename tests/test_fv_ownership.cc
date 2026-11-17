// FV ownership and value-semantics checks. Frozen values were obtained from
// the pre-repair headers/library on 2026-11-16; they test numerical preservation,
// not an independent derivation of the SM form factors.
#include "griffin/sm/EWPOZ2.h"
#include "griffin/sm/SMval.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

// Count only explicitly measured construction/destruction intervals. Input
// setup, numerical evaluation and test reporting run with counting disabled.
namespace allocation_check {
bool enabled = false;
std::size_t created = 0, released = 0;
int allocationsBeforeFailure = -1;
void start() { created = released = 0; enabled = true; }
void stop() { enabled = false; allocationsBeforeFailure = -1; }
}

void* operator new(std::size_t size)
{
  if(allocation_check::enabled && allocation_check::allocationsBeforeFailure == 0)
    throw std::bad_alloc();
  void* pointer = std::malloc(size ? size : 1);
  if(!pointer) throw std::bad_alloc();
  if(allocation_check::enabled)
  {
    ++allocation_check::created;
    if(allocation_check::allocationsBeforeFailure > 0)
      --allocation_check::allocationsBeforeFailure;
  }
  return pointer;
}

void operator delete(void* pointer) noexcept
{
  if(pointer && allocation_check::enabled) ++allocation_check::released;
  std::free(pointer);
}

#if defined(__cpp_sized_deallocation)
void operator delete(void* pointer, std::size_t) noexcept
{
  ::operator delete(pointer);
}
#endif

using namespace griffin;

namespace {

static_assert(std::has_virtual_destructor<psobs>::value,
              "Owned form factors need polymorphic destruction");

void require(bool condition, const char* message)
{
  if(!condition) throw std::runtime_error(message);
}

void close(Cplx actual, Cplx expected, const char* message)
{
  require(std::isfinite(real(actual)) && std::isfinite(imag(actual))
      && std::abs(actual-expected) <= 2e-16+2e-11*std::abs(expected), message);
}

SMval makeInput()
{
  SMval input;
  input.set(InputPar::ZMassExperiment,91.1876);
  input.set(InputPar::WMassExperiment,80.377);
  input.set(InputPar::ZWidthExperiment,2.4952);
  input.set(InputPar::WWidthExperiment,2.085);
  input.set(InputPar::alpha,1.0/137.03599976);
  input.set(InputPar::alphaS,0.1179);
  input.set(InputPar::HiggsMass,125.1);
  input.set(InputPar::topMass,172.5);
  input.set(InputPar::bottomMass,2.87);
  input.set(InputPar::charmMass,1.27);
  input.set(InputPar::tauMass,1.777);
  input.set(InputPar::strangeMass,0.342);
  input.set(InputPar::upMass,0.246);
  input.set(InputPar::downMass,0.246);
  input.set(InputPar::muonMass,0.1057);
  input.set(InputPar::electronMass,0.000511);
  input.set(InputPar::DeltaAlpha,0.059);
  return input;
}

struct Reference {
  Fermion fermion;
  double value[3], error[3]; // LO, NLO, NNLO
};

const Reference references[] = {
  {Fermion::electron,
   {0.00038602189009942354,0.00017673208841195633,0.00020083399294513571},
   {0.0,3.9707307755678049e-05,5.0000000000000004e-06}},
  {Fermion::u,
   {0.0054350523662073684,0.004990283204827038,0.0050855376653853191},
   {0.0,0.00010375895962669773,5.0000000000000004e-06}},
  {Fermion::b,
   {0.016332373659482518,0.016073758822839689,0.016254516792663341},
   {0.0,0.00012119190517150161,1.1e-05}}
};

template<class FV>
FV copyFromTemporaryOwner(const inval& input)
{
  const FV original(Fermion::electron,input);
  const FV copy(original);
  return copy;
}

template<class FV>
void checkValuesAndCopies(const SMval& input, unsigned order)
{
  for(const Reference& point : references)
  {
    const FV value(point.fermion,input);
    close(value.result(),point.value[order],"FV result changed after ownership repair");
    close(value.errest(),point.error[order],"FV error estimate changed after ownership repair");
  }

  SMval alternate = input;
  alternate.set(InputPar::WMassExperiment,80.1);
  FV original(Fermion::electron,input);
  const Cplx originalValue = original.result(), originalError = original.errest();
  FV copy(original);
  close(copy.result(),originalValue,"Copy lost the form factors' perturbative order");
  close(copy.errest(),originalError,"Copy lost the form factors' error estimate");
  copy.setinput(alternate);
  copy.setftype(Fermion::u);
  close(copy.result(),FV(Fermion::u,alternate).result(),"Copied FV setters did not update its children");
  close(original.result(),originalValue,"Changing a copy changed the original FV");

  FV assigned(Fermion::b,alternate);
  assigned = original;
  close(assigned.result(),originalValue,"Copy assignment lost the form factors' order");
  close(assigned.errest(),originalError,"Copy assignment lost the error estimate");
  assigned.setinput(alternate);
  assigned.setftype(Fermion::b);
  close(assigned.result(),FV(Fermion::b,alternate).result(),"Assigned FV setters failed");
  close(original.result(),originalValue,"Changing an assigned copy changed its source");
  FV& alias = original;
  original = alias;
  close(original.result(),originalValue,"Self-copy assignment changed FV");
  close(copyFromTemporaryOwner<FV>(input).result(),originalValue,
        "A copied FV depended on its destroyed source");

  FV moved(std::move(copy));
  close(moved.result(),FV(Fermion::u,alternate).result(),"Move construction lost the owned form factors");
  copy = original; // A moved-from object remains reassignable.
  close(copy.result(),originalValue,"Copy assignment could not restore a moved-from FV");
  allocation_check::start();
  assigned = std::move(moved);
  allocation_check::stop();
  require(allocation_check::created == 0 && allocation_check::released == 2,
          "Move assignment must release the target's old pair and transfer ownership");
  close(assigned.result(),FV(Fermion::u,alternate).result(),"Move assignment changed FV");
  moved = FV(Fermion::electron,input);
  close(moved.result(),originalValue,"Move assignment could not restore a moved-from FV");

  SMval borrowedInput = input;
  {
    FV borrowing(Fermion::electron,borrowedInput);
    borrowedInput.set(InputPar::WMassExperiment,79.9);
    close(borrowing.result(),FV(Fermion::electron,borrowedInput).result(),
          "FV no longer observes changes to its borrowed input");
    require(std::abs(borrowing.result()-originalValue)>1e-12,
            "Borrowed-input check must change the numerical result");
  }
  require(borrowedInput.get(InputPar::WMassComplex)>0.0,
          "Destroying FV damaged its borrowed input");
}

template<class FV>
void checkAllocations(const SMval& input)
{
  for(unsigned repeat=0; repeat<32; ++repeat)
  {
    allocation_check::start();
    {
      FV value(Fermion::electron,input);
      allocation_check::stop();
      require(std::isfinite(real(value.result())),"FV construction produced a nonfinite result");
      allocation_check::enabled = true;
    }
    allocation_check::stop();
    require(allocation_check::created == 2 && allocation_check::released == 2,
            "FV must create only its selected FA/SW pair and release both");
  }

  bool threw = false;
  allocation_check::start();
  allocation_check::allocationsBeforeFailure = 1;
  try { FV value(Fermion::electron,input); }
  catch(const std::bad_alloc&) { threw = true; }
  allocation_check::stop();
  require(threw,"FV construction did not exercise the second-allocation failure");
  require(allocation_check::created == 1 && allocation_check::released == 1,
          "Failure to allocate the second form factor leaked the first");

  psobs* erased = new FV(Fermion::electron,input);
  require(std::isfinite(real(erased->result())),"Polymorphic FV result is nonfinite");
  allocation_check::start();
  delete erased;
  allocation_check::stop();
  require(allocation_check::created == 0 && allocation_check::released == 3,
          "Deleting FV through psobs must release FV and its two children");
}

struct Lifetimes {
  int axial = 0, weak = 0;
  bool failWeakClone = false;
};

class CountedFA : public FA_SMLO {
  Lifetimes& counts_;
protected:
  std::unique_ptr<FA_SMLO> clone() const override
  { return std::unique_ptr<FA_SMLO>(new CountedFA(*this)); }
public:
  CountedFA(Fermion type, const inval& input, Lifetimes& counts)
    : FA_SMLO(type,input), counts_(counts) { ++counts_.axial; }
  CountedFA(const CountedFA& other)
    : FA_SMLO(other), counts_(other.counts_) { ++counts_.axial; }
  ~CountedFA() override { --counts_.axial; }
};

class CountedSW : public SW_SMLO {
  Lifetimes& counts_;
protected:
  std::unique_ptr<SW_SMLO> clone() const override
  {
    if(counts_.failWeakClone) throw std::runtime_error("requested clone failure");
    return std::unique_ptr<SW_SMLO>(new CountedSW(*this));
  }
public:
  CountedSW(Fermion type, const inval& input, Lifetimes& counts)
    : SW_SMLO(type,input), counts_(counts) { ++counts_.weak; }
  CountedSW(const CountedSW& other)
    : SW_SMLO(other), counts_(other.counts_) { ++counts_.weak; }
  ~CountedSW() override { --counts_.weak; }
};

class CountedFV : public FV_SMLO {
  // Finish ownership within each call before evaluating another constructor
  // argument, including under C++11's unspecified argument evaluation order.
  static std::unique_ptr<FA_SMLO> makeAxial(Fermion type, const inval& input,
                                           Lifetimes& counts)
  { return std::unique_ptr<FA_SMLO>(new CountedFA(type,input,counts)); }
  static std::unique_ptr<SW_SMLO> makeWeak(Fermion type, const inval& input,
                                         Lifetimes& counts)
  { return std::unique_ptr<SW_SMLO>(new CountedSW(type,input,counts)); }
public:
  CountedFV(Fermion type, const inval& input, Lifetimes& counts)
    : FV_SMLO(type,input,
        makeAxial(type,input,counts), makeWeak(type,input,counts)) {}
};

void checkExceptionSafetyAndDestruction(const SMval& input)
{
  Lifetimes sourceCounts, targetCounts;
  {
    CountedFV source(Fermion::electron,input,sourceCounts), target(Fermion::u,input,targetCounts);
    const Cplx before = target.result();
    sourceCounts.failWeakClone = true;
    bool threw = false;
    try { target = source; }
    catch(const std::runtime_error&) { threw = true; }
    require(threw,"Assignment did not attempt to clone the actual SW type");
    close(target.result(),before,"A failed clone partially changed the target");
    require(sourceCounts.axial == 1 && sourceCounts.weak == 1
        && targetCounts.axial == 1 && targetCounts.weak == 1,
        "Failed assignment leaked a clone or destroyed the target's children");
    sourceCounts.failWeakClone = false;
    target = source;
    require(sourceCounts.axial == 2 && sourceCounts.weak == 2
        && targetCounts.axial == 0 && targetCounts.weak == 0,
        "Successful assignment did not replace the target's owned children");
  }
  require(sourceCounts.axial == 0 && sourceCounts.weak == 0,
          "Derived FA/SW destructors were not called through their base pointers");
}

} // namespace

int main()
{
  try
  {
    const SMval input = makeInput();
    checkValuesAndCopies<FV_SMLO>(input,0);
    checkValuesAndCopies<FV_SMNLO>(input,1);
    checkValuesAndCopies<FV_SMNNLO>(input,2);
    checkAllocations<FV_SMLO>(input);
    checkAllocations<FV_SMNLO>(input);
    checkAllocations<FV_SMNNLO>(input);
    checkExceptionSafetyAndDestruction(input);
    std::cout << "FV ownership, copy/move and numerical regression checks passed\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    allocation_check::stop();
    std::cerr << error.what() << '\n';
    return 1;
  }
}
