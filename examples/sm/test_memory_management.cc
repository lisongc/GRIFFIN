/*
  File: examples/sm/test_memory_management.cc
  Purpose: Show, with printed messages and runtime checks, how FV_SMLO owns
           and releases its FA/SW helpers during copies, moves, replacement,
           repeated calculations and an interrupted calculation.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-17

  Start reading at main(). The small classes above it only make the real
  FV_SMLO ownership operations visible; they do not change the SM formulas.
  Ordinary users do NOT need these tracing classes in their calculations.

  Build/run from the project root:
    cmake -S . -B build/release -DGRIFFIN_BUILD_EXAMPLES=ON
    cmake --build build/release --target test_memory_management -j 4
    ./build/release/examples/test_memory_management

  A memory leak means an allocated object was not released when it was no
  longer needed. A raw pointer does not automatically release its object:

    // UNSAFE illustration only -- deliberately NOT executed here:
    // FA_SMLO* fa = new FA_SMLO(Fermion::electron, sm);
    // fa = new FA_SMLO(Fermion::electron, sm); // first object is lost!
    // Leaving the scope would not delete the second object either.

  Today's FV_SMLO instead owns its helpers with std::unique_ptr. Replacing
  an owner releases its previous object; destroying it releases its current
  object. No manual delete is needed in the calculation below.

  The counters measure only the instrumented FA/SW objects, not all memory
  allocations in the process. No intentional leak is run. The more detailed
  tests/sm/test_fv_ownership.cc also checks allocation and clone failures.
*/

#include <griffin/sm/SMval.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

using namespace griffin;

namespace {

// Teaching support: keep a visible record of helper lifetimes.
struct HelperCounts {
  int created = 0;
  int destroyed = 0;

  void create(const char* name)
  {
    ++created;
    std::cout << "    create " << name << '\n';
  }
  void destroy(const char* name)
  {
    ++destroyed;
    std::cout << "    destroy " << name << '\n';
  }
  void check(int expectedCreated, int expectedDestroyed, const char* stage) const
  {
    std::cout << "  " << stage << ": created=" << created
              << ", destroyed=" << destroyed
              << ", still alive=" << created - destroyed << '\n';
    // Unlike assert(), this check remains active in a Release build.
    if(created != expectedCreated || destroyed != expectedDestroyed)
      throw std::runtime_error("Unexpected helper lifetime count");
  }
};

// These are the ordinary LO calculations, plus lifetime messages.
// The counters are borrowed and must outlive the helpers.
class TracedFA : public FA_SMLO {
  HelperCounts& counts_;
protected:
  std::unique_ptr<FA_SMLO> clone() const override
  { return std::unique_ptr<FA_SMLO>(new TracedFA(*this)); }
public:
  TracedFA(Fermion type, const inval& input, HelperCounts& counts)
    : FA_SMLO(type, input), counts_(counts) { counts_.create("FA"); }
  TracedFA(const TracedFA& other)
    : FA_SMLO(other), counts_(other.counts_) { counts_.create("FA copy"); }
  ~TracedFA() override { counts_.destroy("FA"); }
};

class TracedSW : public SW_SMLO {
  HelperCounts& counts_;
protected:
  std::unique_ptr<SW_SMLO> clone() const override
  { return std::unique_ptr<SW_SMLO>(new TracedSW(*this)); }
public:
  TracedSW(Fermion type, const inval& input, HelperCounts& counts)
    : SW_SMLO(type, input), counts_(counts) { counts_.create("SW"); }
  TracedSW(const TracedSW& other)
    : SW_SMLO(other), counts_(other.counts_) { counts_.create("SW copy"); }
  ~TracedSW() override { counts_.destroy("SW"); }
};

// This example-only wrapper installs the tracing helpers using the real
// protected FV_SMLO constructor. FV_SMLO itself still does all ownership,
// copying, moving and physics evaluation. It is not a new library interface.
class TracedFV : public FV_SMLO {
  static std::unique_ptr<FA_SMLO> makeAxialFF(
      Fermion type, const inval& input, HelperCounts& counts)
  { return std::unique_ptr<FA_SMLO>(new TracedFA(type, input, counts)); }
  static std::unique_ptr<SW_SMLO> makeSWtheff(
      Fermion type, const inval& input, HelperCounts& counts)
  { return std::unique_ptr<SW_SMLO>(new TracedSW(type, input, counts)); }
public:
  TracedFV(Fermion type, const inval& input, HelperCounts& counts)
    : FV_SMLO(type, input,
              makeAxialFF(type, input, counts), makeSWtheff(type, input, counts)) {}
  // Each factory above returns an already-owned allocation. That is safe
  // whichever constructor argument C++11 evaluates first (FA or SW).
};

void checkSameResult(const FV_SMLO& calculation, Cplx expected)
{
  const Cplx actual = calculation.result();
  if(!std::isfinite(real(actual)) || !std::isfinite(imag(actual))
      || std::abs(actual - expected) > 1e-15)
    throw std::runtime_error("Ownership operation changed the FV result");
}

struct SimulatedInterruption {};

} // namespace

// START HERE: the application, with one step for each ownership operation.
int main()
{
  try
  {
    HelperCounts counts;
    SMval sm;
    sm.set(InputPar::ZMassExperiment, 91.1876);
    sm.set(InputPar::ZWidthExperiment, 2.4952);
    sm.set(InputPar::WMassExperiment, 80.377);
    sm.set(InputPar::WWidthExperiment, 2.085);
    sm.set(InputPar::alpha, 1.0 / 137.03599976);
    std::cout << std::setprecision(10);

    // sm and counts were created first: both outlive all helpers below.
    counts.check(0, 0, "Before any FV object");
    {
      std::cout << "\n1. Construct an electron FV: it owns two helpers.\n";
      TracedFV originalFV(Fermion::electron, sm, counts);
      const Cplx originalValue = originalFV.result();
      std::cout << "  Electron FV = " << originalValue << '\n';
      checkSameResult(originalFV, FV_SMLO(Fermion::electron, sm).result());
      counts.check(2, 0, "After construction");

      std::cout << "\n2. Copy: create another, independently owned pair.\n";
      FV_SMLO copiedFV(originalFV);
      checkSameResult(copiedFV, originalValue);
      counts.check(4, 0, "After copying");

      std::cout << "\n3. Move: transfer the copied pair; do not create another.\n";
      FV_SMLO movedFV(std::move(copiedFV));
      checkSameResult(movedFV, originalValue);
      counts.check(4, 0, "After moving");
      // copiedFV now has no helpers. Do not evaluate it or call its setters.
      // It can safely be destroyed, or assigned a usable configuration.

      std::cout << "\n4. Assign: replace movedFV's pair, releasing its old pair.\n";
      movedFV = originalFV;
      checkSameResult(movedFV, originalValue);
      checkSameResult(originalFV, originalValue);
      counts.check(6, 2, "After replacing the pair");

      std::cout << "\n5. Leave this scope: all remaining owned helpers are freed.\n";
    } // Automatic destruction happens here -- no manual delete.
    counts.check(6, 6, "After leaving the scope");

    std::cout << "\n6. A small W-input scan: no accumulation between iterations.\n";
    const double wMassInputs[] = {80.3, 80.4, 80.5};
    for(int point = 0; point < 3; ++point)
    {
      sm.set(InputPar::WMassExperiment, wMassInputs[point]);
      {
        TracedFV scanFV(Fermion::electron, sm, counts);
        checkSameResult(scanFV, FV_SMLO(Fermion::electron, sm).result());
        std::cout << "  W input = " << wMassInputs[point]
                  << " GeV, FV = " << scanFV.result() << '\n';
        counts.check(8 + 2*point, 6 + 2*point, "Inside this scan step");
      } // The pair is freed before we start the next iteration.
      counts.check(8 + 2*point, 8 + 2*point, "After this scan step");
    }
    // This is an SM input scan, not a SMEFT W-mass-shift prescription.

    std::cout << "\n7. Simulate an error after construction: cleanup still happens.\n";
    try
    {
      TracedFV interruptedFV(Fermion::electron, sm, counts);
      counts.check(14, 12, "Before the simulated interruption");
      throw SimulatedInterruption();
    } // C++ destroys interruptedFV while leaving this scope due to the error.
    catch(const SimulatedInterruption&)
    {
      std::cout << "  Caught the simulated interruption.\n";
    }
    counts.check(14, 14, "After the interrupted calculation");

    // FV owns FA/SW, NOT sm. Destroying all FVs does not destroy the input.
    if(sm.get(InputPar::WMassExperiment) != wMassInputs[2])
      throw std::runtime_error("The borrowed SM input was unexpectedly changed");
    std::cout << "\n  The borrowed SM input is still usable: W input = "
              << sm.get(InputPar::WMassExperiment) << " GeV.\n"
              << "PASS: all 14 traced helpers were destroyed; none remain alive.\n"
              << "These are helper counts, not a process-wide leak scan.\n";
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
