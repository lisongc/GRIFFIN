/*
  File: sm/test_matel_ownership.cc
  Purpose: Verify owning/borrowed matrix-element lifetimes, independent copies,
           moves, exception safety and unchanged v1.1.1 numerical results.
           Adapted from the local v1.2 ownership tests, without v1.2 physics.
  Authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-17
*/

// Numeric-mode values were independently frozen from pre-repair v1.1.1 on 2026-11-17.
// These checks preserve behavior; they are not an independent physics derivation.
#include "griffin/sm/SMval.h"
#include "griffin/sm/xscaas.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

// Count construction/destruction only, excluding numerical work and reporting.
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

void require(bool condition, const char* message)
{
  if(!condition) throw std::runtime_error(message);
}

void close(Cplx actual, Cplx expected, const char* message)
{
  require(std::isfinite(real(actual)) && std::isfinite(imag(actual))
      && std::abs(actual-expected) <= 2e-17+2e-11*std::abs(expected), message);
}

void counts(std::size_t created, std::size_t released, const char* message)
{
  allocation_check::stop();
  require(allocation_check::created == created
      && allocation_check::released == released, message);
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

struct State {
  const psobs* form[4];
  const inval* input;
  Fermion initial, final;
  Current initialCurrent, finalCurrent;
  double s, cosine;
};

// Inspect existing protected views without depending on ownership storage names.
template<typename Matrix>
class Inspect : public Matrix {
public:
  using Matrix::Matrix;
  State state() const
  {
    return {{this->FAi,this->FAo,this->SWi,this->SWo},this->ival,
      this->it,this->ot,this->iff,this->off,this->s,this->cost};
  }
  void changeOwnedValues()
  {
    const psobs* views[] = {this->FAi,this->FAo,this->SWi,this->SWo};
    const double values[] = {0.041,0.042,0.22,0.221};
    for(unsigned i=0; i<4; ++i)
      const_cast<psobsfix*>(static_cast<const psobsfix*>(views[i]))->setval(values[i]);
  }
};

void unchanged(const State& actual, const State& expected)
{
  for(unsigned i=0; i<4; ++i)
    require(actual.form[i] == expected.form[i],"Failed assignment changed a form-factor binding");
  require(actual.input == expected.input && actual.initial == expected.initial
      && actual.final == expected.final && actual.initialCurrent == expected.initialCurrent
      && actual.finalCurrent == expected.finalCurrent && actual.s == expected.s
      && actual.cosine == expected.cosine,"Failed assignment changed matrix state");
}

class TrackedFixed : public psobsfix {
  int& destroyed_;
public:
  TrackedFixed(double value, const inval& input, int& destroyed)
    : psobsfix(value,input), destroyed_(destroyed) {}
  ~TrackedFixed() override { ++destroyed_; }
  const inval* inputAddress() const { return ival; }
};

struct Forms {
  TrackedFixed fai, fao, swi, swo;
  Forms(const inval& input, int& destroyed)
    : fai(0.034,input,destroyed), fao(0.035,input,destroyed),
      swi(0.231,input,destroyed), swo(0.232,input,destroyed) {}
};

template<typename Matrix>
Matrix owned(const inval& input, bool alternate=false)
{
  return Matrix(Fermion::electron,alternate ? Fermion::c : Fermion::muon,
      alternate ? Current::axial : Current::vector,
      alternate ? Current::axial : Current::vector,
      0.034,0.035,0.231,0.232,alternate ? 40000.0 : 10000.0,
      alternate ? -0.2 : 0.3,input);
}

template<typename Matrix>
Matrix borrowed(const inval& input, const Forms& forms)
{
  return Matrix(Fermion::electron,Fermion::muon,Current::vector,Current::vector,
      forms.fai,forms.fao,forms.swi,forms.swo,10000.0,0.3,input);
}

const Cplx references[][2] = {
  {{9.2797675150058188e-06,-1.4741640190865942e-08},
   {-1.0884653179618408e-06,7.8089771575640232e-09}},
  {{9.8677861540393033e-06,-1.6855664212862493e-07},
   {-1.0733508111883957e-06,7.9733000407563842e-08}},
  {{9.8677752969976972e-06,-1.7266930073539898e-07},
   {-1.0735184237437252e-06,8.0464587426664389e-08}}
};

template<typename Matrix>
void checkLifetimeAndValues(const SMval& input, unsigned order)
{
  static_assert(std::is_nothrow_move_constructible<Matrix>::value,
                "Matrix moves must transfer ownership without allocating");
  static_assert(std::is_nothrow_move_assignable<Matrix>::value,
                "Matrix move assignment must not throw");
  for(unsigned point=0; point<2; ++point)
    close(owned<Matrix>(input,point != 0).result(),references[order][point],
          "Numeric constructor changed a frozen matrix-element result");

  for(unsigned repeat=0; repeat<16; ++repeat)
  {
    allocation_check::start();
    { Matrix value=owned<Matrix>(input); }
    counts(4,4,"Numeric matrix construction/destruction must release all four allocations");
  }

  int destroyed=0;
  {
    Forms forms(input,destroyed);
    allocation_check::start();
    { Matrix value=borrowed<Matrix>(input,forms); }
    counts(0,0,"Borrowed matrix construction/destruction must not allocate or delete forms");
    require(destroyed == 0,"Matrix destroyed a borrowed form factor");
    close(borrowed<Matrix>(input,forms).result(),references[order][0],
          "Borrowed and numeric constructors disagree");
  }
  require(destroyed == 4,"External form-factor owners did not retain destruction responsibility");

  psobs* erased=new Matrix(owned<Matrix>(input));
  close(erased->result(),references[order][0],"Polymorphic matrix result changed");
  allocation_check::start();
  delete erased;
  counts(0,5,"Deletion through psobs must release the matrix and its four owned objects");
}

template<typename Matrix>
void checkCopiesAndMoves(const SMval& input)
{
  typedef Inspect<Matrix> M;
  M original=owned<M>(input);
  const Cplx expected=original.result();
  M copy(original);
  for(unsigned i=0; i<4; ++i)
    require(copy.state().form[i] != original.state().form[i],
            "Owning copies must have independent form factors");
  copy.changeOwnedValues();
  require(std::abs(copy.result()-expected)>1e-10,"Copy mutation did not exercise the owned values");
  close(original.result(),expected,"Changing copied form factors changed the original");

  std::unique_ptr<M> survivor;
  {
    M temporary=owned<M>(input);
    survivor.reset(new M(temporary));
  }
  close(survivor->result(),expected,"Copy retained form factors owned by a destroyed source");

  M assigned=owned<M>(input,true);
  allocation_check::start();
  assigned=original;
  counts(4,4,"Owning copy assignment must replace and release the old four objects");
  close(assigned.result(),expected,"Copy assignment changed the result");
  assigned.changeOwnedValues();
  close(original.result(),expected,"Assigned copy shared owned form factors with its source");

  M& self=original;
  allocation_check::start();
  original=self;
  original=std::move(self);
  counts(0,0,"Self-copy/move must preserve owned storage");
  close(original.result(),expected,"Self-assignment changed the result");

  const State before=copy.state();
  const Cplx changed=copy.result();
  allocation_check::start();
  M moved(std::move(copy));
  counts(0,0,"Move construction must transfer existing allocations");
  unchanged(moved.state(),before);
  close(moved.result(),changed,"Move construction changed the result");
  copy=original;
  close(copy.result(),expected,"Moved-from matrix could not be restored by copy assignment");
  allocation_check::start();
  assigned=std::move(moved);
  counts(0,4,"Move assignment must release the target's previous owned objects");
  close(assigned.result(),changed,"Move assignment changed the result");
  moved=owned<M>(input);
  close(moved.result(),expected,"Moved-from matrix could not be restored by move assignment");

  SMval alternate=input;
  alternate.set(InputPar::WMassExperiment,79.9);
  original.setinput(alternate);
  close(original.result(),owned<M>(alternate).result(),"setinput changed fixed form-factor values");
}

template<typename Matrix>
void checkBorrowingAndModeChanges(const SMval& input)
{
  typedef Inspect<Matrix> M;
  int destroyed=0;
  {
    Forms forms(input,destroyed);
    M source=borrowed<M>(input,forms);
    M copy(source);
    unchanged(copy.state(),source.state());
    forms.fai.setval(0.041);
    const Cplx updated=source.result();
    close(copy.result(),updated,"Borrowed copies must observe external form-factor updates");
    require(std::abs(updated-owned<M>(input).result())>1e-10,
            "Borrowed-value mutation did not exercise the result");

    std::unique_ptr<M> survivor;
    {
      M temporary=borrowed<M>(input,forms);
      survivor.reset(new M(temporary));
    }
    close(survivor->result(),updated,"Borrowing copy depended on the source matrix's lifetime");

    M target=owned<M>(input,true);
    allocation_check::start();
    target=source;
    counts(0,4,"Copying borrowed state over owned state must release the old objects");
    unchanged(target.state(),source.state());
    M owning=owned<M>(input);
    allocation_check::start();
    target=owning;
    counts(4,0,"Copying owned state over borrowed state must not delete external objects");
    close(target.result(),owning.result(),"Borrowed-to-owned assignment changed the result");
    require(destroyed == 0,"Assignment deleted an external form factor");

    allocation_check::start();
    target=std::move(source);
    counts(0,4,"Moving borrowed state over owned state must release owned objects only");
    close(target.result(),updated,"Moving borrowed state changed the result");
    source=borrowed<M>(input,forms);
    M movedBorrow(std::move(source));
    close(movedBorrow.result(),updated,"Borrowed move construction changed the result");
    allocation_check::start();
    movedBorrow=std::move(owning);
    counts(0,0,"Moving owned state over borrowed state must transfer ownership without deletion");
    close(movedBorrow.result(),owned<M>(input).result(),"Borrowed-to-owned move changed the result");
    source=borrowed<M>(input,forms);
    M& self=source;
    source=std::move(self);
    close(source.result(),updated,"Borrowed self-move changed the result");

    SMval alternate=input;
    alternate.set(InputPar::WMassExperiment,79.9);
    source.setinput(alternate);
    close(source.result(),borrowed<M>(alternate,forms).result(),"Borrowed setinput changed form factors");
    require(forms.fai.inputAddress() == &input && forms.fao.inputAddress() == &input
        && forms.swi.inputAddress() == &input && forms.swo.inputAddress() == &input,
        "setinput must not mutate external form factors");
    require(destroyed == 0,"A borrowed form factor was destroyed before its external owner");
  }
  require(destroyed == 4,"Borrowed objects were not destroyed exactly once by their external owner");
}

template<typename Matrix>
void checkAllocationFailures(const SMval& input)
{
  typedef Inspect<Matrix> M;
  M source=owned<M>(input);
  SMval alternate=input;
  alternate.set(InputPar::WMassExperiment,79.9);
  int destroyed=0;
  Forms forms(alternate,destroyed);
  for(int fail=0; fail<4; ++fail)
  {
    bool threw=false;
    allocation_check::start();
    allocation_check::allocationsBeforeFailure=fail;
    try { M value=owned<M>(input); }
    catch(const std::bad_alloc&) { threw=true; }
    counts(fail,fail,"Failed construction leaked an earlier allocation");
    require(threw,"Construction did not reach the requested allocation failure");

    threw=false;
    allocation_check::start();
    allocation_check::allocationsBeforeFailure=fail;
    try { M copy(source); }
    catch(const std::bad_alloc&) { threw=true; }
    counts(fail,fail,"Failed copy construction leaked an earlier allocation");
    require(threw,"Copy construction did not reach the requested allocation failure");

    for(unsigned mode=0; mode<2; ++mode)
    {
      M target=mode ? borrowed<M>(alternate,forms) : owned<M>(alternate,true);
      const State before=target.state();
      const Cplx beforeValue=target.result();
      threw=false;
      allocation_check::start();
      allocation_check::allocationsBeforeFailure=fail;
      try { target=source; }
      catch(const std::bad_alloc&) { threw=true; }
      counts(fail,fail,"Failed copy assignment leaked or deleted existing target storage");
      require(threw,"Copy assignment did not reach the requested allocation failure");
      unchanged(target.state(),before);
      close(target.result(),beforeValue,"Failed copy assignment changed the target result");
      require(destroyed == 0,"Failed assignment deleted borrowed objects");
    }
  }
}

template<typename Matrix>
void checkMatrix(const SMval& input, unsigned order)
{
  checkLifetimeAndValues<Matrix>(input,order);
  checkCopiesAndMoves<Matrix>(input);
  checkBorrowingAndModeChanges<Matrix>(input);
  checkAllocationFailures<Matrix>(input);
}

} // namespace

int main()
{
  try {
    const SMval input=makeInput();
    checkMatrix<matel>(input,0);
    checkMatrix<mat_SMNNLO>(input,1);
    checkMatrix<mat_SMaas>(input,2);
    std::cout << "matel ownership, copy/move and numerical regression checks passed\n";
    return 0;
  }
  catch(const std::exception& error) {
    allocation_check::stop();
    std::cerr << error.what() << '\n';
    return 1;
  }
}
