/*
  File: classes.h
  Purpose: Define typed SM input, fermion and current labels and the base
           input, form-factor and scattering classes. Internally created form
           factors are owned and copied safely; external inputs remain borrowed.
           Integer labels and header-dependent mass macros are not supported.
  Merged-version authors: Jackson Wallace and Lisong Chen
  Last updated: 2026-11-17
*/

#ifndef __classes__
#define __classes__

#include <vector>
#include <memory>
#include <utility>
#include <math.h>
#include <iostream>
using namespace std;
#include "griffin/core/Cplx.h"

namespace griffin {

/******************************************************************************
 inval: base class for input parameters
 There can be derived versions that, e.g., compute MW from Gmu.
******************************************************************************/
#define SIZE1 100 // default for data array

enum class InputPar {
  WMassComplex = 0,
  ZMassComplex = 1,
  HiggsMass = 2,
  electronMass = 3,
  muonMass = 4,
  tauMass = 5,
  downMass = 6,
  strangeMass = 7,
  bottomMass = 8,
  upMass = 9,
  charmMass = 10,
  topMass = 11,
  alpha = 12,
  alphaS = 13,
  DeltaAlpha = 14,
  DeltaAlphaHad = 15,
  fermiConstant = 16,
  WWidthComplex = 17,
  ZWidthComplex = 18,
  WMassExperiment = 20,
  ZMassExperiment = 21,
  WWidthExperiment = 22,
  ZWidthExperiment = 23
};

inline int index(InputPar par)
{
  return static_cast<int>(par);
}

class inval {
protected:
  std::vector<double> data;
  virtual void compute(void) {};  // compute some input pars. in terms of others

public:
  inval(const int sizealloc) : data(sizealloc, NAN) {};
  inval(void) : inval(SIZE1) {};
  inval(const inval& copyfrom) : data(copyfrom.data) {};
 
  unsigned int getsize() const
  { return data.size(); }
     
  void set(const InputPar par, const double val)
  {
    const int idx = index(par);
    if(idx < data.size())
      data[idx] = val;
    else
    {
      cerr << "Input value index " << idx << " outside of range" << endl;
      exit(1);
    }
    compute();
  }

  double get(const InputPar par) const
  {
    const int idx = index(par);
    if(idx < data.size())
    {
      if(isfinite(data[idx]))
        return(data[idx]);
    }
    cerr << "Invalid or undefined input value for index " << idx << endl;
    exit(1);
  }

};

/******************************************************************************
 psobs: base class for a pseudo-observable or form factor
******************************************************************************/
class psobs {
protected:
  // Borrowed input; the caller must keep it alive.
  const inval* ival;
public:
  psobs(const inval& input) { ival = &input; }

  // Release derived objects correctly when deleted through psobs/FA/SW pointers.
  virtual ~psobs() = default;
  virtual void setinput(const inval& input) { ival = &input; }
  virtual Cplx result(void) const = 0;
  virtual Cplx errest(void) const { return 0; }
};

// dummy object that returns a fixed value
class psobsfix : public psobs {
protected:
  double myval;
public:
  psobsfix(const double value, const inval& input) : psobs(input)
  { myval = value; }
  
  void setval(const double value) { myval = value; }
  Cplx result(void) const override { return myval; }
};

// Particle codes used by the internal charge/isospin tables.

enum class Fermion {
  d = 1,
  u = 2,
  s = 3,
  c = 4,
  b = 5,
  electron = 11,
  nuElectron = 12,
  muon = 13,
  nuMuon = 14,
  tau = 15,
  nuTau = 16
};

inline int index(Fermion type)
{
  return static_cast<int>(type);
}

inline bool isNeutrino(Fermion type)
{
  return type == Fermion::nuElectron
      || type == Fermion::nuMuon
      || type == Fermion::nuTau;
}

// Scalar and pseudoscalar currents are relevant only for Bhabha scattering.

enum class Current {
  vector = 0,
  axial = 1,
  scalar = 2,
  pseudoscalar = 3
};

inline int index(Current current)
{
  return static_cast<int>(current);
}

// object for matrix element computation
class matel : public psobs {
private:
  // Present only for numerical-value construction. External form factors
  // supplied through the other constructor are never owned or deleted here.
  std::unique_ptr<psobsfix> ownedFAi, ownedFAo, ownedSWi, ownedSWo;
protected:
  Fermion it, ot;
  Current iff, off;
  double s, cost;
  // Non-owning views of either our fixed values or external form factors.
  const psobs *FAi, *FAo, *SWi, *SWo;
public:
  // Own fixed form factors made from the supplied numerical values.
  matel(const Fermion intype, const Fermion outtype,
        const Current inform, const Current outform,
        const double FAin, const double FAout, const double SWin, const double SWout,
        const double sval, const double costheta, const inval& input)
    : psobs(input),
      ownedFAi(new psobsfix(FAin,input)), ownedFAo(new psobsfix(FAout,input)),
      ownedSWi(new psobsfix(SWin,input)), ownedSWo(new psobsfix(SWout,input)),
      it(intype), ot(outtype), iff(inform), off(outform), s(sval), cost(costheta),
      FAi(ownedFAi.get()), FAo(ownedFAo.get()),
      SWi(ownedSWi.get()), SWo(ownedSWo.get()) {}

  // Borrow existing form factors; they and input must outlive this object.
  matel(const Fermion intype, const Fermion outtype,
        const Current inform, const Current outform,
        const psobs& FAin, const psobs& FAout, const psobs& SWin, const psobs& SWout,
        const double sval, const double costheta, const inval& input)
    : psobs(input), it(intype), ot(outtype), iff(inform), off(outform),
      s(sval), cost(costheta), FAi(&FAin), FAo(&FAout), SWi(&SWin), SWo(&SWout) {}

  // Copies duplicate owned values, while external form factors remain borrowed.
  matel(const matel& other)
    : psobs(other),
      ownedFAi(other.ownedFAi ? new psobsfix(*other.ownedFAi) : nullptr),
      ownedFAo(other.ownedFAo ? new psobsfix(*other.ownedFAo) : nullptr),
      ownedSWi(other.ownedSWi ? new psobsfix(*other.ownedSWi) : nullptr),
      ownedSWo(other.ownedSWo ? new psobsfix(*other.ownedSWo) : nullptr),
      it(other.it), ot(other.ot), iff(other.iff), off(other.off),
      s(other.s), cost(other.cost),
      FAi(ownedFAi ? ownedFAi.get() : other.FAi),
      FAo(ownedFAo ? ownedFAo.get() : other.FAo),
      SWi(ownedSWi ? ownedSWi.get() : other.SWi),
      SWo(ownedSWo ? ownedSWo.get() : other.SWo) {}

  matel& operator=(const matel& other)
  {
    if(this != &other)
    {
      // Complete all allocations before changing the destination.
      matel copy(other);
      *this = std::move(copy);
    }
    return *this;
  }

  // Moves transfer ownership. Reassign the source before evaluating it.
  matel(matel&& other) noexcept
    : psobs(other),
      ownedFAi(std::move(other.ownedFAi)), ownedFAo(std::move(other.ownedFAo)),
      ownedSWi(std::move(other.ownedSWi)), ownedSWo(std::move(other.ownedSWo)),
      it(other.it), ot(other.ot), iff(other.iff), off(other.off),
      s(other.s), cost(other.cost),
      FAi(other.FAi), FAo(other.FAo), SWi(other.SWi), SWo(other.SWo)
  {
    other.FAi = other.FAo = other.SWi = other.SWo = nullptr;
  }

  matel& operator=(matel&& other) noexcept
  {
    if(this != &other)
    {
      psobs::operator=(other);
      ownedFAi = std::move(other.ownedFAi);
      ownedFAo = std::move(other.ownedFAo);
      ownedSWi = std::move(other.ownedSWi);
      ownedSWo = std::move(other.ownedSWo);
      it = other.it; ot = other.ot; iff = other.iff; off = other.off;
      s = other.s; cost = other.cost;
      FAi = other.FAi; FAo = other.FAo; SWi = other.SWi; SWo = other.SWo;
      other.FAi = other.FAo = other.SWi = other.SWo = nullptr;
    }
    return *this;
  }

  void setinput(const inval& input) override
  {
    psobs::setinput(input);
    if(ownedFAi) ownedFAi->setinput(input);
    if(ownedFAo) ownedFAo->setinput(input);
    if(ownedSWi) ownedSWi->setinput(input);
    if(ownedSWo) ownedSWo->setinput(input);
  }

  // Change initial-state and final-state currents.
  void setform(const Current inform, const Current outform)
  {
    iff = inform; off = outform;
  }

  // change initial-state and final-state fermion types
  void setftype(const Fermion intype, const Fermion outtype)
  {
    it = intype; ot = outtype;
  }

  // change kinematic variables
  void setkinvar(const double sval, const double costheta)
  {
    s = sval; cost = costheta;
  }

  /* in the base version of matel, the following functions implement tree-level
     expressions; see classes.cc for the code */
  Cplx coeffR(void) const;  // R coefficient of complex pole expansion
  Cplx coeffS(void) const;  // S coefficient of complex pole expansion
  Cplx coeffSp(void) const; // S' coefficient of complex pole expansion
  Cplx resoffZ(void) const; // off-resonance contribution, M^{noexp} - M^{R,mz^2}
  Cplx result(void) const override; // total result (complex pole exp.+ off-res. piece)
};


// form factor for sw_eff predicted in the SM (at LO)
class SW_SMLO : public psobs {
protected:
  Fermion ftyp;
  friend class FV_SMLO;
  // FV copies its helper without losing the derived correction order.
  // The copy still borrows the same SM input.
  virtual std::unique_ptr<SW_SMLO> clone() const
  { return std::unique_ptr<SW_SMLO>(new SW_SMLO(*this)); }
public:
  SW_SMLO(const Fermion type, const inval& input) : psobs(input)
  {
    ftyp = type;
  }
  void setftype(const Fermion type) // change fermion type
  {
    ftyp = type;
  }
  
  Cplx result(void) const override;  // see classes.cc for code
};

// form factor for F_A predicted in the SM (at LO)
class FA_SMLO : public psobs {
protected:
  Fermion ftyp;
  friend class FV_SMLO;
  // Independent helper, same borrowed SM input; derived versions override this.
  virtual std::unique_ptr<FA_SMLO> clone() const
  { return std::unique_ptr<FA_SMLO>(new FA_SMLO(*this)); }
public:
  FA_SMLO(const Fermion type, const inval& input) : psobs(input)
  {
    ftyp = type;
  }
  void setftype(const Fermion type) // change fermion type
  {
    ftyp = type;
  }
    
  Cplx result(void) const override;  // see classes.cc for code
};

// F_V owns its FA/SW helpers; unique_ptr releases them automatically.
// Copies get independent helpers of the same correction order, while the
// SM input remains borrowed. Moves transfer ownership of the helpers.
class FV_SMLO : public psobs {
protected:
  Fermion ftyp;
  std::unique_ptr<FA_SMLO> fa;
  std::unique_ptr<SW_SMLO> sw;

  // Create a form-factor object, and return an owner that will automatically
  // delete it when it is no longer needed.
  template<typename FormFactor>
  static std::unique_ptr<FormFactor> makeFormFactor(const Fermion type,
                                                   const inval& input)
  { return std::unique_ptr<FormFactor>(new FormFactor(type, input)); }

  // Accept the final-order pair directly, avoiding unused lower-order pairs.
  FV_SMLO(const Fermion type, const inval& input,
          std::unique_ptr<FA_SMLO> axiaff, std::unique_ptr<SW_SMLO> swtheff)
    : psobs(input), ftyp(type), fa(std::move(axiaff)), sw(std::move(swtheff)) {}
public:

  FV_SMLO(const Fermion type, const inval& input)
    : FV_SMLO(type, input,
              makeFormFactor<FA_SMLO>(type, input),  // Create the owned FA helper.
              makeFormFactor<SW_SMLO>(type, input)) {}

  // Make a new FV from other, with independent copies of its FA/SW helpers.
  FV_SMLO(const FV_SMLO& other)
    : psobs(other), ftyp(other.ftyp),
      fa(other.fa ? other.fa->clone() : nullptr),   // A missing helper stays empty.
      sw(other.sw ? other.sw->clone() : nullptr) {}
  // Update an FV object that already exists.
  FV_SMLO& operator=(const FV_SMLO& other)
  {
    if(this != &other)
    {
      // Finish both copies first; a failure leaves the destination unchanged.
      std::unique_ptr<FA_SMLO> axiaff = other.fa ? other.fa->clone() : nullptr;
      std::unique_ptr<SW_SMLO> swtheff = other.sw ? other.sw->clone() : nullptr;
      psobs::operator=(other);
      ftyp = other.ftyp;
      fa = std::move(axiaff);
      sw = std::move(swtheff);
    }
    return *this;
  }
  // Reassign a moved-from FV before evaluating it or calling its setters.
  FV_SMLO(FV_SMLO&&) noexcept = default;
  FV_SMLO& operator=(FV_SMLO&&) noexcept = default;

  void setinput(const inval& input) override
  {
    psobs::setinput(input);
    fa->setinput(input);
    sw->setinput(input);
  }
  void setftype(const Fermion type) // change fermion type
  {
    ftyp = type;
    fa->setftype(type);
    sw->setftype(type);
  }
  
  Cplx result(void) const override;  // see classes.cc for code
  Cplx errest(void) const override;  // see classes.cc for code
};

} // namespace

#endif // __classes__
