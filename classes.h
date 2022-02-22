/* classes.h: header file for classes in classes.cc */

#ifndef __classes__
#define __classes__

#include <vector>
#include <math.h>
#include <iostream>
using namespace std;
#include "Cplx.h"
#include "s2lseinline.h"

/******************************************************************************
 inval: base class for input parameters
 There can be derived versions that, e.g., compute MW from Gmu.
******************************************************************************/
#define SIZE1 100 // default for data array
#define MW 0
#define MZ 1
#define MH 2
#define ME 3 // 3-10 are currently not used but included for completeness
#define MM 4
#define ML 5
#define MD 6
#define MS 7  
#define MB 8
#define MU 9
#define MC 10
#define MT 11
#define al 12
#define als 13
#define Delal 14
#define Delalhad 15
#define Gmu 16
#define GamW 17
#define GamZ 18

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
     
  void set(const int idx, const double val)  
  {
    if(idx < data.size())
      data[idx] = val;
    else
    {
      cerr << "Input value index " << idx << " outside of range" << endl;
      exit(1);
    }
    compute();
  }
  
  double get(const int idx) const
  {
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
  const inval* ival;
public:
  psobs(const inval& input) { ival = &input; }
  virtual void setinput(const inval& input) { ival = &input; }
  virtual Cplx result(void) const = 0;
  virtual double errest(void) { return 0; }  // should this be complex ??
};

// use basic relative factor al/(Pi*SWs)*(3*3+3) for error estimate?
// -> could be later replaced by something better...


// dummy object that returns a fixed value
class psobsfix : public psobs {
protected:
  double myval;
public:
  psobsfix(const double value, const inval& input) : psobs(input)
  { myval = value; }
  
  Cplx result(void) const { return myval; }
};

// possible values for intype/outtype below
#define LEP 0
#define NEU 1
#define UQU 2
#define DQU 3
#define BQU 4

// possible values for inform/outform below
#define VEC 0
#define AXV 1

// objects for matrix element coefficients R, S, S', ... in Z-pole expansion
class coeffR : public psobs {
protected:
  int it, ot, iff, off;
  const psobs *FAi, *FAo, *SWi, *SWo;
public:
  // this constructor initiates the object with user-provided values for the 
  // Zff vertex form factor -- it can be used when, e.g., trying to make a fit
  // of these form factors to the data
  coeffR(const int intype, const int outtype, const int inform, const int outform, 
  	 const double FAin, const double FAout, const double SWin, const double SWout,
	 const inval& input) : psobs(input)
  {
    it = intype; ot = outtype; iff = inform; off = outform;
    FAi = new psobsfix(FAin,input);
    FAo = new psobsfix(FAout,input);
    SWi = new psobsfix(SWin,input);
    SWo = new psobsfix(SWout,input);
  }

  // this constructor takes objects are inputs for the Zff vertex form factors;
  // these objects are supposed to compute predictions for the form factors
  // within the SM or in some BSM model
  coeffR(const int intype, const int outtype, const int inform, const int outform, 
  	 const psobs& FAin, const psobs& FAout, const psobs& SWin, const psobs& SWout,
	 const inval& input) : psobs(input)
  {
    it = intype; ot = outtype; iff = inform; off = outform;
    FAi = &FAin; FAo = &FAout; SWi = &SWin; SWo = &SWout;
  }
  
  Cplx result(void) const;  // see classes.cc for code
};

class coeffS : public psobs {
protected:
  int it, ot, iff, off;
public:
  coeffS(const int intype, const int outtype, const int inform, const int outform, 
  	 const inval& input) : psobs(input)
  {
    it = intype; ot = outtype; iff = inform; off = outform;
  }
  
  Cplx result(void) const;  // see classes.cc for code
};

class coeffSp : public psobs {
protected:
  int it, ot, iff, off;
public:
  coeffSp(const int intype, const int outtype, const int inform, const int outform, 
  	  const inval& input) : psobs(input)
  {
    it = intype; ot = outtype; iff = inform; off = outform;
  }
  
  Cplx result(void) const;  // see classes.cc for code
};


// squared matrix element directly computed rather than from R,S,S'
// -> see section 4 of my notes
class matelsq : public psobs {
protected:
  int it, ot, if1, of1, if2, of2;
  double s, cost;
  const psobs *FAi, *FAo, *SWi, *SWo;
public:
  // this constructor initiates the object with user-provided values for the 
  // Zff vertex form factor -- it can be used when, e.g., trying to make a fit
  // of these form factors to the data
  matelsq(const int intype, const int outtype, 
          const int inform1, const int outform1, 
          const int inform2, const int outform2, 
  	  const double FAin, const double FAout, const double SWin, const double SWout,
	  const double sval, const double costheta, const inval& input) : 
	  psobs(input)
  {
    it = intype; ot = outtype; 
    if1 = inform1; of1 = outform1; if2 = inform2; of2 = outform2; 
    FAi = new psobsfix(FAin,input);
    FAo = new psobsfix(FAout,input);
    SWi = new psobsfix(SWin,input);
    SWo = new psobsfix(SWout,input);
    s = sval; cost = costheta;
  }

  // this constructor takes objects are inputs for the Zff vertex form factors;
  // these objects are supposed to compute predictions for the form factors
  // within the SM or in some BSM model
  matelsq(const int intype, const int outtype, 
          const int inform1, const int outform1, 
          const int inform2, const int outform2, 
  	  const psobs& FAin, const psobs& FAout, const psobs& SWin, const psobs& SWout,
	  const double sval, const double costheta, const inval& input) : 
	  psobs(input)
  {
    it = intype; ot = outtype; 
    if1 = inform1; of1 = outform1; if2 = inform2; of2 = outform2; 
    FAi = &FAin; FAo = &FAout; SWi = &SWin; SWo = &SWout;
    s = sval; cost = costheta;
  }
  
  Cplx result(void) const;  // see classes.cc for code
};


// the form factors predicted in the SM (at LO)
class SW_SMLO : public psobs {
protected:
  int ftyp;
public:
  SW_SMLO(const int type, const inval& input) : psobs(input)
  {
    ftyp = type;
  }
  void setftype(const int type)
  {
    ftyp = type;
  }
  
  Cplx result(void) const;  // see classes.cc for code
};

class FA_SMLO : public psobs {
protected:
  int ftyp;
public:
  FA_SMLO(const int type, const inval& input) : psobs(input)
  {
    ftyp = type;
  }
  void setftype(const int type)
  {
    ftyp = type;
  }
    
  Cplx result(void) const;  // see classes.cc for code
};

class FV_SMLO : public psobs {
protected:
  int ftyp;
  FA_SMLO *fa;
  SW_SMLO *sw;
public:
  FV_SMLO(const int type, const inval& input) : psobs(input)
  {
    ftyp = type;
    fa = new FA_SMLO(type, input);
    sw = new SW_SMLO(type, input);
  }
  void setinput(const inval& input)
  {
    psobs::setinput(input);
    fa->setinput(input);
    sw->setinput(input);
  }
  void setftype(const int type)
  {
    ftyp = type;
    fa->setftype(type);
    sw->setftype(type);
  }
  
  Cplx result(void) const;  // see classes.cc for code
};

#endif // __classes__
