

#include <math.h>
#include <iostream>
using namespace std;
#include "Cplx.h"

/******************************************************************************
 inval: base class for input parameters
 There can be derived versions that, e.g., compute MW from Gmu.
******************************************************************************/
class inval
{
  protected:
   int size;
    double *data;                 // me:protected data members can be accessed by directly drived classes
  // using "use std :: vector<>" here? which uses heap, and deletes properly when going out of scope
   virtual void compute(void){}; // compute some input pars. in terms of others

public:
  //int size;
  //double *data;
  //virtual void compute(void){};

  inval(const int sizealloc) //me:constructor initiates input params
  {
    int i;
    size = sizealloc;
    data = new double[size]; //allocate memory of a size-element array to the pointer data, return an address of the first element
    for (i = 0; i <= size; i++)
      data[i] = NAN; //initialize the block of memory allocated by NAN
  }

  ~inval() { delete[] data; }
  //define getter and setter of the class to access the private data
  void set(const int idx, const double val)
  {
    if (idx < size)
      data[idx] = val; //assign value to each element by input val from object
    else
    {
      cerr << "Input value index " << idx << " outside of range" << endl;
      exit(1);
    }
  }

  double get(const int idx) const //what is the difference between get() and set()?
                                  //get() returns a const variable, set() assigns values to data[].?
  {
    if (idx < size)
    {
      if (isnormal(data[idx]))
        return (data[idx]);
    }
    cerr << "Invalid or undefined input value for index " << idx << endl;
    exit(1);
  }
};

#define StanMod 20 // for sizealloc parameter
#define MW 0
#define MZ 1
#define MH 2
#define ME 3 // 3-10 are currently not used but included for completeness
#define MM 4
#define ML 5
#define MD 6
#define MU 7
#define MS 8
#define MC 9
#define MB 10
#define MT 11
#define al 12
#define als 13
#define Delal 14
#define Delalhad 15
#define Gmu 16
#define GamZ 17

/******************************************************************************
 psobs: base class for a pseudo-observable or form factor
 There can be derived versions that, e.g., compute MW from Gmu.

 me: this is an abstract base class that one cannot define an object of it. (
 the pure virtual functions return zeros.) What it mainly does is to define pointers pointing to
 objects of its derived classed?
******************************************************************************/
class psobs
{
  //protected:
public:
  const inval *ival; // define a pointer ival pointing to inval constant, this is protected, may not be necessary??
                     //public:
  // const inval *ival;
  psobs(const inval &input) { ival = &input; }
  virtual Cplx result(void) const = 0;
  virtual double errest(void) { return 0; } // should this be complex ??
};

// use basic relative factor al/(Pi*SWs)*(3*3+3) for error estimate?
// -> could be later replaced by something better...

// dummy object that returns a fixed value. me:I am still a bit of confused abt what this cls does...
class psobsfix : public psobs
{
protected:
  double myval;

public:
  psobsfix(const double value, const inval &input) : psobs(input)
  {
    myval = value;
  } //me: initialize psobsfix

  Cplx result(void) const { return myval; } // me: why a result() has to be defined in here?
                                            //(unlike other polymorphically defined in classes.cc)
};

// possible values for intype/outtype below
#define LEP 0
#define NEU 1
#define UQU 2
#define DQU 3
#define BQU 4
#define CQU 5
#define SQU 6

// possible values for inform/outform below
#define VEC 0
#define AXV 1
// possible vaules for polarized basis
#define VV 0
#define VA 1
#define AV 2
#define AA 3
//possible values for in/out gauge boson in self energy or boxes
#define ZZWW 0
#define GZ 1
#define GG 2
//possible values for perturbative order
#define LO 0
#define NLO 1
#define NNLO 2

//possible values for amp



// objects for matrix element coefficients R, S, S', ... in Z-pole expansion
class coeffR : public psobs
{
  int it, ot, iff, off;
  const psobs *FAi, *FAo, *SWi, *SWo;

public:
  // this constructor initiates the object with user-provided values for the
  // Zff vertex form factor -- it can be used when, e.g., trying to make a fit
  // of these form factors to the data
  coeffR(const int intype, const int outtype, const int inform, const int outform,
         const double FAin, const double FAout, const double SWin, const double SWout,
         const inval &input) : psobs(input)
  {
    it = intype;
    ot = outtype;
    iff = inform;
    off = outform;
    FAi = new psobsfix(FAin, input);
    FAo = new psobsfix(FAout, input);
    SWi = new psobsfix(SWin, input);
    SWo = new psobsfix(SWout, input);
  }

  // this constructor takes objects are inputs for the Zff vertex form factors;
  // these objects are supposed to compute predictions for the form factors
  // within the SM or in some BSM model
  coeffR(const int intype, const int outtype, const int inform, const int outform,
         const psobs &FAin, const psobs &FAout, const psobs &SWin, const psobs &SWout,
         const inval &input) : psobs(input)
  {
    it = intype;
    ot = outtype;
    iff = inform;
    off = outform;
    FAi = &FAin;
    FAo = &FAout;
    SWi = &SWin;
    SWo = &SWout;
  }

  Cplx result(void) const; // what this actually does will be written elsewhere
};

class coeffS : public psobs
{
  int it, ot, iff, off;

public:
  coeffS(const int intype, const int outtype, const int inform, const int outform,
         const inval &input) : psobs(input)
  {
    it = intype;
    ot = outtype;
    iff = inform;
    off = outform;
  }

  Cplx result(void) const; // what this actually does will be written elsewhere
};

class coeffSp : public psobs
{
  int it, ot, iff, off;

public:
  coeffSp(const int intype, const int outtype, const int inform, const int outform,
          const inval &input) : psobs(input)
  {
    it = intype;
    ot = outtype;
    iff = inform;
    off = outform;
  }

  Cplx result(void) const; // what this actually does will be written elsewhere
};

// squared matrix element directly computed rather than from R,S,S'
// -> see section 4 of my notes
class matelsq : public psobs
{
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
          const double sval, const double costheta, const inval &input) : psobs(input)
  {
    it = intype;
    ot = outtype;
    if1 = inform1;
    of1 = outform1;
    if2 = inform2;
    of2 = outform2;
    FAi = new psobsfix(FAin, input);
    FAo = new psobsfix(FAout, input);
    SWi = new psobsfix(SWin, input);
    SWo = new psobsfix(SWout, input);
    s = sval;
    cost = costheta;
  }

  // this constructor takes objects are inputs for the Zff vertex form factors;
  // these objects are supposed to compute predictions for the form factors
  // within the SM or in some BSM model
  matelsq(const int intype, const int outtype,
          const int inform1, const int outform1,
          const int inform2, const int outform2,
          const psobs &FAin, const psobs &FAout, const psobs &SWin, const psobs &SWout,
          const double sval, const double costheta, const inval &input) : psobs(input)
  {
    it = intype;
    ot = outtype;
    if1 = inform1;
    of1 = outform1;
    if2 = inform2;
    of2 = outform2;
    FAi = &FAin;
    FAo = &FAout;
    SWi = &SWin;
    SWo = &SWout;
    s = sval;
    cost = costheta;
  }

  Cplx result(void) const; // what this actually does will be written elsewhere
};

// the form factors predicted in the SM (at LO)
class FA_SMLO : public psobs
{
  int ftyp;

public:
  FA_SMLO(const int type, const inval &input) : psobs(input)
  {
    ftyp = type;
  }

  Cplx result(void) const; // what this actually does will be written elsewhere
};

/*******lisong's try***************************************************************/

// the form factors of Z, Z', and Z'' predicted in the SM (at NLO)
class Z_SMNLO : public psobs
{
  int ftyp, off;

public:
  Z_SMNLO(const int type, const int formt, const inval &input) : psobs(input)
  {
    off = formt;
    ftyp = type;
  }

  Cplx result(void) const; // computing form factor of Z
};

class Zp_SMNLO : public psobs
{
  int ftyp, off;

public:
  Zp_SMNLO(const int type, const int formt, const inval &input) : psobs(input)
  {
    off = formt;
    ftyp = type;
  }

  Cplx result(void) const; // computing form factor of Z
};

// the form factors of G, G',G'',... predicted in the SM (at NLO)
class G_SMNLO : public psobs
{
  int ftyp, off, odd;

public:
  G_SMNLO(const int ord, const int type, const int formt, const inval &input) : psobs(input)
  {
    off = formt;
    ftyp = type;
    odd = ord;
  }
  Cplx result(void) const; //computing form factor of G
};

class Gp_SMNLO : public psobs
{
  int ftyp, off;

public:
  Gp_SMNLO(const int type, const int formt, const inval &input) : psobs(input)
  {
    off = formt;
    ftyp = type;
  }
  Cplx result(void) const; //computing form factor of G
};

// The Box contributions to R, S, and S'.

class b_SM : public psobs
{
  int ot, iff, off;

public:
  b_SM(const int outtype, const int inform, const int outform, const inval &input) : psobs(input), off(outform), iff(inform), ot(outtype) {}

  Cplx result(void) const; //computing box contributions
};

class SW_SMLO : public psobs
{
  int ftyp;

public:
  SW_SMLO(const int type, const inval &input) : psobs(input), ftyp(type) {}

  Cplx result(void) const; // what this actually does will be written elsewhere
};

class FA_SMNLO : public FA_SMLO
{
  int ftp;

public:
  FA_SMNLO(const int type, const inval &input) : FA_SMLO(type, input), ftp(type)
  {
  }

  Cplx result(void) const;
};

class SW_SMNLO : public SW_SMLO
{ 
  int ftp;
  public:
  SW_SMNLO(const int type, const inval &input) : SW_SMLO(type, input), ftp(type){}
  
  Cplx result(void) const;
};



