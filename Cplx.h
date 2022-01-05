#ifndef __Cplx__
#define __Cplx__

#include<complex>
#include<iostream>
#include<string>
using namespace std;

#define Double double  // for use with B0 and C0 functions from TVID
typedef std::complex<double> Cplx;

inline Cplx operator + (int x, const Cplx& y)
{ return(Cplx(x)+y); }
inline Cplx operator + (const Cplx& x,int y)
{ return(x+Cplx(y)); }

inline Cplx operator - (int x, const Cplx& y)
{ return(Cplx(x)-y); }
inline Cplx operator - (const Cplx& x,int y)
{ return(x-Cplx(y)); }

inline Cplx operator * (int x, const Cplx& y)
{ return(Cplx(x)*y); }
inline Cplx operator * (const Cplx& x,int y)
{ return(x*Cplx(y)); }

inline Cplx operator / (int x, const Cplx& y)
{ return(Cplx(x)/y); }
inline Cplx operator / (const Cplx& x,int y)
{ return(x/Cplx(y)); }

inline bool operator == (const Cplx& x, int y)
{ return(x==Cplx(y)); }
inline bool operator == (int x, const Cplx& y)
{ return(Cplx(x)==y); }
inline bool operator != (const Cplx& x, int y)
{ return(x!=Cplx(y)); }
inline bool operator != (int x, const Cplx& y)
{ return(Cplx(x)!=y); }

inline Double SDouble(const string& s1) { return(stod(s1,NULL)); }
inline Cplx SCplx(const string& s1, const string& s2) 
{ return(Cplx(stod(s1,NULL),stod(s2,NULL))); }

#define Pi 3.141592653589793
#define TwoPi 6.283185307179587
#define PiSonSix 1.644934066848226  // pi^2/6
#define EULERGAMMA 0.5772156649015329
#define EulerGamma 0.5772156649015329
static const Cplx PiI(0.,Pi);
static const Cplx TwoPiI(0.,TwoPi);
static const Cplx I(0.,1.);

//inline Double sqr(const Double x) { return x*x; }
inline Cplx sqr(const Cplx& x) { return x*x; }

#endif // __Cplx__
