/* ------------------------------------------------------------- */
/* B0.h, B0.cc                                                   */
/* Stefan Bauberger (stefan@bauberger.de)                        */
/* last revision 4.12.95                                         */
/* ------------------------------------------------------------- */
/* scalar one-loop function B0, some derivatives and the         */
/* discontinuity, for use with the package s2lse                 */
/* ------------------------------------------------------------- */

#include "B0.h"
#include "Cplx.h"


#define Errstream cerr

// one-loop self-energy,  B0A = B0 -1/delta +EulerGamma -log(4 Pi mu^2) -2
Cplx B0A(Double ps, Double m1s, Double m2s)
{
  if (m1s==Double(0)) interchange(&m1s,&m2s);
  if ((ps==Double(0))&&(m2s==Double(0))) return(-1-log(m1s));
  if ((ps==Double(0))&&(m1s!=m2s))
    return(m2s*log(m2s)/(m1s-m2s)+m1s*log(m1s)/(m2s-m1s)-1);
  if ((ps==Double(0))&&(m1s==m2s))
    return(-2-log(m1s));
  if ((m2s==Double(0))&&(m1s==Double(0)))
    return(-clog2(-ps));
  if ((m2s==Double(0))&&(m1s==ps))
    return(-clog1(ps));
  if (m2s==Double(0))
    return(-(m1s/ps)*log(m1s)+((m1s-ps)/ps)*clog2(m1s-ps));
  if (abs(ps)*ASYMP_LIMIT_B0S < (m1s+m2s) )
    return(m2s*log(m2s)/(m1s-m2s)+m1s*log(m1s)/(m2s-m1s)-1
           +ps*DB0(0.,m1s,m2s));
  Cplx r=-log(m1s*m2s)/2+(m1s-m2s)/ps*log(m2s/m1s)/2;
  Double s,d,ss,ds,la,r1,r2;
  s=sqrt(m1s)+sqrt(m2s);
  ss=s*s;
  d=sqrt(m1s)-sqrt(m2s);
  ds=d*d;
  if (ps<=ds)
  {
    la=sqrt((ss-ps)*(ds-ps));
    r+=la*log((m1s+m2s-ps+la)/(2*sqrt(m1s*m2s)))/ps;
  }
  else if (ps>=ss)
    {
      la=sqrt((ss-ps)*(ds-ps));
      r-=la*(log((-m1s-m2s+ps+la)/(2*sqrt(m1s*m2s)))-PiI)/ps;
    }
    else
    {
      r1=sqrt(ps-ds);
      r2=sqrt(ss-ps);
      r-=2.*r1*r2*atan(r1/r2)/ps;
    };
  return(r);
}

Cplx B0A(Cplx ps, Double m1s, Double m2s)
{
  Cplx r=-log(m1s*m2s)/2+(m1s-m2s)/ps*log(m2s/m1s)/2;
  Double s,d,ss,ds,r1,r2;
  s=sqrt(m1s)+sqrt(m2s);
  ss=s*s;
  d=sqrt(m1s)-sqrt(m2s);
  ds=d*d;
  Cplx Cla=sqrt(ss-ps)*sqrt(ds-ps);
  r+=Cla*log((m1s+m2s-ps+Cla)/(2*sqrt(m1s*m2s)))/ps;
  return(r);
}

// with one vanishing mass:
Cplx B0MZ(Double ps, Double ms)
{
  if (ps==Double(0))        // TO IMPROVE: make this condition weaker!
    return(-1.-log(ms));
  else
    return(-(ms/ps)*log(ms)+((ms-ps)/ps)*clog2(ms-ps));
}

// subtracted B0(ps,m1s,m2s)-B0(ps,m1s,0):
Cplx B0N1(Double ps, Double m1s, Double m2s)
{
  if ((m1s>abs(ps)*ASYMP_LIMIT_B0N1)&&(ps != Double(0)))
    return(log(m1s-ps)-
           ((m1s-ps)/(m1s-ps-m2s))*log(m1s-ps)+
           ((m2s)/(m1s-ps-m2s))*log(m2s));
  else return(B0A(ps,m1s,m2s)-B0MZ(ps,m1s));
}

static Double fsqrt(Double x)  // = sqrt(1-x)-1
{
  if (abs(x)<.15) {
    Double add=-x/2; Double r=add;
    add *= x/4; r+=add;
    add *= x/2; r+=add;
    add *= 5*x/8; r+=add;
    add *= 7*x/10; r+=add;
    add *= 9*x/12; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 11*x/14; r+=add;
    add *= 13*x/16; r+=add;
    add *= 15*x/18; r+=add;
    add *= 17*x/20; r+=add;
    add *= 19*x/22; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 21*x/24; r+=add;
    add *= 23*x/26; r+=add;
    add *= 25*x/28; r+=add;
    add *= 27*x/30; r+=add;
    add *= 29*x/32; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 31*x/34; r+=add;
    add *= 33*x/36; r+=add;
    add *= 35*x/38; r+=add;
    add *= 37*x/40; r+=add;
    add *= 39*x/42; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 41*x/44; r+=add;
    add *= 43*x/46; r+=add;
    add *= 45*x/48; r+=add;
    add *= 47*x/50; r+=add;
    add *= 49*x/52; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 51*x/54; r+=add;
    add *= 53*x/56; r+=add;
    add *= 55*x/58; r+=add;
    add *= 57*x/60; r+=add;
    add *= 59*x/62; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 61*x/64; r+=add;
    add *= 63*x/66; r+=add;
    add *= 65*x/68; r+=add;
    add *= 67*x/70; r+=add;
    add *= 69*x/72; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 71*x/74; r+=add;
    add *= 73*x/76; r+=add;
    add *= 75*x/78; r+=add;
    add *= 77*x/80; r+=add;
    add *= 79*x/82; r+=add;
    return(r); }
  else
    return(sqrt(1-x)-1);
}
 
static Double log1Mx(Double x) // log(1-x)
{ if (abs(x)<.15) {
    Double add=x; Double r=x;
    add *= x; r+=add/2;
    add *= x; r+=add/3;
    add *= x; r+=add/4;
    add *= x; r+=add/5;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(-r);
    add *= x; r+=add/6;
    add *= x; r+=add/7;
    add *= x; r+=add/8;
    add *= x; r+=add/9;
    add *= x; r+=add/10;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(-r);
    add *= x; r+=add/11;
    add *= x; r+=add/12;
    add *= x; r+=add/13;
    add *= x; r+=add/14;
    add *= x; r+=add/15;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(-r);
    add *= x; r+=add/16;
    add *= x; r+=add/17;
    add *= x; r+=add/18;
    add *= x; r+=add/19;
    add *= x; r+=add/20;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(-r);
    add *= x; r+=add/21;
    add *= x; r+=add/22;
    add *= x; r+=add/23;
    add *= x; r+=add/24;
    add *= x; r+=add/25;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(-r);
    add *= x; r+=add/26;
    add *= x; r+=add/27;
    add *= x; r+=add/28;
    add *= x; r+=add/29;
    add *= x; r+=add/30;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(-r);
    add *= x; r+=add/31;
    add *= x; r+=add/32;
    add *= x; r+=add/33;
    add *= x; r+=add/34;
    add *= x; r+=add/35;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(-r);
    add *= x; r+=add/36;
    add *= x; r+=add/37;
    add *= x; r+=add/38;
    add *= x; r+=add/39;
    add *= x; r+=add/40;
    return(-r); };
  return(log(1-x));
}

static Double Mlog1Mx1(Double x)  // -log(1-x)-x
{
    Double add = x*x; Double r=add/2;
    add *= x; r+=add/3;
    add *= x; r+=add/4;
    add *= x; r+=add/5;
    add *= x; r+=add/6;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/7;
    add *= x; r+=add/8;
    add *= x; r+=add/9;
    add *= x; r+=add/10;
    add *= x; r+=add/11;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/12;
    add *= x; r+=add/13;
    add *= x; r+=add/14;
    add *= x; r+=add/15;
    add *= x; r+=add/16;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/17;
    add *= x; r+=add/18;
    add *= x; r+=add/19;
    add *= x; r+=add/20;
    add *= x; r+=add/21;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/22;
    add *= x; r+=add/23;
    add *= x; r+=add/24;
    add *= x; r+=add/25;
    add *= x; r+=add/26;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/27;
    add *= x; r+=add/28;
    add *= x; r+=add/29;
    add *= x; r+=add/30;
    add *= x; r+=add/31;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/32;
    add *= x; r+=add/33;
    add *= x; r+=add/34;
    add *= x; r+=add/35;
    add *= x; r+=add/36;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/37;
    add *= x; r+=add/38;
    add *= x; r+=add/39;
    add *= x; r+=add/40;
    add *= x; r+=add/41;
    return(r);
}

static Double B0Bpart1(Double x1, Double x2, Double ps, Double m1s, Double m2s, Double la)
        // a product of the two series where the leading constant for
        // m2^2 -> Infinity has been subtracted.
{
  Double s1=fsqrt(x1),s2=Mlog1Mx1(x2);
  return((m2s/ps)*((1+s1)*s2+s1*x2)
         +4*m1s*m2s/(sqr(m2s+la)-sqr(m1s-ps)));
}

Cplx B0B(Double ps, Double m1s, Double m2s) // B0 - Delta - log(u^2) -1 + log(m2^2) 
{
  if (m2s==Double(0)) { Errstream << "B0B called with wrong argument\n";
                      return(Cplx(0)); };
  if ((ps==Double(0))&&(m1s==Double(0))) return(Cplx(0));
  if ((ps==Double(0))&&(m1s!=m2s))
    return(m1s*log(m1s/m2s)/(m2s-m1s));
  if ((ps==Double(0))&&(m1s==m2s))
    return(Cplx(-1));
  if ((m1s==Double(0))&&(m2s==ps))
    return(Cplx(1));
  if (m1s==Double(0))
    return(-((m2s-ps)/ps)*clog1(m2s/(m2s-ps))+1);
  Double s,d,ss,ds,la,r1,r2;
  s=sqrt(m1s)+sqrt(m2s);
  ss=s*s;
  d=sqrt(m1s)-sqrt(m2s);
  ds=d*d;
  if (ps<ds/2)
  {
    Double x1=(2*m2s*(m1s+ps)-sqr(m1s-ps))/sqr(m2s);
    la=m2s*(1+fsqrt(x1));
    Double x2=2*ps/(m2s-m1s+ps+la);
    if ((x1>.2)||(x2>.2))
      return(-(2*m1s/(m2s-ps-m1s+la))*log(m2s/m1s)
             +1+(la/ps)*log1Mx(2*ps/(m2s-m1s+ps+la)));
    else
      return(-(2*m1s/(m2s-ps-m1s+la))*log(m2s/m1s)
             -B0Bpart1(x1,x2,ps,m1s,m2s,la));
  };
  Cplx r=1+log(m2s)-log(m1s*m2s)/2+(m1s-m2s)/ps*log(m2s/m1s)/2;
  if (ps<=ds)
  {
    la=sqrt((ss-ps)*(ds-ps));
    r+=la*log((m1s+m2s-ps+la)/(2*sqrt(m1s*m2s)))/ps;
  }
  else if (ps>=ss)
    {
      la=sqrt((ss-ps)*(ds-ps));
      r-=la*(log((-m1s-m2s+ps+la)/(2*sqrt(m1s*m2s)))-PiI)/ps;
    }
    else
    {
      r1=sqrt(ps-ds);
      r2=sqrt(ss-ps);
      r-=2.*r1*r2*atan(r1/r2)/ps;
    };
  return(r);
}



static Cplx F11(Double x)
{ if (abs(x)>15.) {
    Double add=1/x; Double r=add/2;
    add/=x; r+=add/3;
    add/=x; r+=add/4;
    add/=x; r+=add/5;
    add/=x; r+=add/6;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/7;
    add/=x; r+=add/8;
    add/=x; r+=add/9;
    add/=x; r+=add/10;
    add/=x; r+=add/11;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/12;
    add/=x; r+=add/13;
    add/=x; r+=add/14;
    add/=x; r+=add/15;
    add/=x; r+=add/16;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/17;
    add/=x; r+=add/18;
    add/=x; r+=add/19;
    add/=x; r+=add/20;
    add/=x; r+=add/21;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/22;
    add/=x; r+=add/23;
    add/=x; r+=add/24;
    add/=x; r+=add/25;
    add/=x; r+=add/26;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/27;
    add/=x; r+=add/28;
    add/=x; r+=add/29;
    add/=x; r+=add/30;
    add/=x; r+=add/31;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/32;
    add/=x; r+=add/33;
    add/=x; r+=add/34;
    add/=x; r+=add/35;
    add/=x; r+=add/36;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/37;
    add/=x; r+=add/38;
    add/=x; r+=add/39;
    add/=x; r+=add/40;
    add/=x; r+=add/41;
    return(r); };
 return(-x*clog1(1-1/x)-1);
}

inline Cplx F12(Double x)
{ if (abs(x)>15.) {
    Double add=1/x; Double r=add/2;
    add/=x; r+=add/3;
    add/=x; r+=add/4;
    add/=x; r+=add/5;
    add/=x; r+=add/6;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/7;
    add/=x; r+=add/8;
    add/=x; r+=add/9;
    add/=x; r+=add/10;
    add/=x; r+=add/11;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/12;
    add/=x; r+=add/13;
    add/=x; r+=add/14;
    add/=x; r+=add/15;
    add/=x; r+=add/16;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/17;
    add/=x; r+=add/18;
    add/=x; r+=add/19;
    add/=x; r+=add/20;
    add/=x; r+=add/21;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/22;
    add/=x; r+=add/23;
    add/=x; r+=add/24;
    add/=x; r+=add/25;
    add/=x; r+=add/26;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/27;
    add/=x; r+=add/28;
    add/=x; r+=add/29;
    add/=x; r+=add/30;
    add/=x; r+=add/31;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/33;
    add/=x; r+=add/33;
    add/=x; r+=add/34;
    add/=x; r+=add/35;
    add/=x; r+=add/36;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/37;
    add/=x; r+=add/38;
    add/=x; r+=add/39;
    add/=x; r+=add/40;
    add/=x; r+=add/41;
    return(r); };
 return(-x*clog2(1-1/x)-1);
}

inline Cplx F11(Cplx x)
{
  Double r=real(x);
  if (abs(imag(x)) < CPLX_LIMIT_F*abs(r)) return(F11(r));
  return(-x*log(1-1/x)-1.);
}

inline Cplx F12(Cplx x)
{
  Double r=real(x);
  if (abs(imag(x)) < CPLX_LIMIT_F*abs(r)) return(F12(r));
  return(-x*log(1-1/x)-1.);
}

static Cplx F21(Double x)
{ if (abs(x)>15.) {
    Double add=1/x; Double r=add/3;
    add/=x; r+=add/4;
    add/=x; r+=add/5;
    add/=x; r+=add/6;
    add/=x; r+=add/7;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/8;
    add/=x; r+=add/9;
    add/=x; r+=add/10;
    add/=x; r+=add/11;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/12;
    add/=x; r+=add/13;
    add/=x; r+=add/14;
    add/=x; r+=add/15;
    add/=x; r+=add/16;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/17;
    add/=x; r+=add/18;
    add/=x; r+=add/19;
    add/=x; r+=add/20;
    add/=x; r+=add/21;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/22;
    add/=x; r+=add/23;
    add/=x; r+=add/24;
    add/=x; r+=add/25;
    add/=x; r+=add/26;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/27;
    add/=x; r+=add/28;
    add/=x; r+=add/29;
    add/=x; r+=add/30;
    add/=x; r+=add/31;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/32;
    add/=x; r+=add/33;
    add/=x; r+=add/34;
    add/=x; r+=add/35;
    add/=x; r+=add/36;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/37;
    add/=x; r+=add/38;
    add/=x; r+=add/39;
    add/=x; r+=add/40;
    add/=x; r+=add/41;
    return(r); };
 return(-sqr(x)*clog1(1-1/x)-x-.5);
}

static Cplx F22(Double x)
{ if (abs(x)>10.) {
    Double add=1/x; Double r=add/3;
    add/=x; r+=add/4;
    add/=x; r+=add/5;
    add/=x; r+=add/6;
    add/=x; r+=add/7;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/8;
    add/=x; r+=add/9;
    add/=x; r+=add/10;
    add/=x; r+=add/11;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/12;
    add/=x; r+=add/13;
    add/=x; r+=add/14;
    add/=x; r+=add/15;
    add/=x; r+=add/16;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/17;
    add/=x; r+=add/18;
    add/=x; r+=add/19;
    add/=x; r+=add/20;
    add/=x; r+=add/21;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/22;
    add/=x; r+=add/23;
    add/=x; r+=add/24;
    add/=x; r+=add/25;
    add/=x; r+=add/26;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/27;
    add/=x; r+=add/28;
    add/=x; r+=add/29;
    add/=x; r+=add/30;
    add/=x; r+=add/31;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/32;
    add/=x; r+=add/33;
    add/=x; r+=add/34;
    add/=x; r+=add/35;
    add/=x; r+=add/36;
    if (abs(add*SERIES_LIMIT)<abs(r)) return(r);
    add/=x; r+=add/37;
    add/=x; r+=add/38;
    add/=x; r+=add/39;
    add/=x; r+=add/40;
    add/=x; r+=add/41;
    return(r); };
 return(-sqr(x)*clog2(1-1/x)-x-.5);
}

// the derivative of B0 with respect to the momentum
Cplx DB0(Double ps, Double m1s, Double m2s)
{
  if (m1s==Double(0)) interchange(&m1s,&m2s);
  if ((ps==Double(0))&&(m2s==Double(0))&&(m1s !=Double(0) )) {
    return( .5/m1s ); };
  if ((ps==Double(0))&&(m1s != m2s)) {
    Double x=m1s-m2s;
    return( .5/(x*x*x) * (m1s*m1s-m2s*m2s+2.*m1s*m2s*log(m2s/m1s))); };
  if ((ps==Double(0))&&(m1s == m2s)) return(1./(6.*m1s));
  if ((m2s==Double(0))&&(m1s==Double(0))) return(-1./ps);
  if ((m2s==Double(0))&&(m1s != Double(0))) {
    if (abs(ps*10)>m1s) return((-1.+(m1s/ps)*clog1(m1s/(m1s-ps)))/ps);
    else return(Mlog1Mx1(ps/m1s)*m1s/sqr(ps)); };

  Double la1=la(ps,m1s,m2s);

  if (abs(la1) < ASYMP_LIMIT_DB01 *(ps+m1s+m2s) ) {
    if (ps>m1s+m2s) {
      Errstream << "Error: DB0fin diverges on threshold:\n";
      Errstream << "  ps = " << ps << " , m1s = " << m1s << " , m2s = "
	        << m2s << "\n";
      return(Cplx(0));
    };
    Double x=(ps+m2s-m1s)/(2*ps);
    return((-2+(1-2*x)*log((x-1)/x))/ps);
  };

  if (la1>0) {
    Double sla=sqrt(la1);
    Double x1=(ps+m2s-m1s+sla)/(2*ps);
    Double x2=(ps+m2s-m1s-sla)/(2*ps);
    return(-1/sla * ((1-x1)*F11(x1)-(1-x2)*F12(x2))); };

  Cplx rps=.5*((m1s-m2s)/ps)*log(m1s/m2s);
  Double s,d,ss,ds,la,r1,r2,sm=2.*sqrt(m1s*m2s);
  Cplx rr,lo;
  s=sqrt(m1s)+sqrt(m2s);
  ss=s*s;
  d=sqrt(m1s)-sqrt(m2s);
  ds=d*d;
  r1=sqrt(ps-ds);
  r2=sqrt(ss-ps);
  rr=(m1s+m2s-ps-I*r1*r2)/sm;
  lo=-I*2.*atan(r1/r2);
  rps+=(.5*sm/ps)*(1/rr-rr)*lo-1.-((rr*rr+1)/(rr*rr-1))*lo;
  return(rps/ps);
}


// the second derivative of B0 with respect to the momentum
Cplx DDB0(Double ps, Double m1s, Double m2s)
{

  Double la1=la(ps,m1s,m2s);

//  if (abs(la1) < ASYMP_LIMIT_DB01 *(ps+m1s+m2s) ) {
//    if (ps>m1s+m2s) {
//      Errstream << "Error: DDB0fin diverges on threshold:\n";
//      Errstream << "  ps = " << ps << " , m1s = " << m1s << " , m2s = "
//	        << m2s << "\n";
//      return(Cplx(0));
//    };
//    Double x=(ps+m2s-m1s)/(2*ps);
//    return((-2+(1-2*x)*log((x-1)/x))/ps);
//  };

  // TO IMPROVE: very unstable formula:
  if (la1>0) {
    Double sla=sqrt(la1);
    Double x1=(ps+m2s-m1s+sla)/(2*ps);
    Double x2=(ps+m2s-m1s-sla)/(2*ps);
    return(1/la1 * (-(ps-m1s-m2s)*DB0(ps,m1s,m2s)
		   -2+x1+x2
		   +(1-x1)*(2*x1-1)*F11(x1)+(1-x2)*(2*x2-1)*F12(x2))); }
  else {
    Cplx sla=csqrt1(la1);
    Cplx x1=(ps+m2s-m1s+sla)/(2*ps);
    Cplx x2=(ps+m2s-m1s-sla)/(2*ps);
    return(1/la1 * (-(ps-m1s-m2s)*DB0(ps,m1s,m2s)
		   -2+x1+x2
		   +(1-x1)*(2*x1-1)*F11(x1)+(1-x2)*(2*x2-1)*F12(x2))); }
}



static Double fsqrt1(Double x)  // = sqrt(1-x)-1+x/2
{
  if (abs(x)<.15) {
    Double add=-x/2; Double r=0;
    add *= x/4; r+=add;
    add *= x/2; r+=add;
    add *= 5*x/8; r+=add;
    add *= 7*x/10; r+=add;
    add *= 9*x/12; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 11*x/14; r+=add;
    add *= 13*x/16; r+=add;
    add *= 15*x/18; r+=add;
    add *= 17*x/20; r+=add;
    add *= 19*x/22; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 21*x/24; r+=add;
    add *= 23*x/26; r+=add;
    add *= 25*x/28; r+=add;
    add *= 27*x/30; r+=add;
    add *= 29*x/32; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 31*x/34; r+=add;
    add *= 33*x/36; r+=add;
    add *= 35*x/38; r+=add;
    add *= 37*x/40; r+=add;
    add *= 39*x/42; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 41*x/44; r+=add;
    add *= 43*x/46; r+=add;
    add *= 45*x/48; r+=add;
    add *= 47*x/50; r+=add;
    add *= 49*x/52; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 51*x/54; r+=add;
    add *= 53*x/56; r+=add;
    add *= 55*x/58; r+=add;
    add *= 57*x/60; r+=add;
    add *= 59*x/62; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 61*x/64; r+=add;
    add *= 63*x/66; r+=add;
    add *= 65*x/68; r+=add;
    add *= 67*x/70; r+=add;
    add *= 69*x/72; r+=add;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= 71*x/74; r+=add;
    add *= 73*x/76; r+=add;
    add *= 75*x/78; r+=add;
    add *= 77*x/80; r+=add;
    add *= 79*x/82; r+=add;
    return(r); }
  else
    return(sqrt(1-x)-1+x/2);
}



// the derivative of B0 with respect to the momentum,
// with 1/(2 m2^2) subtracted
Cplx DB0D(Double ps, Double m1s, Double m2s)
{
  if (m2s==Double(0))  { Errstream << "B0D called with wrong argument\n";
                 return(Cplx(0)); };
  if ((ps==Double(0))&&(m1s==Double(0)))
    return(Cplx(0));
  if ((ps==Double(0))&&(m1s == m2s)) return(-1/(3.*m1s));
  if ((ps==Double(0))&&(m1s != m2s)) {
    if (abs((m1s-m2s)/(m1s+m2s))<ASYMP_LIMIT_DB01)
      return(-2/(3*(m1s+m2s)));
    Double x=m1s-m2s;
    return( .5/(x*x*x*m2s) * m1s 
            * (4*m1s*m2s-sqr(m1s)-3*sqr(m2s)+2*sqr(m2s)*log(m2s/m1s) )); };
  if (m1s==Double(0)) {
    if (abs(ps)*10<abs(m2s))
      return(F21(m2s/ps)/m2s);
    else
      return((-1+(m2s/ps)*clog1(m2s/(m2s-ps))-ps/(2*m2s))/ps); };

  Double la1=la(ps,m1s,m2s);
  if (la1>0) {
    if ((abs(ps)+abs(m1s))*10<m2s) {
      Double xi=(2*m2s*(m1s+ps)-sqr(m1s-ps))/sqr(m2s);
      Double fsqrtx=fsqrt(xi);
      Double sla=m2s*(1+fsqrtx);
      Double x1=(ps+m2s-m1s+sla)/(2*ps);
      Double x2=(ps-m1s-m2s*fsqrtx)/(2*ps);
      Double x2M1=(-sqr(ps-m1s)/(2*m2s)-m2s*fsqrt1(xi))/(2*ps);
      return(1/sla
              *(x2M1+x2*x2M1*log(x2M1/x2)-F11(x1)+F21(x1)
                -fsqrtx/2)); }
    else {
      Double sla=sqrt(la1);
      Double x1=(ps+m2s-m1s+sla)/(2*ps);
      Double x2=(ps+m2s-m1s-sla)/(2*ps);
      return(-1/(2*m2s)-1/sla * ((1-x1)*F11(x1)-(1-x2)*F12(x2))); }; };

  Cplx rps=-ps/(2*m2s)+.5*((m1s-m2s)/ps)*log(m1s/m2s);
  Double s,d,ss,ds,sla,r1,r2,sm=2*sqrt(m1s*m2s);
  Cplx rr,lo;
  s=sqrt(m1s)+sqrt(m2s);
  ss=s*s;
  d=sqrt(m1s)-sqrt(m2s);
  ds=d*d;
  r1=sqrt(ps-ds);
  r2=sqrt(ss-ps);
  rr=(m1s+m2s-ps-I*r1*r2)/sm;
  lo=-I*2*atan(r1/r2);
  rps+=(.5*sm/ps)*(1/rr-rr)*lo-1-((rr*rr+1)/(rr*rr-1))*lo;
  return(rps/ps);
}

// the derivative of B0 with respect to the first (squared) mass
Cplx DM1B0(Double ps, Double m1s, Double m2s)
{
  if (m1s==Double(0)) {
    Errstream << "DM1B0 is infrared singular for m_1^2 = 0\n";
    return(Cplx(0)); };
  if (m2s==Double(0)) {
    if (ps==Double(0)) return(-1/m1s);
    else return(clog2((m1s-ps)/m1s)/ps); };
  if ((ps==Double(0))&&(m1s!=m2s))
    return( m2s/sqr(m1s-m2s) * log(m1s/m2s) + 1/(m2s-m1s));
  if (ps==Double(0))
    { if (m1s==m2s) return(-.5/m1s);
      else
        if (m2s==0) return(-1/m1s);
        else return((-m1s + m2s + m2s*log(m1s) 
                     - m2s*log(m2s))/sqr(-m1s + m2s)); };
  Double la1=la(ps,m1s,m2s);
  if (la1<0.) {
    Cplx sla=csqrt2(la1);
    Cplx x1=(ps+m2s-m1s+sla)/(2*ps);
    Cplx x2=(ps+m2s-m1s-sla)/(2*ps);
    return((F11(x1)-F12(x2))/sla); }
  else {
    Double sla=sqrt(la1);
    Double x1=(ps+m2s-m1s+sla)/(2*ps);
    Double x2=(ps+m2s-m1s-sla)/(2*ps);
    return((F11(x1)-F12(x2))/sla); };
}

// the derivative of B0 with respect to the first (squared) mass
// and the momentum
Cplx DDM1B0(Double ps, Double m1s, Double m2s)
{
  if (m1s==Double(0)) {
    Errstream << "DDM1B0 is infra-red singular for m_1^2 = 0\n";
    return(Cplx(0)); };
  if ((m2s==Double(0))&&(ps!=Double(0))) {
    if (ps==m1s) {
      Errstream << "DDM1B0 not implemented for m_2^2 = 0 , p^2 = m_1^2\n";
      return(Cplx(0)); };
    return((clog1(m1s/(m1s-ps))+ps/(ps-m1s))/sqr(ps)); };
  if (ps==Double(0))
    { if (m2s==0) return(-.5/sqr(m1s));
      else 
        if (m2s==m1s) return(-1/12/sqr(m1s));
        else return((-4*m1s*m2s + 4*m1s*m2s*log(m1s) - 4*m1s*m2s*log(m2s)
                     - sqr(m1s) + 
                     5*sqr(m2s) + 2*log(m1s)*sqr(m2s) - 2*log(m2s)*sqr(m2s))/
		    (2*powint(-m1s + m2s,4))); };
  Double la1=la(ps,m1s,m2s);
  if (la1<0.) {
    Cplx sla=csqrt2(la1);
    Cplx x1=(ps+m2s-m1s+sla)/(2*ps);
    Cplx x2=(ps+m2s-m1s-sla)/(2*ps);
    return((2+(1-x1)*F11(x1)+(1-x2)*F12(x2)+(m2s+m1s-ps)*(F11(x1)-F12(x2))/sla)
	   /la1); }
  else {
    Double sla=sqrt(la1);
    Double x1=(ps+m2s-m1s+sla)/(2*ps);
    Double x2=(ps+m2s-m1s-sla)/(2*ps);
    return((2+(1-x1)*F11(x1)+(1-x2)*F12(x2)+(m2s+m1s-ps)*(F11(x1)-F12(x2))/sla)
	   /la1); };
}

// B0C = B0 - Delta - log(u^2) -1 + log(m2^2) - (m1^2/m2^2) log(m1^2/m2^2)
Cplx B0C(Double ps, Double m1s, Double m2s)
{
  if ((m2s==Double(0))||(m1s==Double(0))) {
    Errstream << "B0C called with wrong argument\n";
    return(Cplx(0)); };
  if ((ps==Double(0))&&(m1s!=m2s))
    return(sqr(m1s)*log(m1s/m2s)/(m2s-m1s)/m2s);
  if ((ps==Double(0))&&(m1s==m2s))
    return(Cplx(-1));
  Double s,d,ss,ds,sla,r1,r2;
  s=sqrt(m1s)+sqrt(m2s);
  ss=s*s;
  d=sqrt(m1s)-sqrt(m2s);
  ds=d*d;
  if (abs(ps)<ds/2)
  {
    Double x1=(2*m2s*(m1s+ps)-sqr(m1s-ps))/sqr(m2s);
    sla=m2s*(1+fsqrt(x1));
    Double x2=2*ps/(m2s-m1s+ps+sla);
    if ((abs(x1)>.2)||(abs(x2)>.2))
      return(4 *(m1s/m2s) *(m2s*(m1s+ps)+m1s*ps) * log(m1s/m2s)
               /((m2s-ps-m1s+sla)*(m2s+ps+m1s+sla))
             +1+(sla/ps)*log1Mx(2*ps/(m2s-m1s+ps+sla)));
    else
      return(4 *m1s/m2s *(m2s*(m1s+ps)+m1s*ps)
               /(m2s-ps-m1s+sla) /(m2s+ps+m1s+sla) * log(m1s/m2s)
             -B0Bpart1(x1,x2,ps,m1s,m2s,sla));
  };
  Cplx r=1+log(m2s)-log(m1s*m2s)/2+(m1s-m2s)/ps*log(m2s/m1s)/2
          - (m1s/m2s)* log(m1s/m2s);
  if (ps<=ds)
  {
    sla=sqrt((ss-ps)*(ds-ps));
    r+=sla*log((m1s+m2s-ps+sla)/(2*sqrt(m1s*m2s)))/ps;
  }
  else if (ps>=ss)
    {
      sla=sqrt((ss-ps)*(ds-ps));
      r-=sla*(log((-m1s-m2s+ps+sla)/(2*sqrt(m1s*m2s)))-PiI)/ps;
    }
    else
    {
      r1=sqrt(ps-ds);
      r2=sqrt(ss-ps);
      r-=2*r1*r2*atan(r1/r2)/ps;
    };
  return(r);
}

static Double Mlog1Mx2(Double x)  // -log(1-x)-x-x^2/2
{
    Double add = x*x; Double r=0;
    add *= x; r+=add/3;
    add *= x; r+=add/4;
    add *= x; r+=add/5;
    add *= x; r+=add/6;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/7;
    add *= x; r+=add/8;
    add *= x; r+=add/9;
    add *= x; r+=add/10;
    add *= x; r+=add/11;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/12;
    add *= x; r+=add/13;
    add *= x; r+=add/14;
    add *= x; r+=add/15;
    add *= x; r+=add/16;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/17;
    add *= x; r+=add/18;
    add *= x; r+=add/19;
    add *= x; r+=add/20;
    add *= x; r+=add/21;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/22;
    add *= x; r+=add/23;
    add *= x; r+=add/24;
    add *= x; r+=add/25;
    add *= x; r+=add/26;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/27;
    add *= x; r+=add/28;
    add *= x; r+=add/29;
    add *= x; r+=add/30;
    add *= x; r+=add/31;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/32;
    add *= x; r+=add/33;
    add *= x; r+=add/34;
    add *= x; r+=add/35;
    add *= x; r+=add/36;
    if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
    add *= x; r+=add/37;
    add *= x; r+=add/38;
    add *= x; r+=add/39;
    add *= x; r+=add/40;
    add *= x; r+=add/41;
    return(r);
}


static Double B0Dpart1(Double x1, Double x2, Double ps, Double m1s, Double m2s, Double sla)
        // like B0Bpart1, but additional subtraction p^2/(2 m2^2)
{
  Double s1=fsqrt(x1),s2=Mlog1Mx1(x2),s3=fsqrt1(x1),s4=Mlog1Mx2(x2);
  return((m2s/ps)*(s1*s2+s3*x2+s4)
   + (-3*sla*m1s*ps - powint(m1s,3) + 2*powint(ps,3) + sla*sqr(m1s) + 
       m2s*sqr(m1s) + 4*ps*sqr(m1s) + 2*sla*sqr(ps) - 
       5*m1s*sqr(ps) - 3*m2s*sqr(ps))
      /(m2s*sqr(sla - m1s + m2s + ps)) 
   - 3*m2s*ps*s1 / sqr(sla - m1s + m2s + ps) );
}




// B0D = B0 - Delta - log(mu^2) -1 + log(m2^2) - (m1^2/m2^2) log(m1^2/m2^2)
//        - ps/(2*m2s)
Cplx B0D(Double ps, Double m1s, Double m2s)
{
  if (m2s==Double(0)) { Errstream << "B0D called with wrong argument\n";
                return(Cplx(0)); };
  if (m1s==Double(0)) {
    if (ps==Double(0)) return(Cplx(0));
// following line added 09 Nov 99:
    if (ps==m2s) return(Cplx(0.5));
    if (abs(ps*10)>m2s) return((1-m2s/ps)*clog1(m2s/(m2s-ps))+1-ps/(2*m2s));
    else {
      Double x=ps/m2s;
      Double add = x*x; Double r=add/6;
      add *= x; r+=add/12;
      add *= x; r+=add/20;
      add *= x; r+=add/30;
      add *= x; r+=add/42;
      if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
      add *= x; r+=add/56;
      add *= x; r+=add/72;
      add *= x; r+=add/90;
      add *= x; r+=add/110;
      add *= x; r+=add/132;
      if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
      add *= x; r+=add/156;
      add *= x; r+=add/182;
      add *= x; r+=add/210;
      add *= x; r+=add/240;
      add *= x; r+=add/272;
      if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
      add *= x; r+=add/306;
      add *= x; r+=add/342;
      add *= x; r+=add/380;
      add *= x; r+=add/420;
      add *= x; r+=add/462;
      if (abs(add)*SERIES_LIMIT<abs(r)) return(r);
      add *= x; r+=add/506;
      add *= x; r+=add/552;
      add *= x; r+=add/600;
      add *= x; r+=add/650;
      add *= x; r+=add/702;
      return(r); }; };
			       
  if ((ps==Double(0))&&(m1s!=m2s))
    return(sqr(m1s)*log(m1s/m2s)/(m2s-m1s)/m2s);
  if ((ps==Double(0))&&(m1s==m2s))
    return(Cplx(-1));
  Double s,d,ss,ds,sla,r1,r2;
  s=sqrt(m1s)+sqrt(m2s);
  ss=s*s;
  d=sqrt(m1s)-sqrt(m2s);
  ds=d*d;
  if (abs(ps)<ds/2)
  {
    Double x1=(2*m2s*(m1s+ps)-sqr(m1s-ps))/sqr(m2s);
    sla=m2s*(1+fsqrt(x1));
    Double x2=2*ps/(m2s-m1s+ps+sla);
    if ((abs(x1)>.2)||(abs(x2)>.2))
      return(4 *(m1s/m2s) *(m2s*(m1s+ps)+m1s*ps) * log(m1s/m2s)
               /((m2s-ps-m1s+sla)*(m2s+ps+m1s+sla))-ps/(2*m2s)
             +1+(sla/ps)*log1Mx(2*ps/(m2s-m1s+ps+sla)));
    else
      return(4 *m1s/m2s *(m2s*(m1s+ps)+m1s*ps)
               /(m2s-ps-m1s+sla) /(m2s+ps+m1s+sla) * log(m1s/m2s)
             -B0Dpart1(x1,x2,ps,m1s,m2s,sla));
  };
  Cplx r=1+log(m2s)-log(m1s*m2s)/2+(m1s-m2s)/ps*log(m2s/m1s)/2
          - (m1s/m2s)* log(m1s/m2s)-ps/(2*m2s);
  if (ps<=ds)
  {
    sla=sqrt((ss-ps)*(ds-ps));
    r+=sla*log((m1s+m2s-ps+sla)/(2*sqrt(m1s*m2s)))/ps;
  }
  else if (ps>=ss)
    {
      sla=sqrt((ss-ps)*(ds-ps));
      r-=sla*(log((-m1s-m2s+ps+sla)/(2*sqrt(m1s*m2s)))-PiI)/ps;
    }
    else
    {
      r1=sqrt(ps-ds);
      r2=sqrt(ss-ps);
      r-=2*r1*r2*atan(r1/r2)/ps;
    };
  return(r);
}

