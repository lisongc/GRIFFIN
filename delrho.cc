/*-----------------------------------------------------------------------------
delrho.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 9 Nov 2021
-------------------------------------------------------------------------------
higher-order corrections to the rho parameter
-----------------------------------------------------------------------------*/

#include "li.h"
#include "delrho.h"

// from CERN 95-03, p.187, eq.(45) 
double getalphas(double musqval, double mzsqval, double alphasmz)
{
  static double b0 = 1.9166666666666667, b1 = 2.4166666666666665,
  		b2 = 2.826678240740741;
  double l,lold=0,ll;
  for(l = Pi/b0/alphasmz; abs(l-lold) > 1e-12*l; )
  {
    lold = l;
    ll = log(l);
    l = Pi/b0/alphasmz * (1. - 1./b0/l*b1*ll/b0
  	  + 1./sqr(b0*l)*(sqr(b1/b0)*(ll*ll-ll-1)+b2/b0));
  }
  l += log(musqval/mzsqval);
  ll = log(l);
  return(Pi/b0/l * (1. - 1./b0/l*b1*ll/b0
  	  + 1./sqr(b0*l)*(sqr(b1/b0)*(ll*ll-ll-1)+b2/b0)));
}

inline double Cl2(double x)
{
  return(imag(li2(exp(I*x))));
}

double phi(double z)
{
  if(z <= 1)
    return(4*sqrt(z/(1-z)) * Cl2(2*asin(sqrt(z))));
  else
  {
    double u = (1-sqrt(1-1/z))/2.;
    return(1./sqrt(1-1/z) * (-4*real(li2(u)) + 2*sqr(log(u)) - sqr(log(4*z)) 
           + 2*PiSonSix));
  }
}

double g(double x)
{
  if(x <= 4)
    return(sqrt(4-x) * (Pi - 2*asin(sqrt(x/4))));
  else
    return(2*sqrt(x/4-1) * log((1-sqrt(1-4/x))/(1+sqrt(1-4/x))));
}

double Rmt4(double ht)
{
  return( 25 - 4*ht + (0.5 - 1/ht)*Pi*Pi + (ht-4) * sqrt(ht)*g(ht)/2. +
         (-6 - 6*ht + sqr(ht)/2.) * log(ht) +
         (6/ht - 15 + 12*ht - 3*sqr(ht)) * real(li2(1-ht)) +
         1.5*(-10 + 6*ht - sqr(ht))*phi(ht/4.) );
}

double delrho2a2(const inval* ival)   // yt^4 corrections to \Delta\rho
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MHs = sqr(ival->get(MH)),
         MTs = sqr(ival->get(MT)),
         AL = ival->get(al) * DRHOSCHEME,
         ALS = ival->get(als);
  return(3*sqr(AL/(1-MWs/MZs) * MTs/(16*Pi*MWs)) * Rmt4(MHs/MTs));
}

double delrho2aas(const inval* ival)  // yt^2*as corrections to \Delta\rho
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         AL = ival->get(al) * DRHOSCHEME,
         ALS = ival->get(als);
  return(-3*AL/(1-MWs/MZs) * MTs/(16*Pi*MWs) *
          getalphas(MTs,MZs,ALS)/Pi * 4/3. *(.5 + PiSonSix));
}

double delrho3a2as(const inval* ival) // yt^4*as corrections to \Delta\rho
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
	 r = ival->get(MH)/ival->get(MT),
         AL = ival->get(al) * DRHOSCHEME,
         ALS = ival->get(als);
// use asymptotic formulas from hep-ph/0302275 and Pade approximant for 
// intermediate region
  double lr = log(sqr(r)), rr;
  if(r < 1)
    rr = 157.295 + 112*(-1 + r) - 3.52*powint(-1 + r,4) + 2.06*powint(-1 + r,5) - 24.73*sqr(-1 + r) + 7.39*trip(-1 + r);
  else {
    if(r > 6.32456)
      rr = 79.73 - 47.77*log(4/sqr(r)) + 42.07*sqr(log(4/sqr(r))) + (4*(225.16 - 179.74*log(4/sqr(r)) + 70.22*sqr(log(4/sqr(r))) - 19.22*trip(log(4/sqr(r)))))/sqr(r) + 16*powint(r,-4)*(-76.07 + 25.33*log(4/sqr(r)) - 9.17*sqr(log(4/sqr(r))) - 5.57*trip(log(4/sqr(r)))) + 64*powint(r,-6)*(-10.1 - 24.69*log(4/sqr(r)) - 0.3*sqr(log(4/sqr(r))) - 5.46*trip(log(4/sqr(r)))) + 256*powint(r,-8)*(-4.52 - 32.85*log(4/sqr(r)) + 0.72*sqr(log(4/sqr(r))) - 5.25*trip(log(4/sqr(r)))) + 1024*powint(r,-10)*(-2.55 - 36.61*log(4/sqr(r)) + 1.06*sqr(log(4/sqr(r))) - 5.14*trip(log(4/sqr(r)))) + 9.*trip(log(4/sqr(r)));
    else
      rr = (1.6130162464188158e6 + 0.24095899337779886*powint(-40. + sqr(r),4) - 0.000276067198925092*powint(-40. + sqr(r),5) + 135496.63717141392*(-40. + sqr(r)) + 4172.995142333796*sqr(-40. + sqr(r)) + 54.922842086270755*trip(-40. + sqr(r)))/(1. + 0.0005356504084810248*powint(-1. + sqr(r),4) + 2.5454367938286495e-8*powint(-1. + sqr(r),5) + 1.2259290669021898*(-1. + sqr(r)) + 0.3897464368845946*sqr(-1. + sqr(r)) + 0.03201481261896207*trip(-1. + sqr(r)));
  }
  return(sqr(AL/(1-MWs/MZs) * MTs/(16*Pi*MWs)) * getalphas(MTs,MZs,ALS)/Pi * rr);
}

double delrho3a3(const inval* ival)   // yt^6 corrections to \Delta\rho
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
	 r = ival->get(MH)/ival->get(MT),
         AL = ival->get(al) * DRHOSCHEME;
  double lr = log(sqr(r)), rr;
// use asymptotic formulas from hep-ph/0302275 and Pade approximant for 
// intermediate region
  if(r < 1)
    rr = 95.92 - 111.98*(-1 + r) + 7.27*powint(-1 + r,4) - 15.6*powint(-1 + r,5) + 8.099*sqr(-1 + r) + 9.36*trip(-1 + r);
  else {
    if(r > 6.32456)
      rr = -189.93 - 231.48*log(4/sqr(r)) + ((-3.17 - 83.25*log(4/sqr(r)))*sqr(r))/4. - 142.06*sqr(log(4/sqr(r))) + 2.75*trip(log(4/sqr(r))) + 16*powint(r,-4)*(227.55 - 510.55*log(4/sqr(r)) + 87.77*sqr(log(4/sqr(r))) + 6.41*trip(log(4/sqr(r)))) + 64*powint(r,-6)*(-58.4 - 329.18*log(4/sqr(r)) + 20.42*sqr(log(4/sqr(r))) + 14.54*trip(log(4/sqr(r)))) + 256*powint(r,-8)*(-36.14 - 381.88*log(4/sqr(r)) + 18.63*sqr(log(4/sqr(r))) + 15.04*trip(log(4/sqr(r)))) + 1024*powint(r,-10)*(-39.08 - 416.36*log(4/sqr(r)) + 13.76*sqr(log(4/sqr(r))) + 17.19*trip(log(4/sqr(r)))) + (4*(-332.34 + 77.71*log(4/sqr(r)) - 68.67*sqr(log(4/sqr(r))) + 51.79*trip(log(4/sqr(r)))))/sqr(r);
    else
      rr = (867128.685562442 + 1.3139429394527067*powint(-40. + sqr(r),4) + 0.006216850000329791*powint(-40. + sqr(r),5) + 102284.64023019791*(-40. + sqr(r)) + 4796.812370243386*sqr(-40. + sqr(r)) + 112.15129331129012*trip(-40. + sqr(r)))/(1. + 0.00007543456213005985*powint(-1. + sqr(r),4) - 6.719145520243815e-8*powint(-1. + sqr(r),5) + 0.8347606037834723*(-1. + sqr(r)) + 0.09484010402472164*sqr(-1. + sqr(r)) + 0.00532696624331143*trip(-1. + sqr(r)));
  }
  return(trip(AL/(1-MWs/MZs) * MTs/(16*Pi*MWs)) * rr);
}

double delrho3aas2(const inval* ival) // yt^2*as^2 corrections to \Delta\rho
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         AL = ival->get(al) * DRHOSCHEME,
         ALS = ival->get(als);
  return(3*AL/(1-MWs/MZs) * MTs/(16*Pi*MWs) *
          sqr(getalphas(MTs,MZs,ALS)/Pi) * (-14.5940));
}

double delrho4aas3(const inval* ival) // yt^2*as^3 corrections to \Delta\rho
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         AL = ival->get(al) * DRHOSCHEME,
         ALS = ival->get(als);
  return(3*AL/(1-MWs/MZs) * MTs/(16*Pi*MWs) *
          trip(getalphas(MTs,MZs,ALS)/Pi) * (-93.15));
}
