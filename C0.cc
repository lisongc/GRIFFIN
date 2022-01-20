/* ------------------------------------------------------------- */
/* C0.cc                                                         */
/* Ayres Freitas (afreitas@pitt.edu)			         */
/* last revision 06.08.19                                        */
/* ------------------------------------------------------------- */
/* scalar one-loop function C0                     	         */
/* ------------------------------------------------------------- */

#include "oneloop.h"
#include "li.h"

#define Errstream cerr


#ifdef QUADPRECISION
static doubledouble ZERO_LIMIT = SDouble("1e-28");
static doubledouble ZEROS_LIMIT = SDouble("1e-15");
#else
#define ZERO_LIMIT 1e-14
#define ZEROS_LIMIT 1e-8
#endif


Double sign0(Double var)
{
  return ((var >= 0.0) ? 1.0 : -1.0);
}

Cplx Eta(Cplx c1, Cplx c2)
{
  Double im1 = imag(c1), im2 = imag(c2), im12 = imag(c1*c2);
  
  if((im1 == Double(0) && real(c1) < -ZERO_LIMIT/*Double(0)*/) ||
     (im2 == Double(0) && real(c2) < -ZERO_LIMIT/*Double(0)*/) ||
     (im12== Double(0) && real(c1*c2) < -ZERO_LIMIT/*Double(0)*/))
  {
//    Errstream << "Eta function on cut;  c1 = " << c1 << "  c2 = " << c2
//	    	<< "  c1*c2 = " << c1*c2 << endl;
//    return(Cplx(HUGE_VAL));
    return(Cplx(0));
  }
  if(im1 < Double(0) && im2 < Double(0) && im12 > Double(0))
    return(TwoPiI);
  if(im1 > Double(0) && im2 > Double(0) && im12 < Double(0))
    return(-TwoPiI);
  return Cplx(0);
}

#ifdef QUADPRECISION
static const Cplx IEPS(Double(0), SDouble("1e-30"));
#else
static const Cplx IEPS(0, 1e-15);
#endif

// the one-loop scalar triangle function
//   implementation according to A. Denner, Fortschr. Phys. 41 (1993) 307
Cplx C0(Double p10, Double p12, Double p20, Double m0, Double m1, Double m2)
{
  Cplx c0, al, al0, al1, al2, y0[3], x[3][2], y[3][2];
  Double m, p, thp[3], thy[3];
  int i, j;
// Special case of (al least) one zero momentum

  if(abs(p10) <= ZERO_LIMIT || abs(p12) <= ZERO_LIMIT || abs(p20) <= ZERO_LIMIT)
  {
    // Permutate until p10==0
    while(abs(p10) > ZERO_LIMIT)
    {
      m = m0;  m0 = m1;  m1 = m2;  m2 = m;
      p = p10; p10= p12; p12= p20; p20= p;
    }
    if(abs(p12-p20) <= ZERO_LIMIT)
    {
// Regular C0 function with p10==0, p12 == p20
      if(abs(m0-m1) <= ZERO_LIMIT)
	return DM1B0(p12,m0,m2);
      else
	return (B0A(p12,m0,m2) - B0A(p12,m1,m2))/(m0-m1);
    }
    
    if(m0*m1 == Double(0))
    {
      // Permutate until m0==0
      if(m1 == Double(0))
      {
	m1 = m0;  m0 = Double(0);
	p = p12;  p12 = p20;  p20 = p;
      }      	
      if(m1 == Double(0))
      {
	Errstream << "IR divergency for p1=0 and m0,m1=0 cannot be handled"
		<< endl;
	return(Cplx(0));
      }
// Regular C0 function with p10==0, p12 != p20, m0==0
      al = 1. + (p12-p20)/(-m1 - IEPS*m1);
//      al1 = sqrt( la(p12,m1,m2) + IEPS*p12*(p12-m1-m2) );
      al1 = sqrt( la(p12,m1,m2) * (1.+ IEPS*sign0(p12)) );
      if(al1 == Cplx(0)) al1 = IEPS*abs(al);
      x[0][0] = m2/(p20-m2+IEPS*m1);
      x[1][0] = (p12-m1-m2+al1)/(2*m1);
      x[1][1] = (p12-m1-m2-al1)/(2*m1);
      for(j = 0; j < 2; j++)
	c0 += li2(1.+x[1][j]) - li2(1.+x[1][j]/al)
		- Eta(1./al,-x[1][j]) * log(1.+x[1][j]/al);
      c0 += li2(1.+x[0][0]/al) - li2(1.+x[0][0])
	      + Eta(1./al,-x[0][0]) * log(1.+x[0][0]/al);
      c0 += log(al) * log((m2-p20)/m1 - IEPS) - sqr(log(al))/2.;
      c0 /= (p20-p12);
      return c0;
    }
    else
    {
// Regular C0 function with p10==0, p12 != p20, m0,m1 != 0
      al = 1. + (p12-p20)/(m0-m1 - IEPS*(m0+m1));
//      al0 = sqrt( la(p20,m2,m0) + IEPS*p20*(p20-m0-m2) );
//      al1 = sqrt( la(p12,m1,m2) + IEPS*p12*(p12-m1-m2) );
      al0 = sqrt( la(p20,m2,m0) * (1.+ IEPS*sign0(p20)) );
      al1 = sqrt( la(p12,m1,m2) * (1.+ IEPS*sign0(p12)) );
      if(al0 == Cplx(0)) al0 = IEPS*abs(p20);
      if(al1 == Cplx(0)) al1 = IEPS*abs(p12);
      x[0][0] = (p20-m0-m2+al0)/(2*m0);
      x[0][1] = (p20-m0-m2-al0)/(2*m0);
      x[1][0] = (p12-m1-m2+al1)/(2*m1);
      x[1][1] = (p12-m1-m2-al1)/(2*m1);
      c0 = Cplx(0);
      for(i = 0; i < 2; i++)
      {
	for(j = 0; j < 2; j++)
	  c0 += (1-2*i)*( li2(1.+x[i][j]/al) - li2(1.+x[i][j])
		  		+ Eta(1./al,-x[i][j]) * log(1.+x[i][j]/al) );
      }
      c0 += log(al) * log(m0/m1);
      c0 /= (p20-p12);
      return c0;
    }
  }

// Regular C0 function with p10,p12,p20 != 0
  al = la(p10,p12,p20);
  if(abs(al) <= ZERO_LIMIT)
    al = Cplx(0);
  al = sqrt(al);
//  al = sla(p10,p12,p20);
  al0 = sqrt( la(p12,m1,m2) * (1.+ IEPS*sign0(p12)) );
  al1 = sqrt( la(p20,m2,m0) * (1.+ IEPS*sign0(p20)) );
  al2 = sqrt( la(p10,m1,m0) * (1.+ IEPS*sign0(p10)) );

  if(al0 == Cplx(0)) al0 = IEPS*abs(al);
  if(al1 == Cplx(0)) al1 = IEPS*abs(al);
  if(al2 == Cplx(0)) al2 = IEPS*abs(al);
	  
  for(j = 0; j < 2; j++)
  {
    x[0][j] = (p12-m1+m2 + (2*j-1)*al0)/(2*p12);
    x[1][j] = (p20-m2+m0 + (2*j-1)*al1)/(2*p20);
    x[2][j] = (p10-m0+m1 + (2*j-1)*al2)/(2*p10);
  } 


  if(abs(al) <= ZEROS_LIMIT)
  {
// special case of threshold momentum configuration
    y0[0] = (p12*(p12-p10-p20+2*m0-m1-m2) - (p20-p10)*(m1-m2))/(2*p12);
    y0[1] = (p20*(p20-p12-p10+2*m1-m2-m0) - (p10-p12)*(m2-m0))/(2*p20);
    y0[2] = (p10*(p10-p20-p12+2*m2-m0-m1) - (p12-p20)*(m0-m1))/(2*p10);
    if(p12 <= Double(0) || p20 <= Double(0) || p10 <= Double(0))
    {
      Errstream << "Tachyonic external momenta in threshold configuration encountered! Can this be?" << endl;
      return(Cplx(HUGE_VAL));
    }
    c0 = Cplx(0);
    for(i = 0; i < 3; i++)
    {
      for(j = 0; j < 2; j++)
	c0 += ( (1-x[i][j])*log((1-x[i][j])/y0[i])
		+ x[i][j]*log(-x[i][j]/y0[i])
	  /*      - Eta(1.-x[i][j], 1/y0[i]) * (1.-x[i][j])
	        - Eta(-x[i][j], 1/y0[i]) * x[i][j]*/ ) / y0[i];
      c0 += Eta(-x[i][0],-x[i][1]) / y0[i];
    }
    return c0;
  }
  
// general momentum configuration

  y0[0] = (p12*(p12-p10-p20+2*m0-m1-m2) - (p20-p10)*(m1-m2)
	    + al*(p12-m1+m2))/(2*al*p12);
  y0[1] = (p20*(p20-p12-p10+2*m1-m2-m0) - (p10-p12)*(m2-m0)
	    + al*(p20-m2+m0))/(2*al*p20);
  y0[2] = (p10*(p10-p20-p12+2*m2-m0-m1) - (p12-p20)*(m0-m1)
	    + al*(p10-m0+m1))/(2*al*p10);
  
  for(i = 0; i < 3; i++)
  {
    for(j = 0; j < 2; j++)
      y[i][j] = y0[i] - x[i][j];
    thp[i] = Double(0);
    thy[i] = Double(0);
    if(imag(y[i][0]*y[i][1]) <= Double(0))
      thy[i] = Double(1);
  }
  if(p12 <= Double(0)) thp[0] = Double(1);
  if(p20 <= Double(0)) thp[1] = Double(1);
  if(p10 <= Double(0)) thp[2] = Double(1);
  
  c0 = Cplx(0);
  for(i = 0; i < 3; i++)
  {
    for(j = 0; j < 2; j++)
//    {
      c0 += li2((y0[i]-1.)/y[i][j]) - li2(y0[i]/y[i][j])
	  /*    + Eta(1.-x[i][j], 1/y[i][j]) * log((y0[i]-1.)/y[i][j])
	      - Eta(  -x[i][j], 1/y[i][j]) * log(y0[i]/y[i][j])*/;
//    }
    if(abs(real(y0[i]-1.)) < ZEROS_LIMIT){
         c0 -= Cplx(0);}
    else{
         c0 -= log((y0[i]-1.)/y0[i]) * (Eta(-x[i][0],-x[i][1]) - Eta(y[i][0],y[i][1])
	    				- TwoPiI * thp[i] * thy[i]);}
  }

  c0 /= al;
  return c0;
}


#ifdef QUADPRECISION
static const Cplx Idel(Double(0), SDouble("1e-10"));
#else
static const Cplx Idel(0, 1e-6);
#endif

// the derivative of C0 with respect to the third (squared) mass
Cplx DM3C0(Double p10, Double p20, Double p12, Double m1, Double m2, Double m3)
{
  Cplx m3c = m3*(1-Idel);  // add small I*eps term to avoid 0/0 instabilities
  Cplx den1 = (sqr(m1) + sqr(m3c - p12) - 2*m1*(m3c + p12)),
       den2 = (sqr(m2) + sqr(m3c - p20) - 2*m2*(m3c + p20));
  return (((8*m1*(-(powint(m1,3)*p20*(m2 - m3 + p20)) + 
	     p10*(-(sqr(m2 - m3)*(sqr(m1 - m3) - (m1 + m3)*p12)) + 
        	(m2 + m3)*(-m1 + m3 - p12)*(-m1 + m3 + p12)*p20 - 
        	(m1 + m3 - p12)*p12*sqr(p20)) + 
	     sqr(m1)*(m2*(m3 - p20)*(-2*p12 + p20) + sqr(m2)*(p12 + p20) + 
        	(m3 - p20)*(m3*(p12 - 2*p20) - 3*p12*p20)) + 
	     m2*(m3 - p12)*(sqr(m2)*p12 + 
        	(p12 - p20)*(sqr(m3) + m3*p20 - 2*p12*p20) + 
        	m2*(-3*p12*p20 + m3*(-2*p12 + p20))) - 
	     m1*(powint(m2,3)*p12 - sqr(m2)*(m3 - p12)*(p12 - 2*p20) + 
        	(m3 - p20)*(p12 - p20)*(m3*(m3 + p12) - 2*p12*p20) - 
        	m2*(sqr(m3)*(p12 + p20) - 3*p12*p20*(p12 + p20) + 
        	   2*m3*(sqr(p12) + sqr(p20))))))/(den1*den2) +
           (4*sqr(m1)*(m1*(p20 + p12 - p10) + p12*(-2*m2 + p20 - p12 + p10) +
                       m3*(-p20 + p12 + p10))*B0A(0,m1,m1))/den1 +
           (4*m1*m2*(-2*m1*p20 + m3*p20 - sqr(p20) - m3*p12 + p20*p12 +
                     m2*(p20 + p12 - p10) + m3*p10 + p20*p10)*B0A(0,m2,m2))/
           den2 +
           (((-m2 + m3 + p20)*(m1 + m3 - p12)*(sqr(m1) + sqr(m2 - p10) - 2*m1*(m2 + p10)))/den2 +
            ((m1*m1 + (m2 - p10)*(m3 - p12) - m1*(m2 + m3 + p10 + p12 - 2*p20))*
             (2*m1*(-m1 + m3 + p12)/den1 + (m1 + m2 - p10)*(m2 - m3 - p20)/den2)))*B0A(0,m3,m3) +
           (4*m1*p20*(-sqr(m2) + m1*(m2 - m3 + p20) + (m3 - p20)*p10 +
                      m2*(m3 + p20 - 2*p12 + p10))*B0A(p20,m2,m3))/den2 -
           (4*m1*p12*(sqr(m1) + (m3 - p12)*(m2 - p10) -
                     m1*(m2 + m3 - 2*p20 + p12 + p10))*B0A(p12,m1,m3))/den1 +
           4*m1*p10*B0A(p10,m1,m2))/
          (4.*m1*(sqr(m1)*p20 + m2*(m3c - p12)*(p20 - p12) + sqr(m2)*p12 -
                  m2*(m3c + p12)*p10 + p10*((m3c - p20)*(m3c - p12) + m3c*p10) -
                  m1*((m3c - p20)*(p20 - p12) + m2*(p20 + p12 - p10) + (m3c + p20)*p10))));
}

