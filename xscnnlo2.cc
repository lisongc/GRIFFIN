/*-----------------------------------------------------------------------------
xscnnlo2.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 16 Apr 2026
-------------------------------------------------------------------------------
matrix element class needed for description of cross-section at O(Nf^2 al^2) 
and O(Nf al^2) precision away from Z pole
-----------------------------------------------------------------------------*/

#include "xscnnlo2.h"
#include "ff0.h"
#include "ff.h"
#include "oneloop.h"
#include "linex.h"

namespace griffin {

#include "aa2foffz.in"

Cplx mat_SMNNLO2::compoffZ2f(double sval) const
{
  double mz = ival->get(MZ);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  Cplx zie1 = Cplx(rz1fs(it,iff,sval,*ival), iz1fs(it,iff,sval,*ival)),
       zjf1 = Cplx(rz1fs(ot,off,sval,*ival), iz1fs(ot,off,sval,*ival)),
       zie1z = Cplx(rz1f(it,iff,*ival), iz1f(it,iff,*ival)),
       zjf1z = Cplx(rz1f(ot,off,*ival), iz1f(ot,off,*ival)),
       gie1 = Cplx(rg1fs(it,iff,sval,*ival), ig1fs(it,iff,sval,*ival)),
       gjf1 = Cplx(rg1fs(ot,off,sval,*ival), ig1fs(ot,off,sval,*ival)),
       sa1 = Cplx(rsg1fs(sval,*ival), isg1fs(sval,*ival)),
       ztie2 = zt2f(it,iff,sval,*ival),
       ztjf2 = zt2f(ot,off,sval,*ival),
       gie2 = g2f(it,iff,sval,*ival),
       gjf2 = g2f(ot,off,sval,*ival),
       sa2 = sg2f(sval,*ival),
       szt2 = szt2f(sval,*ival);
  Cplx mats1 = (zie0*ztjf2 + ztie2*zjf0 + zie1*zjf1 - zie1z*zjf1z
  	        - zie0*zjf0*szt2/(sval-mz*mz))/(sval-mz*mz)
	       + (gie0*gjf2 + gie2*gjf0 - gie0*gjf0*sa2/sval 
	          + (gie1-gie0*sa1/sval)*(gjf1-gjf0*sa1/sval))/sval;
  return(mats1);
}

Cplx mat_SMNNLO2::resoffZ2f(void) const
{
  if(it==ot)
    cerr << "The case of identical initial and final flavors has not been implemented yet" << endl;
  double mz = ival->get(MZ);
  if(abs(1-s/(mz*mz)) < 4e-4)
  {
    double del = s/(mz*mz)-1;
    Cplx mats1a= compoffZ2f(mz*mz*(1-4.1e-4));
    Cplx mats1b= compoffZ2f(mz*mz*(1+4.1e-4));
    return((mats1b*(del+4.1e-4) - mats1a*(del-4.1e-4))/(2*4.1e-4));
  }
  return compoffZ2f(s);
}


#include "aafboffz.grid"

Cplx linex3dc(double *axis1, int len1, 
              double *axis2, int len2, 
              double *axis3, int len3, 
              double pnt1, double pnt2, double pnt3, 
	      double *datar, double *datai)
{
  return(Cplx(
    linex3d(axis1, len1, axis2, len2, axis3, len3, pnt1, pnt2, pnt3, datar),
    linex3d(axis1, len1, axis2, len2, axis3, len3, pnt1, pnt2, pnt3, datai)));
}
              
Cplx mat_SMNNLO2::compoffZfb(double sval) const
{
  double mz = ival->get(MZ);
  double MZs = sqr(mz),
         MTs = sqr(ival->get(MT)),
         MHs = sqr(ival->get(MH)),
         AL2 = sqr(ival->get(al)),
         deltaAlpha = ival->get(Delal);
  double gie0 = g0(it,iff,*ival), gjf0 = g0(ot,off,*ival),
         zie0 = z0(it,iff,*ival), zjf0 = z0(ot,off,*ival);
  Cplx zie1f = Cplx(rz1fs(it,iff,sval,*ival), iz1fs(it,iff,sval,*ival)),
       zjf1f = Cplx(rz1fs(ot,off,sval,*ival), iz1fs(ot,off,sval,*ival)),
       zie1b = Cplx(rz1bs(it,iff,sval,*ival), iz1bs(it,iff,sval,*ival)),
       zjf1b = Cplx(rz1bs(ot,off,sval,*ival), iz1bs(ot,off,sval,*ival)),
       zie1fz = Cplx(rz1f(it,iff,*ival), iz1f(it,iff,*ival)),
       zjf1fz = Cplx(rz1f(ot,off,*ival), iz1f(ot,off,*ival)),
       zie1bz = Cplx(rz1b(it,iff,*ival), iz1b(it,iff,*ival)),
       zjf1bz = Cplx(rz1b(ot,off,*ival), iz1b(ot,off,*ival)),
       gie1f = Cplx(rg1fs(it,iff,sval,*ival), ig1fs(it,iff,sval,*ival)),
       gjf1f = Cplx(rg1fs(ot,off,sval,*ival), ig1fs(ot,off,sval,*ival)),
       gie1b = Cplx(rg1bs(it,iff,sval,*ival), ig1bs(it,iff,sval,*ival)),
       gjf1b = Cplx(rg1bs(ot,off,sval,*ival), ig1bs(ot,off,sval,*ival)),
       sa1f = Cplx(rsg1fs(sval,*ival), isg1fs(sval,*ival)),
       sa1b = Cplx(rsg1bs(sval,*ival), isg1bs(sval,*ival)),
       sz1b = Cplx(rsz1bs(sval,*ival), isz1bs(sval,*ival)),
       sz1bz = Cplx(rsz1b(*ival), isz1b(*ival)),
       szp1bz = Cplx(rsz1bp(*ival), isz1bp(*ival)),
       ztie2, ztjf2, gie2, gjf2, sa2, szt2;
/*  szt2= linex3d(axis1aafb, sizeof(axis1aafb)/sizeof(double),
  		axis2aafb, sizeof(axis2aafb)/sizeof(double),
  		axis3aafb, sizeof(axis3aafb)/sizeof(double),
		sval/MZs, sqr(log10(MHs)/2), MTs/MZs, &dataaafbSZ[0][0][0]);
  sa2 = linex3d(axis1aafb, sizeof(axis1aafb)/sizeof(double),
  		axis2aafb, sizeof(axis2aafb)/sizeof(double),
  		axis3aafb, sizeof(axis3aafb)/sizeof(double),
		sval/MZs, sqr(log10(MHs)/2), MTs/MZs, &dataaafbSA[0][0][0]);*/
  double *zdatar, *zdatai, *gdatar, *gdatai;
  switch(it) {
   case ELE:
   case MUO:
   case TAU: 
    if(iff==VEC)
    {
      zdatar = &dataraafbVZl[0][0][0];  zdatai = &dataiaafbVZl[0][0][0];
      gdatar = &dataraafbVAl[0][0][0];  gdatai = &dataiaafbVAl[0][0][0];
    }
    else
    {
      zdatar = &dataraafbAZl[0][0][0];  zdatai = &dataiaafbAZl[0][0][0];
      gdatar = &dataraafbAAl[0][0][0];  gdatai = &dataiaafbAAl[0][0][0];
    }
    break;
   case UQU:
   case CQU:
    if(iff==VEC)
    {
      zdatar = &dataraafbVZu[0][0][0];  zdatai = &dataiaafbVZu[0][0][0];
      gdatar = &dataraafbVAu[0][0][0];  gdatai = &dataiaafbVAu[0][0][0];
    }
    else
    {
      zdatar = &dataraafbAZu[0][0][0];  zdatai = &dataiaafbAZu[0][0][0];
      gdatar = &dataraafbAAu[0][0][0];  gdatai = &dataiaafbAAu[0][0][0];
    }
    break;
   case DQU:
   case SQU:
    if(iff==VEC)
    {
      zdatar = &dataraafbVZd[0][0][0];  zdatai = &dataiaafbVZd[0][0][0];
      gdatar = &dataraafbVAd[0][0][0];  gdatai = &dataiaafbVAd[0][0][0];
    }
    else
    {
      zdatar = &dataraafbAZd[0][0][0];  zdatai = &dataiaafbAZd[0][0][0];
      gdatar = &dataraafbAAd[0][0][0];  gdatai = &dataiaafbAAd[0][0][0];
    }
    break;
   default: cerr << "The case of neutrinos or bottom quarks has not been implemented yet" << endl;
    return(0);
  }
  ztie2= linex3dc(axis1aafb, sizeof(axis1aafb)/sizeof(double),
  		  axis2aafb, sizeof(axis2aafb)/sizeof(double),
  		  axis3aafb, sizeof(axis3aafb)/sizeof(double),
		  sval/MZs, sqr(log10(MHs)/2), MTs/MZs, zdatar, zdatai);
  gie2 = linex3dc(axis1aafb, sizeof(axis1aafb)/sizeof(double),
  		  axis2aafb, sizeof(axis2aafb)/sizeof(double),
  		  axis3aafb, sizeof(axis3aafb)/sizeof(double),
		  sval/MZs, sqr(log10(MHs)/2), MTs/MZs, gdatar, gdatai);
  switch(ot) {
   case ELE:
   case MUO:
   case TAU: 
    if(off==VEC)
    {
      zdatar = &dataraafbVZl[0][0][0];  zdatai = &dataiaafbVZl[0][0][0];
      gdatar = &dataraafbVAl[0][0][0];  gdatai = &dataiaafbVAl[0][0][0];
    }
    else
    {
      zdatar = &dataraafbAZl[0][0][0];  zdatai = &dataiaafbAZl[0][0][0];
      gdatar = &dataraafbAAl[0][0][0];  gdatai = &dataiaafbAAl[0][0][0];
    }
    break;
   case UQU:
   case CQU:
    if(off==VEC)
    {
      zdatar = &dataraafbVZu[0][0][0];  zdatai = &dataiaafbVZu[0][0][0];
      gdatar = &dataraafbVAu[0][0][0];  gdatai = &dataiaafbVAu[0][0][0];
    }
    else
    {
      zdatar = &dataraafbAZu[0][0][0];  zdatai = &dataiaafbAZu[0][0][0];
      gdatar = &dataraafbAAu[0][0][0];  gdatai = &dataiaafbAAu[0][0][0];
    }
    break;
   case DQU:
   case SQU:
    if(off==VEC)
    {
      zdatar = &dataraafbVZd[0][0][0];  zdatai = &dataiaafbVZd[0][0][0];
      gdatar = &dataraafbVAd[0][0][0];  gdatai = &dataiaafbVAd[0][0][0];
    }
    else
    {
      zdatar = &dataraafbAZd[0][0][0];  zdatai = &dataiaafbAZd[0][0][0];
      gdatar = &dataraafbAAd[0][0][0];  gdatai = &dataiaafbAAd[0][0][0];
    }
    break;
   default: cerr << "The case of neutrinos or bottom quarks has not been implemented yet" << endl;
    return(0);
  }
  ztjf2= linex3dc(axis1aafb, sizeof(axis1aafb)/sizeof(double),
  		  axis2aafb, sizeof(axis2aafb)/sizeof(double),
  		  axis3aafb, sizeof(axis3aafb)/sizeof(double),
		  sval/MZs, sqr(log10(MHs)/2), MTs/MZs, zdatar, zdatai);
  gjf2 = linex3dc(axis1aafb, sizeof(axis1aafb)/sizeof(double),
  		  axis2aafb, sizeof(axis2aafb)/sizeof(double),
  		  axis3aafb, sizeof(axis3aafb)/sizeof(double),
		  sval/MZs, sqr(log10(MHs)/2), MTs/MZs, gdatar, gdatai);
  Cplx mats1 = zie0*AL2*ztjf2*zjf0 + AL2*ztie2*zie0*zjf0 
                + (zie1f*zjf1b + zie1b*zjf1f 
  		- zie1fz*zjf1bz- zie1bz*zjf1fz 
  	        /*- zie0*zjf0*AL2*szt2/(sval-MZs)*/)/(sval-MZs)
	       + (gie0*(AL2*gjf2*g0(ot,VEC,*ival)) 
	          + (AL2*gie2*g0(it,VEC,*ival))*gjf0 /*- gie0*gjf0*AL2*sa2/sval*/ 
	          + (gie1f-gie0*sa1f/sval)*(gjf1b-gjf0*sa1b/sval)
		  + (gie1b-gie0*sa1b/sval)*(gjf1f-gjf0*sa1f/sval))/sval;
// add deltalpha contribution not contained in grids:
  mats1 += deltaAlpha*((1.5*(zie0*(zjf1b-zjf1bz) + (zie1b-zie1bz)*zjf0) 
  			- 2*zie0*zjf0*((sz1b-sz1bz)/(sval-MZs) - szp1bz))/
		       (sval-MZs)
          	       + (gie0*gjf1b + gie1b*gjf0)/sval);
  return(mats1);
}

#include "aafbbox.grid"
#include "ffw1.in"

Cplx mat_SMNNLO2::compboxfb(double sval) const
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         AL = ival->get(al);
  double SWsq = 1-MWs/MZs, CWsq = MWs/MZs;
  Cplx zzz0, zzz1, reszzz, www, reswww, azz0, azz1, resazz, aza0, aza1, resaza,
       zza0, zza1, reszza, resctv;
  
  Cplx vvbzz = B1s(it,ot,VEC,VEC,s,-s/2*(1-cost),*ival,0,0,1,1,0),
       avbzz = B1s(it,ot,VEC,AXV,s,-s/2*(1-cost),*ival,0,0,1,1,0),
       vabzz = B1s(it,ot,AXV,VEC,s,-s/2*(1-cost),*ival,0,0,1,1,0),
       aabzz = B1s(it,ot,AXV,AXV,s,-s/2*(1-cost),*ival,0,0,1,1,0),
       vvbww = B1s(it,ot,VEC,VEC,s,-s/2*(1-cost),*ival,0,0,1,0,1),
       avbww = B1s(it,ot,VEC,AXV,s,-s/2*(1-cost),*ival,0,0,1,0,1),
       vabww = B1s(it,ot,AXV,VEC,s,-s/2*(1-cost),*ival,0,0,1,0,1),
       aabww = B1s(it,ot,AXV,AXV,s,-s/2*(1-cost),*ival,0,0,1,0,1),
       vvbaz = B1s(it,ot,VEC,VEC,s,-s/2*(1-cost),*ival,0,1,1,0,0),
       avbaz = B1s(it,ot,VEC,AXV,s,-s/2*(1-cost),*ival,0,1,1,0,0),
       vabaz = B1s(it,ot,AXV,VEC,s,-s/2*(1-cost),*ival,0,1,1,0,0),
       aabaz = B1s(it,ot,AXV,AXV,s,-s/2*(1-cost),*ival,0,1,1,0,0);
  inval ival2(*ival);
  ival2.set(MB, 0);
  double ziv1 = rz1f(it,VEC,ival2), 
         zia1 = rz1f(it,AXV,ival2),   
	 zov1 = rz1f(ot,VEC,ival2), 
         zoa1 = rz1f(ot,AXV,ival2),
         w1 = w1f(ival2);
  double ziv0 = vz0(it,ival2), zia0 = az0(it,ival2),
         zov0 = vz0(ot,ival2), zoa0 = az0(ot,ival2);
  Cplx llctv = (vvbzz+avbzz+vabzz+aabzz)
  		*2*((ziv1+zia1)/(ziv0+zia0) + (zov1+zoa1)/(zov0+zoa0))
	     + (vvbaz+avbaz+vabaz+aabaz)
  		*((ziv1+zia1)/(ziv0+zia0) + (zov1+zoa1)/(zov0+zoa0))
	     + (vvbww+avbww+vabww+aabww) * 4*w1,
       lrctv = (vvbzz+avbzz-vabzz-aabzz)
  		*2*((ziv1-zia1)/(ziv0-zia0) + (zov1+zoa1)/(zov0+zoa0))
	     + (vvbaz+avbaz-vabaz-aabaz)
  		*((ziv1-zia1)/(ziv0-zia0) + (zov1+zoa1)/(zov0+zoa0)),
       rlctv = (vvbzz-avbzz+vabzz-aabzz)
  		*2*((ziv1+zia1)/(ziv0+zia0) + (zov1-zoa1)/(zov0-zoa0))
	     + (vvbaz-avbaz+vabaz-aabaz)
  		*((ziv1+zia1)/(ziv0+zia0) + (zov1-zoa1)/(zov0-zoa0)),
       rrctv = (vvbzz-avbzz-vabzz+aabzz)
  		*2*((ziv1-zia1)/(ziv0-zia0) + (zov1-zoa1)/(zov0-zoa0))
	     + (vvbaz-avbaz-vabaz+aabaz)
  		*((ziv1-zia1)/(ziv0-zia0) + (zov1-zoa1)/(zov0-zoa0));
  
//  zzz0 = Cplx(6.869001176e12,-1.144367670e13);
//  zzz1 = Cplx(-3.150199579e13,2.291958786e12);
//  azz0 = Cplx(-2.218979533e5,-1.254285595e5);
//  azz1 = Cplx(-1.227117258e5,1.453078525e4);
//  aza0 = Cplx(-3.611404327e2,6.667992962e1);
//  aza1 = Cplx(-1.186664883e2,9.953365640e1);
//  zza0 = Cplx(5.861748731e11,-1.208892733e11);
//  zza1 = Cplx(3.350031957e11,-5.196671503e10);
  
  zzz0 = linex3dc(axis1aafbzzz, sizeof(axis1aafbzzz)/sizeof(double),
                  axis2aafbzzz, sizeof(axis2aafbzzz)/sizeof(double),
		  axis3aafbzzz, sizeof(axis3aafbzzz)/sizeof(double),
		  s/MZs, acos(cost), MTs/MZs, 
		  &dataraafbzzz0[0][0][0], &dataiaafbzzz0[0][0][0]);
  zzz1 = linex3dc(axis1aafbzzz, sizeof(axis1aafbzzz)/sizeof(double),
                  axis2aafbzzz, sizeof(axis2aafbzzz)/sizeof(double),
		  axis3aafbzzz, sizeof(axis3aafbzzz)/sizeof(double),
		  s/MZs, acos(cost), MTs/MZs, 
		  &dataraafbzzz1[0][0][0], &dataiaafbzzz1[0][0][0]);
  azz0 = linex3dc(axis1aafbazz, sizeof(axis1aafbazz)/sizeof(double),
                  axis2aafbazz, sizeof(axis2aafbazz)/sizeof(double),
		  axis3aafbazz, sizeof(axis3aafbazz)/sizeof(double),
		  s/MZs, acos(cost), MTs/MZs, 
		  &dataraafbazz0[0][0][0], &dataiaafbazz0[0][0][0]);
  azz1 = linex3dc(axis1aafbazz, sizeof(axis1aafbazz)/sizeof(double),
                  axis2aafbazz, sizeof(axis2aafbazz)/sizeof(double),
		  axis3aafbazz, sizeof(axis3aafbazz)/sizeof(double),
		  s/MZs, acos(cost), MTs/MZs, 
		  &dataraafbazz1[0][0][0], &dataiaafbazz1[0][0][0]);
  aza0 = linex3dc(axis1aafbaza, sizeof(axis1aafbaza)/sizeof(double),
                  axis2aafbaza, sizeof(axis2aafbaza)/sizeof(double),
		  axis3aafbaza, sizeof(axis3aafbaza)/sizeof(double),
		  s/MZs, acos(cost), MTs/MZs, 
		  &dataraafbaza0[0][0][0], &dataiaafbaza0[0][0][0]);
  aza1 = linex3dc(axis1aafbaza, sizeof(axis1aafbaza)/sizeof(double),
                  axis2aafbaza, sizeof(axis2aafbaza)/sizeof(double),
		  axis3aafbaza, sizeof(axis3aafbaza)/sizeof(double),
		  s/MZs, acos(cost), MTs/MZs, 
		  &dataraafbaza1[0][0][0], &dataiaafbaza1[0][0][0]);
  zza0 = linex3dc(axis1aafbzza, sizeof(axis1aafbzza)/sizeof(double),
                  axis2aafbzza, sizeof(axis2aafbzza)/sizeof(double),
		  axis3aafbzza, sizeof(axis3aafbzza)/sizeof(double),
		  s/MZs, acos(cost), MTs/MZs, 
		  &dataraafbzza0[0][0][0], &dataiaafbzza0[0][0][0]);
  zza1 = linex3dc(axis1aafbzza, sizeof(axis1aafbzza)/sizeof(double),
                  axis2aafbzza, sizeof(axis2aafbzza)/sizeof(double),
		  axis3aafbzza, sizeof(axis3aafbzza)/sizeof(double),
		  s/MZs, acos(cost), MTs/MZs, 
		  &dataraafbzza1[0][0][0], &dataiaafbzza1[0][0][0]);
  
  if(iff==VEC)
  {
    if(off==VEC)
    {
      resctv = (llctv + lrctv + rlctv + rrctv)/4.;
      reszzz = (powint(CWsq,-3)*powint(SWsq,-3)*
       (2*(zzz0 - 2*cost*zzz1 + zzz0*sqr(cost))*sqr(SWsq)*sqr(Qf[it])*
          (-2*SWsq*I3f[ot]*Qf[ot] + sqr(I3f[ot]) + 2*sqr(SWsq)*sqr(Qf[ot])) - 
         2*SWsq*I3f[it]*Qf[it]*(-2*SWsq*(zzz0 + zzz1)*I3f[ot]*Qf[ot]*
             sqr(-1 + cost) + (zzz0 + zzz1)*sqr(-1 + cost)*sqr(I3f[ot]) + 
            2*(zzz0 - 2*cost*zzz1 + zzz0*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])) + 
         sqr(I3f[it])*(-2*SWsq*(zzz0 + zzz1)*I3f[ot]*Qf[ot]*sqr(-1 + cost) + 
            (zzz0 + zzz1)*sqr(-1 + cost)*sqr(I3f[ot]) + 
            2*(zzz0 - 2*cost*zzz1 + zzz0*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])))*
       trip(AL))/(72.*Pi*sqr(MZs)*sqr(s)*sqr(-1 + sqr(cost)));
      resazz = ((-5 + 8*CWsq)*(Qf[ot]*(-((azz0 + azz1)*I3f[ot]*sqr(-1 + cost)) + 
            2*SWsq*Qf[ot]*(azz0 - 2*azz1*cost + azz0*sqr(cost)))*sqr(I3f[it]) + 
         2*SWsq*sqr(Qf[it])*(-3*SWsq*I3f[ot]*Qf[ot]*
             (azz0 - 2*azz1*cost + azz0*sqr(cost)) + 
            (azz0 - 2*azz1*cost + azz0*sqr(cost))*sqr(I3f[ot]) + 
            SWsq*((-(azz1*(1 + cost*(6 + cost))) + 
                  azz0*(3 + cost*(2 + 3*cost)))*SWsq - 
               (azz0 + azz1)*(-1 + CWsq)*sqr(-1 + cost))*sqr(Qf[ot])) + 
         I3f[it]*Qf[it]*(4*(azz0 + azz1)*SWsq*I3f[ot]*Qf[ot]*sqr(-1 + cost) - 
            (azz0 + azz1)*sqr(-1 + cost)*sqr(I3f[ot]) - 
            6*(azz0 - 2*azz1*cost + azz0*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])))*
       trip(AL))/(18.*Pi*sqr(CWsq)*sqr(s)*sqr(SWsq)*sqr(-1 + sqr(cost)));
      resaza = -((-5 + 8*CWsq)*Qf[it]*Qf[ot]*
        (I3f[ot]*Qf[it] + (I3f[it] - 4*SWsq*Qf[it])*Qf[ot])*
        (aza0 - 2*aza1*cost + aza0*sqr(cost))*trip(AL))/
       (18.*(-1 + cost)*CWsq*Pi*s*SWsq*sqr(1 + cost));
      reszza = (Qf[it]*Qf[ot]*(-2*SWsq*Qf[it]*(I3f[ot] - 2*SWsq*Qf[ot])*
          (zza0 - 2*cost*zza1 + zza0*sqr(cost)) + 
         I3f[it]*((zza0 + zza1)*I3f[ot]*sqr(-1 + cost) - 
            2*SWsq*Qf[ot]*(zza0 - 2*cost*zza1 + zza0*sqr(cost))))*trip(AL))/
       (144.*(-1 + cost)*Pi*s*sqr(1 + cost)*sqr(CWsq)*sqr(MZs)*sqr(SWsq));
    }
    else
    {
      resctv = (llctv + lrctv - rlctv - rrctv)/4.;
      reszzz = (powint(CWsq,-3)*powint(SWsq,-3)*
       (2*I3f[ot]*(I3f[ot] - 2*SWsq*Qf[ot])*
          (zzz0 - 2*cost*zzz1 + zzz0*sqr(cost))*sqr(SWsq)*sqr(Qf[it]) - 
         2*SWsq*I3f[it]*Qf[it]*(-2*SWsq*(zzz0 + zzz1)*I3f[ot]*Qf[ot]*
             sqr(-1 + cost) + (zzz0 + zzz1)*sqr(-1 + cost)*sqr(I3f[ot]) + 
            2*(-2*cost*zzz0 + zzz1 + zzz1*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])) + 
         sqr(I3f[it])*(-2*SWsq*(zzz0 + zzz1)*I3f[ot]*Qf[ot]*sqr(-1 + cost) + 
            (zzz0 + zzz1)*sqr(-1 + cost)*sqr(I3f[ot]) + 
            2*(-2*cost*zzz0 + zzz1 + zzz1*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])))*
       trip(AL))/(72.*Pi*sqr(MZs)*sqr(s)*sqr(-1 + sqr(cost)));
      resazz = ((-5 + 8*CWsq)*(Qf[ot]*(-((azz0 + azz1)*I3f[ot]*sqr(-1 + cost)) + 
            2*SWsq*Qf[ot]*(azz1 - 2*azz0*cost + azz1*sqr(cost)))*sqr(I3f[it]) + 
         2*SWsq*sqr(Qf[it])*(-3*SWsq*I3f[ot]*Qf[ot]*
             (azz0 - 2*azz1*cost + azz0*sqr(cost)) + 
            (azz0 - 2*azz1*cost + azz0*sqr(cost))*sqr(I3f[ot]) + 
            (azz0 + azz1)*SWsq*(-1 + CWsq + SWsq)*sqr(-1 + cost)*sqr(Qf[ot])) + 
         I3f[it]*Qf[it]*(4*(azz0 + azz1)*SWsq*I3f[ot]*Qf[ot]*sqr(-1 + cost) - 
            (azz0 + azz1)*sqr(-1 + cost)*sqr(I3f[ot]) - 
            6*(azz1 - 2*azz0*cost + azz1*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])))*
       trip(AL))/(18.*Pi*sqr(CWsq)*sqr(s)*sqr(SWsq)*sqr(-1 + sqr(cost)));
      resaza = -((-5 + 8*CWsq)*Qf[it]*Qf[ot]*
        (I3f[ot]*Qf[it]*(aza0 - 2*aza1*cost + aza0*sqr(cost)) + 
          I3f[it]*Qf[ot]*(aza1 - 2*aza0*cost + aza1*sqr(cost)))*trip(AL))/
       (18.*(-1 + cost)*CWsq*Pi*s*SWsq*sqr(1 + cost));
      reszza = (Qf[it]*Qf[ot]*(-2*SWsq*I3f[ot]*Qf[it]*
          (zza0 - 2*cost*zza1 + zza0*sqr(cost)) + 
         I3f[it]*((zza0 + zza1)*I3f[ot]*sqr(-1 + cost) - 
            2*SWsq*Qf[ot]*(-2*cost*zza0 + zza1 + zza1*sqr(cost))))*trip(AL))/
       (144.*(-1 + cost)*Pi*s*sqr(1 + cost)*sqr(CWsq)*sqr(MZs)*sqr(SWsq));
    }
  }
  else
  {
    if(off==VEC)
    {
      resctv = (llctv - lrctv + rlctv - rrctv)/4.;
      reszzz = (powint(CWsq,-3)*powint(SWsq,-3)*
       (2*I3f[ot]*(I3f[ot] - 2*SWsq*Qf[ot])*
          (-2*cost*zzz0 + zzz1 + zzz1*sqr(cost))*sqr(SWsq)*sqr(Qf[it]) - 
         2*SWsq*I3f[it]*Qf[it]*(-2*SWsq*(zzz0 + zzz1)*I3f[ot]*Qf[ot]*
             sqr(-1 + cost) + (zzz0 + zzz1)*sqr(-1 + cost)*sqr(I3f[ot]) + 
            2*(zzz0 - 2*cost*zzz1 + zzz0*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])) + 
         sqr(I3f[it])*(-2*SWsq*(zzz0 + zzz1)*I3f[ot]*Qf[ot]*sqr(-1 + cost) + 
            (zzz0 + zzz1)*sqr(-1 + cost)*sqr(I3f[ot]) + 
            2*(zzz0 - 2*cost*zzz1 + zzz0*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])))*
       trip(AL))/(72.*Pi*sqr(MZs)*sqr(s)*sqr(-1 + sqr(cost)));
      resazz = ((-5 + 8*CWsq)*(Qf[ot]*(-((azz0 + azz1)*I3f[ot]*sqr(-1 + cost)) + 
            2*SWsq*Qf[ot]*(azz0 - 2*azz1*cost + azz0*sqr(cost)))*sqr(I3f[it]) + 
         2*SWsq*sqr(Qf[it])*(-3*SWsq*I3f[ot]*Qf[ot]*
             (azz1 - 2*azz0*cost + azz1*sqr(cost)) + 
            (azz1 - 2*azz0*cost + azz1*sqr(cost))*sqr(I3f[ot]) + 
            (azz0 + azz1)*SWsq*(-1 + CWsq + SWsq)*sqr(-1 + cost)*sqr(Qf[ot])) + 
         I3f[it]*Qf[it]*(4*(azz0 + azz1)*SWsq*I3f[ot]*Qf[ot]*sqr(-1 + cost) - 
            (azz0 + azz1)*sqr(-1 + cost)*sqr(I3f[ot]) - 
            6*(azz0 - 2*azz1*cost + azz0*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])))*
       trip(AL))/(18.*Pi*sqr(CWsq)*sqr(s)*sqr(SWsq)*sqr(-1 + sqr(cost)));
      resaza = -((-5 + 8*CWsq)*Qf[it]*Qf[ot]*
        (I3f[it]*Qf[ot]*(aza0 - 2*aza1*cost + aza0*sqr(cost)) + 
          I3f[ot]*Qf[it]*(aza1 - 2*aza0*cost + aza1*sqr(cost)))*trip(AL))/
       (18.*(-1 + cost)*CWsq*Pi*s*SWsq*sqr(1 + cost));
      reszza = (Qf[it]*Qf[ot]*(-2*SWsq*I3f[ot]*Qf[it]*
          (-2*cost*zza0 + zza1 + zza1*sqr(cost)) + 
         I3f[it]*((zza0 + zza1)*I3f[ot]*sqr(-1 + cost) - 
            2*SWsq*Qf[ot]*(zza0 - 2*cost*zza1 + zza0*sqr(cost))))*trip(AL))/
       (144.*(-1 + cost)*Pi*s*sqr(1 + cost)*sqr(CWsq)*sqr(MZs)*sqr(SWsq));
    }
    else
    {
      resctv = (llctv - lrctv - rlctv + rrctv)/4.;
      reszzz = (powint(CWsq,-3)*powint(SWsq,-3)*
       (2*(-2*cost*zzz0 + zzz1 + zzz1*sqr(cost))*sqr(SWsq)*sqr(Qf[it])*
          (-2*SWsq*I3f[ot]*Qf[ot] + sqr(I3f[ot]) + 2*sqr(SWsq)*sqr(Qf[ot])) - 
         2*SWsq*I3f[it]*Qf[it]*(-2*SWsq*(zzz0 + zzz1)*I3f[ot]*Qf[ot]*
             sqr(-1 + cost) + (zzz0 + zzz1)*sqr(-1 + cost)*sqr(I3f[ot]) + 
            2*(-2*cost*zzz0 + zzz1 + zzz1*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])) + 
         sqr(I3f[it])*(-2*SWsq*(zzz0 + zzz1)*I3f[ot]*Qf[ot]*sqr(-1 + cost) + 
            (zzz0 + zzz1)*sqr(-1 + cost)*sqr(I3f[ot]) + 
            2*(-2*cost*zzz0 + zzz1 + zzz1*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])))*
       trip(AL))/(72.*Pi*sqr(MZs)*sqr(s)*sqr(-1 + sqr(cost)));
      resazz = ((-5 + 8*CWsq)*(Qf[ot]*(-((azz0 + azz1)*I3f[ot]*sqr(-1 + cost)) + 
            2*SWsq*Qf[ot]*(azz1 - 2*azz0*cost + azz1*sqr(cost)))*sqr(I3f[it]) + 
         2*SWsq*sqr(Qf[it])*(-3*SWsq*I3f[ot]*Qf[ot]*
             (azz1 - 2*azz0*cost + azz1*sqr(cost)) + 
            (azz1 - 2*azz0*cost + azz1*sqr(cost))*sqr(I3f[ot]) + 
            SWsq*(-((azz0*(1 + cost*(6 + cost)) - azz1*(3 + cost*(2 + 3*cost)))*
                  SWsq) - (azz0 + azz1)*(-1 + CWsq)*sqr(-1 + cost))*sqr(Qf[ot]))
           + I3f[it]*Qf[it]*(4*(azz0 + azz1)*SWsq*I3f[ot]*Qf[ot]*
             sqr(-1 + cost) - (azz0 + azz1)*sqr(-1 + cost)*sqr(I3f[ot]) - 
            6*(azz1 - 2*azz0*cost + azz1*sqr(cost))*sqr(SWsq)*sqr(Qf[ot])))*
       trip(AL))/(18.*Pi*sqr(CWsq)*sqr(s)*sqr(SWsq)*sqr(-1 + sqr(cost)));
      resaza = -((-5 + 8*CWsq)*Qf[it]*Qf[ot]*
        (I3f[ot]*Qf[it] + (I3f[it] - 4*SWsq*Qf[it])*Qf[ot])*
        (aza1 - 2*aza0*cost + aza1*sqr(cost))*trip(AL))/
       (18.*(-1 + cost)*CWsq*Pi*s*SWsq*sqr(1 + cost));
      reszza = (Qf[it]*Qf[ot]*(-2*SWsq*Qf[it]*(I3f[ot] - 2*SWsq*Qf[ot])*
          (-2*cost*zza0 + zza1 + zza1*sqr(cost)) + 
         I3f[it]*((zza0 + zza1)*I3f[ot]*sqr(-1 + cost) - 
            2*SWsq*Qf[ot]*(-2*cost*zza0 + zza1 + zza1*sqr(cost))))*trip(AL))/
       (144.*(-1 + cost)*Pi*s*sqr(1 + cost)*sqr(CWsq)*sqr(MZs)*sqr(SWsq));
    }
  }
  if(I3f[it]==I3f[ot])
//    www = Cplx(2.655607953e5,4.728079144e5);
    www = linex3dc(axis1aafbwww, sizeof(axis1aafbwww)/sizeof(double),
                   axis2aafbwww, sizeof(axis2aafbwww)/sizeof(double),
		   axis3aafbwww, sizeof(axis3aafbwww)/sizeof(double),
		   s/MZs, acos(cost), MTs/MZs, 
		   &dataraafbwww0[0][0][0], &dataiaafbwww0[0][0][0]);
  else
//    www = Cplx(-2.815833256e5,-4.565509878e5);
    www = linex3dc(axis1aafbwww, sizeof(axis1aafbwww)/sizeof(double),
                   axis2aafbwww, sizeof(axis2aafbwww)/sizeof(double),
		   axis3aafbwww, sizeof(axis3aafbwww)/sizeof(double),
		   s/MZs, acos(cost), MTs/MZs, 
		   &dataraafbwww1[0][0][0], &dataiaafbwww1[0][0][0]);
  reswww = (powint(MWs - MZs,-3)*trip(AL)*trip(MZs))/(16.*Pi*sqr(1+cost)*sqr(s))
  	   * www;
  return resctv+resaza+reszza+resazz+reszzz+reswww;
}

Cplx mat_SMNNLO2::resvpol(void) const
{
  double MZs = sqr(ival->get(MZ)),
         MWs = sqr(ival->get(MW)),
         MTs = sqr(ival->get(MT)),
         AL = ival->get(al);
  double SWsq = 1-MWs/MZs, CWsq = MWs/MZs, LMF, fac;
  Cplx aaa0, aaa1, resaaa, aaz0, aaz1, resaaz;
  int i;

//  aaa0 = Cplx(-6.409366571e1,-3.360793073e1);
//  aaa1 = Cplx(-1.288086522,-7.701101623e2);
//  aaa0 = Cplx(-3.436638402e1,2.827378516);
//  aaa1 = Cplx(7.599694173,-2.131747388e2);
//  aaz0 = Cplx(-1.974262937e5,-4.580946509e6);
//  aaz1 = Cplx(1.311778141e6,-1.615769407e5);
//  aaz0 = Cplx(-1.007277008e1,-2.337217606e2);
//  aaz1 = Cplx(6.692745622e1,-8.243721466);

  aaa0 = aaa1 = aaz0 = aaz1 = 0;
  for(i = ME; i <= MT; i++)
  {
    LMF = log(ival->get(i)/ival->get(MZ));  
    if(i>=ME && i<=ML) fac = 1;
    if(i>=MD && i<=MB) fac = 3*(1/3.)*(1/3.);
    if(i>=MU && i<=MT) fac = 3*(2/3.)*(2/3.);
    aaa0 += linex3dc(axis1aafbaaa, sizeof(axis1aafbaaa)/sizeof(double),
                     axis2aafbaaa, sizeof(axis2aafbaaa)/sizeof(double),
		     axis3aafbaaa, sizeof(axis3aafbaaa)/sizeof(double),
		     s/MZs, acos(cost), LMF, 
		     &dataraafbaaa0[0][0][0], &dataiaafbaaa0[0][0][0]) * fac;
    aaa1 += linex3dc(axis1aafbaaa, sizeof(axis1aafbaaa)/sizeof(double),
                     axis2aafbaaa, sizeof(axis2aafbaaa)/sizeof(double),
		     axis3aafbaaa, sizeof(axis3aafbaaa)/sizeof(double),
		     s/MZs, acos(cost), LMF, 
		     &dataraafbaaa1[0][0][0], &dataiaafbaaa1[0][0][0]) * fac;
    aaz0 += linex3dc(axis1aafbaaz, sizeof(axis1aafbaaz)/sizeof(double),
                     axis2aafbaaz, sizeof(axis2aafbaaz)/sizeof(double),
		     axis3aafbaaz, sizeof(axis3aafbaaz)/sizeof(double),
		     s/MZs, acos(cost), LMF, 
		     &dataraafbaaz0[0][0][0], &dataiaafbaaz0[0][0][0]) * fac;
    aaz1 += linex3dc(axis1aafbaaz, sizeof(axis1aafbaaz)/sizeof(double),
                     axis2aafbaaz, sizeof(axis2aafbaaz)/sizeof(double),
		     axis3aafbaaz, sizeof(axis3aafbaaz)/sizeof(double),
		     s/MZs, acos(cost), LMF, 
		     &dataraafbaaz1[0][0][0], &dataiaafbaaz1[0][0][0]) * fac;
  }

  if(iff==VEC)
  {
    if(off==VEC)
    {
      resaaa = (8*(aaa0 - 2*aaa1*cost + aaa0*sqr(cost))*sqr(Qf[it])*
       sqr(Qf[ot])*trip(AL))/(9.*Pi*s*sqr(-1 + sqr(cost)));
      resaaz = (4*Qf[it]*Qf[ot]*(-2*SWsq*(I3f[it] - 2*SWsq*Qf[it])*Qf[ot]*
          (aaz0 - 2*aaz1*cost + aaz0*sqr(cost)) + 
         I3f[ot]*((aaz0 + aaz1)*I3f[it]*sqr(-1 + cost) - 
            2*SWsq*Qf[it]*(aaz0 - 2*aaz1*cost + aaz0*sqr(cost))))*trip(AL))/
       (9.*CWsq*Pi*SWsq*s*sqr(-1 + sqr(cost)));
    }
    else
    {
      resaaa = 0;
      resaaz = (4*Qf[it]*Qf[ot]*(-2*SWsq*I3f[ot]*Qf[it]*
          (aaz0 - 2*aaz1*cost + aaz0*sqr(cost)) + 
         I3f[it]*((aaz0 + aaz1)*I3f[ot]*sqr(-1 + cost) - 
            2*SWsq*Qf[ot]*(aaz1 - 2*aaz0*cost + aaz1*sqr(cost))))*trip(AL))/
       (9.*CWsq*Pi*SWsq*s*sqr(-1 + sqr(cost)));
    }
  }
  else
  {
    if(off==VEC)
    {
      resaaa = 0;
      resaaz = (4*Qf[it]*Qf[ot]*(-2*SWsq*I3f[ot]*Qf[it]*
          (aaz1 - 2*aaz0*cost + aaz1*sqr(cost)) + 
         I3f[it]*((aaz0 + aaz1)*I3f[ot]*sqr(-1 + cost) - 
            2*SWsq*Qf[ot]*(aaz0 - 2*aaz1*cost + aaz0*sqr(cost))))*trip(AL))/
       (9.*CWsq*Pi*SWsq*s*sqr(-1 + sqr(cost)));
    }
    else
    {
      resaaa = (8*(aaa1 - 2*aaa0*cost + aaa1*sqr(cost))*sqr(Qf[it])*
       sqr(Qf[ot])*trip(AL))/(9.*Pi*s*sqr(-1 + sqr(cost)));
      resaaz = (4*Qf[it]*Qf[ot]*(-2*SWsq*(I3f[it] - 2*SWsq*Qf[it])*Qf[ot]*
          (aaz1 - 2*aaz0*cost + aaz1*sqr(cost)) + 
         I3f[ot]*((aaz0 + aaz1)*I3f[it]*sqr(-1 + cost) - 
            2*SWsq*Qf[it]*(aaz1 - 2*aaz0*cost + aaz1*sqr(cost))))*trip(AL))/
       (9.*CWsq*Pi*SWsq*s*sqr(-1 + sqr(cost)));
    }
  }
  return resaaa+resaaz;
}

Cplx mat_SMNNLO2::resoffZfb(void) const
{
  if(it==ot) {
    cerr << "The case of identical initial and final flavors has not been implemented yet" << endl;
    return(0);
  }
  if(it==BQU || it==NUE || it==NUM || it==NUT ||
     ot==BQU || ot==NUE || ot==NUM || ot==NUT)
  {
    cerr << "The case of neutrinos or bottom quarks has not been implemented yet" << endl;
    return(0);
  }
  return compboxfb(s)+compoffZfb(s);
}

Cplx mat_SMNNLO2::result(void) const
{
  double mz = ival->get(MZ),
         gz = ival->get(GamZ);
  Cplx sminuss0(s - mz*mz, mz*gz);
  return(coeffR()/sminuss0 + resoffZ() + resvpol());
}

} // namespace
