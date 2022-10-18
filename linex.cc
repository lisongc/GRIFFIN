/*-----------------------------------------------------------------------------
linex.cc
Ayres Freitas (afreitas@pitt.edu)
last revision: 23 Nov 2021
-------------------------------------------------------------------------------
perform linear interpolation and extrapolation for 1D, 2D and 3D grids
-----------------------------------------------------------------------------*/

#include <iostream>
using namespace std;

double linex1d(double *axis1, int len1, double x, double *data)
{
  double x0, x1;
  int i;
  
  for(i = 1; i < len1-1 && x > axis1[i]; i++);
  x0 = axis1[i-1];
  x1 = axis1[i];
  return(data[i-1] + (x-x0)*(data[i]-data[i-1])/(x1-x0));
}

double linex2d(double *axis1, int len1, double *axis2, int len2, 
               double x, double y, double *data)
{
  double x0, x1, y0, y1;
  int i, j;
  
  for(i = 1; i < len1-1 && x > axis1[i]; i++);
  x0 = axis1[i-1];
  x1 = axis1[i];
  for(j = 1; j < len2-1 && y > axis2[j]; j++);
  y0 = axis2[j-1];
  y1 = axis2[j];
  
  return ( (x1-x)*(y1-y)*data[(i-1)*len2+(j-1)]
          +(x1-x)*(y-y0)*data[(i-1)*len2+j]
          +(x-x0)*(y1-y)*data[i*len2+(j-1)]
          +(x-x0)*(y-y0)*data[i*len2+j] )/((x1-x0)*(y1-y0));
}

double linex3d(double *axis1, int len1, double *axis2, int len2, 
               double *axis3, int len3, double x, double y, double z, 
	       double *data)
{
  double x0, x1, y0, y1, z0, z1;
  int i, j, k;
  
  for(i = 1; i < len1-1 && x > axis1[i]; i++);
  x0 = axis1[i-1];
  x1 = axis1[i];
  for(j = 1; j < len2-1 && y > axis2[j]; j++);
  y0 = axis2[j-1];
  y1 = axis2[j];
  for(k = 1; k < len3-1 && z > axis3[k]; k++);
  z0 = axis3[k-1];
  z1 = axis3[k];
  
  return ( (x1-x)*(y1-y)*(z1-z)*data[(i-1)*len2*len3+(j-1)*len3+(k-1)]
          +(x1-x)*(y1-y)*(z-z0)*data[(i-1)*len2*len3+(j-1)*len3+k]
	  +(x1-x)*(y-y0)*(z1-z)*data[(i-1)*len2*len3+j*len3+(k-1)]
          +(x1-x)*(y-y0)*(z-z0)*data[(i-1)*len2*len3+j*len3+k]
	  +(x-x0)*(y1-y)*(z1-z)*data[i*len2*len3+(j-1)*len3+(k-1)]
          +(x-x0)*(y1-y)*(z-z0)*data[i*len2*len3+(j-1)*len3+k]
	  +(x-x0)*(y-y0)*(z1-z)*data[i*len2*len3+j*len3+(k-1)]
          +(x-x0)*(y-y0)*(z-z0)*data[i*len2*len3+j*len3+k] )/
	 ((x1-x0)*(y1-y0)*(z1-z0));
}
