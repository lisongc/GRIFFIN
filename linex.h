/* linex.h: header file for linex.cc */

double linex1d(double *axis1, int len1, double pnt1, double *data);
double linex2d(double *axis1, int len1, double *axis2, int len2, 
               double pnt1, double pnt2, double *data);
double linex3d(double *axis1, int len1, 
               double *axis2, int len2, 
               double *axis3, int len3, 
               double pnt1, double pnt2, double pnt3, double *data);
