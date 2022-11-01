/* linex.h: header file for linex.cc */

// perform linear interpolation of 1D grid
// data tuples (xi, fi) are given as follows:
// axis1: array of x1, x2, ...
// len1: length of axis1 array
// data: array of f1=f(x1), f2=f(x2), ...
// pnt: variable x
// linex1d returns f(x)
double linex1d(double *axis1, int len1, double pnt1, double *data);

// perform linear interpolation of 2D grid
// data tuples (xi, yj, fij) are given as follows:
// axis1: array of x1, x2, ...
// len1: length of axis1 array
// axis2: array of y1, y2, ...
// len2: length of axis2 array
// data: array of f(x1,y1), f(x1,y2), ..., f(x1,y_len2), f(x2,y1), ...
// pnt1, pnt2: variables x, y
// linex2d returns f(x,y)
double linex2d(double *axis1, int len1, double *axis2, int len2, 
               double pnt1, double pnt2, double *data);

// perform linear interpolation of 3D grid
// data tuples (xi, yj, zk, fijk) are given as follows:
// axis1: array of x1, x2, ...
// len1: length of axis1 array
// axis2: array of y1, y2, ...
// len2: length of axis2 array
// axis3: array of z1, z2, ...
// len3: length of axis3 array
// data: array of f(x1,y1,z1), f(x1,y1,z2), ..., f(x1,y1,z_len3), 
//		  f(z1,y2,z1), ..., f(x1,y2,z_len3), 
//		  ...,
//		  f(x1,y_len2,z1), ..., f(x1,y_len2,z_len3), 
//		  f(x2,y1,z1), ...
// pnt1, pnt2, pnt3: variables x, y, z
// linex2d returns f(x,y,z)
double linex3d(double *axis1, int len1, 
               double *axis2, int len2, 
               double *axis3, int len3, 
               double pnt1, double pnt2, double pnt3, double *data);
