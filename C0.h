/* ------------------------------------------------------------- */
/* C0.h, C0.cc                                                   */
/* Ayres Freitas (afreitas@pitt.edu)			         */
/* last revision 06.08.19                                        */
/* ------------------------------------------------------------- */
/* scalar one-loop function C0                     	         */
/* ------------------------------------------------------------- */

#include "s2lseinline.h"
	
// the one-loop scalar triangle function
Cplx C0Q(Double p10, Double p20, Double p12, Double m1, Double m2, Double m3);

// the derivative of C0 with respect to the third (squared) mass
Cplx DM3C0(Double p10, Double p20, Double p12, Double m1, Double m2, Double m3);
