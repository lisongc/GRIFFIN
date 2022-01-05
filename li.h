/* ------------------------------------------------------------- */
/* li.h, li.cc                                                   */
/* Stefan Bauberger (stefan@bauberger.de)                        */
/* last revision 4.12.95                                         */
/* ------------------------------------------------------------- */
/* Calculation of poly-logarithms Li_2 and Li_3 for real or      */
/* complex values with simple taylor series                      */
/* adapted for use with double or quad-precision                 */ 
/* ------------------------------------------------------------- */

// Li_2: li2(Cplx x), li2(Double x)=li2(x-I*eps), li2c(Double x)=li2(x+I*eps)
// Li_3: li3(Cplx x), li3(Double x)=li3(x-I*eps), li3c(Double x)=li3(x+I*eps)

#ifndef __li__
#define __li__

#include "Cplx.h"

Cplx li2(Double x);
Cplx li2c(Double x);
Cplx li2(Cplx x);
Cplx li3(Double x);
Cplx li3c(Double x);
Cplx li3(Cplx x);

#endif

