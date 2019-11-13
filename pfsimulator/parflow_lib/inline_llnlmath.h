#ifndef _INLINE_LLNLMATH_H
#define _INLINE_LLNLMATH_H

/* Clone of llnlmath.c in header form for easier build inlining */

#include "llnltyps.h"
#include <math.h>

BEGIN_EXTERN_C

#ifndef ZERO
#define ZERO RCONST(0.0)
#endif

#ifndef ONE
#define ONE RCONST(1.0)
#endif

#ifndef TWO
#define TWO RCONST(2.0)
#endif

/******************************************************************
*                                                                *
* Macros : MIN, MAX, ABS, SQR                                    *
*----------------------------------------------------------------*
* MIN(A, B) returns the minimum of A and B.                      *
*                                                                *
* MAX(A, B) returns the maximum of A and B.                      *
*                                                                *
* ABS(A) returns the absolute value of A.                        *
*                                                                *
* SQR(A) returns the square of A.                                *
*                                                                *
******************************************************************/
#ifndef MIN
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef ABS
#define ABS(A)    ((A < 0) ? -(A) : (A))
#endif

#ifndef SQR
#define SQR(A)    ((A)*(A))
#endif

/******************************************************************
*                                                                *
* Function : UnitRoundoff                                        *
* Usage    : real uround;                                        *
*            uround = UnitRoundoff();                            *
*----------------------------------------------------------------*
* UnitRoundoff returns the unit roundoff u for real floating     *
* point arithmetic, where u is defined to be the smallest        *
* positive real such that 1.0 + u != 1.0.                        *
*                                                                *
******************************************************************/

inline
real UnitRoundoff(void)
{
  real u;
  volatile real one_plus_u;
  u = ONE;
  one_plus_u = ONE + u;
  while (one_plus_u != ONE)
  {
    u /= TWO;
    one_plus_u = ONE + u;
  }
  u *= TWO;
  return (u);
}


/******************************************************************
*                                                                *
* Function : RPowerI                                             *
* Usage    : int exponent;                                       *
*            real base, ans;                                     *
*            ans = RPowerI(base,exponent);                       *
*----------------------------------------------------------------*
* RPowerI returns the value base^exponent, where base is a real  *
* and exponent is an int.                                        *
*                                                                *
******************************************************************/

inline
real RPowerI(real base, int exponent)
{
  int i, expt;
  real prod;

  prod = ONE;
  expt = ABS(exponent);
  for (i = 1; i <= expt; i++)
    prod *= base;
  if (exponent < 0)
    prod = ONE / prod;
  return(prod);
}


/******************************************************************
*                                                                *
* Function : RPowerR                                             *
* Usage    : real base, exponent, ans;                           *
*            ans = RPowerR(base,exponent);                       *
*----------------------------------------------------------------*
* RPowerR returns the value base^exponent, where both base and   *
* exponent are reals. If base < 0.0, then RPowerR returns 0.0.   *
*                                                                *
******************************************************************/

inline
real RPowerR(real base, real exponent)
{
  if (base <= ZERO)
    return(ZERO);

  return((real)pow((double)base, (double)exponent));
}


/******************************************************************
*                                                                *
* Function : RSqrt                                               *
* Usage    : real sqrt_x;                                        *
*            sqrt_x = RSqrt(x);                                  *
*----------------------------------------------------------------*
* RSqrt(x) returns the square root of x. If x < 0.0, then RSqrt  *
* returns 0.0.                                                   *
*                                                                *
******************************************************************/

inline
real RSqrt(real x)
{
  if (x <= ZERO)
    return(ZERO);

  return((real)sqrt((double)x));
}

END_EXTERN_C

#endif // _INLINE_LLNLMATH_H
