#ifndef _M_MATH_DEF_H_
#define _M_MATH_DEF_H_

#include <stdio.h>
#include <stdlib.h>

#ifndef m_abs
#define m_abs(a)	(((a) > 0) ? (a) : (-(a)))
#endif

#ifndef m_max
#define m_max(a,b)	(((a) > (b)) ? (a) : (b))
#endif

#ifndef m_min
#define m_min(a,b)  (((a) < (b)) ? (a) : (b))
#endif

#ifndef m_sgn
#define m_sgn(a)	( ((a) > 0.0) ? (1.0) : (-1.0) )
#endif

#ifndef m_SIGN
#define m_SIGN(a,b) ( ((b) >= 0.0) ? (m_abs(a)) : (-m_abs(a)) )
#endif

#ifndef m_max_3
#define m_max_3(a,b,c)  ( m_max( m_max(a, b), c ) )
#endif

#ifndef m_min_3
#define m_min_3(a,b,c)  ( m_min( m_min(a, b), c ) )
#endif

#ifndef PI
#define PI		3.1415926535897932384626433832795
#endif

#ifndef PI2
#define PI2		6.283185307179586476925286766559
#endif

#ifndef PI_180
#define PI_180	0.017453292519943295769236907684886
#endif

#ifndef RAD2DEG
#define RAD2DEG	57.295779513082320876798154814105
#endif

/* Define _USE_MATH_DEFINES before including math.h to expose these macro
 * definitions for common math constants.  These are placed under an #ifdef
 * since these commonly-defined names are not part of the C/C++ standards.
 */

/* Definitions of useful mathematical constants
 * M_E        - e
 * M_LOG2E    - log2(e)
 * M_LOG10E   - log10(e)
 * M_LN2      - ln(2)
 * M_LN10     - ln(10)
 * M_PI       - pi
 * M_PI_2     - pi/2
 * M_PI_4     - pi/4
 * M_1_PI     - 1/pi
 * M_2_PI     - 2/pi
 * M_2_SQRTPI - 2/sqrt(pi)
 * M_SQRT2    - sqrt(2)
 * M_SQRT1_2  - 1/sqrt(2)
 */

#define M_E        2.71828182845904523536
#define M_LOG2E    1.44269504088896340736
#define M_LOG10E   0.434294481903251827651
#define M_LN2      0.693147180559945309417
#define M_LN10     2.30258509299404568402
#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#define M_PI_4     0.785398163397448309616
#define M_1_PI     0.318309886183790671538
#define M_2_PI     0.636619772367581343076
#define M_2_SQRTPI 1.12837916709551257390
#define M_SQRT2    1.41421356237309504880
#define M_SQRT1_2  0.707106781186547524401


inline double m_dsqr( double a ) { return a*a; }

#include <math.h>

// Computes (a^2 + b^2)^(1/2) without destructive underflow or overflow.
inline double m_norm( double a, double b )
{
	double absa = m_abs(a);
	double absb = m_abs(b);
	if( absa>absb ) return absa*sqrt( 1.0+m_dsqr(absb/absa) );
	else return ( absb==0.0 ? 0.0 : absb*sqrt( 1.0+m_dsqr(absa/absb) ) );
}

//static double dswaptemp;
//#define DSWAP(a,b) {dswaptemp=(a);(a)=(b);(b)=dswaptemp;}
inline void m_dswap( double & a, double & b )
{
	double dtmp = a;
	a = b;
	b = dtmp;
}

/* standard error handler */
inline void m_error( char * error_text )
{
	//fprintf(stderr,"Run-time error...\n");
	fprintf( stderr, "%s\n", error_text );
	//fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


#endif