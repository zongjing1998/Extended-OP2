#ifndef _LIMITER_H_
#define _LIMITER_H_

#include "m_math_def.h"

inline double minmod3( double a, double b, double c )
{
	if( ((a*b)>0.0)&&((b*c)>0.0) )
		return m_sgn(a) * m_min_3( m_abs(a), m_abs(b), m_abs(c) );
	else
		return 0.0;
}

inline double minmod2( double a, double b )
{
	//return (a*b) < 0 ? 0.0 : ( fabs(a) < fabs(b) ? a : b);
	
	if( (a*b)<0.0 )	return 0.0;
	else if( m_abs(a)<m_abs(b) ) return a;
	else return b;
}

#define epsi_vande	0.00001
inline double vande( double a, double b ,double c )
{
	if( ((a*b)>0)||((a*c)>0.0) )
		return m_sgn(a)*m_min_3( m_abs(a), m_abs(b), 2.0*b*c/(m_abs(a)+m_abs(c)+epsi_vande) );
	else
		return 0.0;
}
#undef epsi_vande


// limiter
inline double limiter2( double a, double b )
{
	return minmod2( a, b );
}

inline double limiter3( double a, double b, double c )
{
	return minmod3( a, b, c );
}

#endif
