#ifndef _GRIDPARAM_
#define _GRIDPARAM_

#include "def2d.h"
#include <math.h>

// cal side para dl norm check something

void CalSidePara( 
	unsigned int nele, unsigned int nsid, unsigned int nnod, 
	elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa );

void CalElemPara(
	unsigned int nele, unsigned int nsid, unsigned int nnod, 
	elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa );

void PeriodicBC( 
	unsigned int nele, unsigned int nsid, unsigned int nnod, 
	elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa, 
	unsigned int  &nbsd, unsigned int   nisd, 
	unsigned int *&pbsd, unsigned int *&pisd );


// pt0, pt1, pt2 in a anti-clockwise order
inline double area_of_tri( double x0, double y0, double x1, double y1, double x2, double y2 )
{
	return 0.5 * (  ( x1 - x0 ) * ( y2 - y0 ) - ( y1 - y0 ) * ( x2 - x0 )  );
}

// pt0, pt1, pt2 in a anti-clockwise order
inline double area_of_tri( const double * xy0, const double * xy1, const double * xy2 )
{
	return 0.5 * (  ( xy1[0] - xy0[0] ) * ( xy2[1] - xy0[1] ) - ( xy1[1] - xy0[1] ) * ( xy2[0] - xy0[0] )  );
}

// pt0, pt1, pt2 in a anti-clockwise order
inline double area_of_tri( const double * x, const double * y )
{
	return 0.5 * (  ( x[1] - x[0] ) * ( y[2] - y[0] ) - ( y[1] - y[0] ) * ( x[2] - x[0] )  );
}

inline void normalize_vector( const double * xy0, const double * xy1, double & nx, double & ny )
{
	double tnx = xy1[0] - xy0[0];
	double tny = xy1[1] - xy0[1];
	double tmp = sqrt( tnx * tnx + tny * tny );
	nx = tnx / tmp;
	ny = tny / tmp;
}

inline void normalize_vector( const double * xy0, const double * xy1, double * norm )
{
	double tnx = xy1[0] - xy0[0];
	double tny = xy1[1] - xy0[1];
	double tmp = sqrt( tnx * tnx + tny * tny );
	norm[0] = tnx / tmp;
	norm[1] = tny / tmp;
}

inline double normal_of_vector( const double * xy0, const double * xy1 )
{
	double tnx = xy1[0] - xy0[0];
	double tny = xy1[1] - xy0[1];
	return sqrt( tnx * tnx + tny * tny );
}

inline double vector_mul_vector( const double * v0, const double * v1 )
{
	return v0[0] * v1[0] + v0[1] * v1[1];
}

#endif