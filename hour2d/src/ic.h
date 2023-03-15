#ifndef _IC_
#define _IC_

#include "def2d.h"

void Init( 
	unsigned int irk, int& it0, double& t0, double& cputm0, 
	unsigned int nele, unsigned int nsid, unsigned int nnod, 
	elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa,
	double *res10, double *res20, double *resx0 );


#if		EQUATIO==ADVECTI
double uIC   ( double x, double y, double t, para *ppa );
double uICdt ( double x, double y, double t, para *ppa );
double uICdt2( double x, double y, double t, para *ppa );
#elif	EQUATIO==BURGERS
double uIC   ( double x, double y, double t, para *ppa );
double uICdx ( double x, double y, double t, para *ppa );
double uICdy ( double x, double y, double t, para *ppa );
#elif	EQUATIO==EULNSEQ
double    RIC( double x, double y, double t, para *ppa );
double    mIC( double x, double y, double t, para *ppa );
double    nIC( double x, double y, double t, para *ppa );
double    eIC( double x, double y, double t, para *ppa );
double    PIC( double x, double y, double t, para *ppa );
double    TIC( double x, double y, double t, para *ppa );
double  mICdx( double x, double y, double t, para *ppa );
double  mICdy( double x, double y, double t, para *ppa );
#endif
///////
#endif
