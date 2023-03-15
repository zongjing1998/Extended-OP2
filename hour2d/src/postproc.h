#ifndef _POSTFILTEER_
#define _POSTFILTEER_

#include "def2d.h"


void PreciAna( 
	unsigned int irk, int it, double t, double cputm, 
	unsigned int nele, unsigned int nsid, unsigned int nnod, 
	elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa );

void DisDetect( 
	unsigned int irk, int it, double t, double cputm, 
	unsigned int nele, unsigned int nsid, unsigned int nnod, 
	elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa );

int PosiCorreProc( 
	unsigned int irk, int it, double t, double cputm, 
	unsigned int nele, unsigned int nsid, unsigned int nnod, 
	elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa );

#endif