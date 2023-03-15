#ifndef _FLUIDDEF_H_
#define _FLUIDDEF_H_

// perfect gas constant
#define gammaa			1.4    	//linux changed
#define gamm1			0.4	//(gamma-1.0)
#define gam1r			2.5 //1.0/gamm1 reciprocal
#define gamrp			3.5 //1.0+1.0/gamm1 == gamma/gamm1

// prandtl's number
#define Pr				0.72
#define Prt				0.9

#define RMIN			0.0000001
#define PMIN			0.0000001
#define RMAX			1000000.0
#define PMAX			1000000.0

#endif
