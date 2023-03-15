#ifndef _POLY1D_NOMIAL_
#define _POLY1D_NOMIAL_

#include "def2d.h"

//                    use 0-MAXSEGQDPTS :MAXSEGQDPTS
extern double QdPtCoe1d[MAXSEGQDPTS + 1][MAXSEGQDPTS][2];
extern double QdPtCoe1d3[MAXSEGQDPTS + 1][MAXSEGQDPTS][3];
extern double QdPtBas1d[MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVORD + 1];
extern double QdPtBasDe1d[MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVORD + 1];
extern int QdNPt2NOrd1d[MAXSEGQDPTS + 1];
extern int QdNOrd2NPt1d[MAXSEGQDPTS + MAXSEGQDPTS];

//                       use 0-MAXSEGQDPTS :MAXSEGQDPTS*MAXSEGQDPTS
extern double QdPtCoeRect[MAXSEGQDPTS + 1][MAXSEGQDPTS * MAXSEGQDPTS][4];
extern double QdPtBasRect[MAXSEGQDPTS + 1][MAXSEGQDPTS * MAXSEGQDPTS][MAXFVDOF];
extern double QdPtBasDxiRect[MAXSEGQDPTS + 1][MAXSEGQDPTS * MAXSEGQDPTS][MAXFVDOF];
extern double QdPtBasDetRect[MAXSEGQDPTS + 1][MAXSEGQDPTS * MAXSEGQDPTS][MAXFVDOF];
//                            use 1-MAXSDQDPTS  :MAXSEGQDPTS
extern double SdQdPtBasRect[8][MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVDOF];
extern double SdQdPtBasDxiRect[8][MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVDOF];
extern double SdQdPtBasDetRect[8][MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVDOF];
extern double SdQdPtRectXi[8][MAXSEGQDPTS + 1][MAXSEGQDPTS];
extern double SdQdPtRectEt[8][MAXSEGQDPTS + 1][MAXSEGQDPTS];

void CalBasArr1d(para *ppa);
void CalBasArrRect();

void Get1dTensor(
    unsigned int nord,
    double *Bas1d_Xi, double *Basde1d_Xi, double *Bas1d_Et, double *Basde1d_Et,
    double *BasRect, double *BasDxiRect, double *BasDetRect);

void Get1dTensorArbi(
    unsigned int nord, double xi, double et,
    double *BasRect, double *BasDxiRect, double *BasDetRect);

double BasisPoly1d(unsigned int nord, double xi);
double BasisPolyDe1d(unsigned int nord, double xi);

// ndof stand for poly order
inline void GetBas1d(unsigned int ndof, double xi, double *ret)
{
	// double basis[MAXFVDOF];
	for (unsigned int i = 0; i < ndof; i++)
		ret[i] = BasisPoly1d(i, xi);
}

// ndof stand for poly order
inline void GetBasDe1d(unsigned int ndof, double xi, double *ret)
{
	// double basisde[MAXFVDOF];
	for (unsigned int i = 0; i < ndof; i++)
		ret[i] = BasisPolyDe1d(i, xi);
}

#endif