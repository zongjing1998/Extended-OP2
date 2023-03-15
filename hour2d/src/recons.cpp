// _RECONSTRUCTION_

#include <stdio.h>
//#include <conio.h>
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>

#include "bc2d.h"
#include "misc.h"
#include "poly1d.h"
#include "poly2d.h"
#include "recons.h"
#include "quadrature.h"

#include "omp.h"

// �ع�wh=(uh,up)
void Recons(
    int irk, int it, double t, double cputm,
    int nele, int nsid, int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa,
    int nbsd, int nisd,
    int *pbsd, int *pisd, cvbd2d *pcvbd)
{
#pragma omp parallel for
	for (int e = 0; e < (int)nele; e++)
	{
		unsigned int pdg = pvar[e].pdg;
		unsigned int pfv = pvar[e].pfv;
		unsigned int ndgdof = DOFofOrd2d[pdg];
		unsigned int nfvdof = DOFofOrd2d[pfv];

		for (unsigned int i = 0; i < ndgdof; i++)
			for (unsigned int j = 0; j < NEQ; j++)
				pvar[e].wh[i][j] = pvar[e].uh[irk][i][j];

		for (unsigned int i = ndgdof; i < nfvdof; i++)
			for (unsigned int j = 0; j < NEQ; j++)
				pvar[e].wh[i][j] = 0.0;
	}
}
