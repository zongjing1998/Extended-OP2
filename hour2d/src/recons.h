#ifndef _RECONSTRUCTION_
#define _RECONSTRUCTION_

#include "def2d.h"

// �ع�wh=(uh,up)
void Recons(
    int irk, int it, double t, double cputm,
    int nele, int nsid, int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa,
    int nbsd, int nisd,
    int *pbsd, int *pisd, cvbd2d *pcvbd);

#endif