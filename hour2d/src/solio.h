#ifndef _SOLUTIONIO_H_
#define _SOLUTIONIO_H_

#include "def2d.h"

// read the data for restart
void ReadData(
    unsigned int irk, int &it0, double &t0, double &cputm0,
    double *res10, double *res20, double *resx0,
    unsigned int nele, uvar *pvar, para *ppa);
void SaveData(
    unsigned int irk, int it, double t, double cputm,
    double *res10, double *res20, double *resx0,
    unsigned int nele, uvar *pvar, para *ppa);

void OutSol(
    unsigned int irk, int it, double t, double cputm,
    unsigned int nele, unsigned int nsid, unsigned int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa);

void OutWallSol(
    int irk, int it, double t, double cputm,
    int ncvbd, int *pcvbd,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa);

#endif
