#ifndef _GRIDIO_H_
#define _GRIDIO_H_

#include "gridparam.h"
#include "def2d.h"
#include "rcm.h"
#define TETREXTETRA 1 // HEDRON
#define TETREXHEXA 2  // HEDRON
#define TETREXPRISM 3
#define TETREXPYRAMID 4
#define TETREXQUAD 5 // RILATERAL
#define TETREXTRI 6  // ANGLE

int *ReadGridTetrex(
    const char *sGridFile,
    int &nele, int &nsid, int &nnod,
    elem *&pele, side *&psid, node *&pnod, uvar *&pvar,
    int &nbsd, int &nisd, int &ncvbd,
    int *&pbsd, int *&pisd, cvbd2d *&pcvbd, int rcm);

// ����߽�����
void ReadBCTetrex(
    FILE *sgrid,
    int nele, int nsid, int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar,
    int nbsd, int nisd, int ncvbd,
    int *&pbsd, int *&pisd, cvbd2d *pcvbd, int *elemind, int rcm);

// �Ե������߼���Ԫ����������ʱ�뷽��
// ���ڱ߽�㣬����ֻ����nNSide=nNElem+1
void SortListofNd(
    int nele, int nsid, int nnod,
    elem *pele, side *psid, node *pnod);

// return 1 error // return 0 ok!
int CheckGrid(
    int nele, int nsid, int nnod,
    const elem *pele, const side *psid, const node *pnod);

#endif
