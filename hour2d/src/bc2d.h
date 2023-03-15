#ifndef _BC2D_H_
#define _BC2D_H_

#include "def2d.h"

#define ISOTWALL 1
#define ADIAWALL 2

#define BC_BLKCONNECTION (-1)  // block connection	С��0Ϊ�����ӱ߽�����
#define BC_NONE 0x0000	       // None------------------------------->Tetrex bc 0 Unspecified
#define BC_INTERBLK 0x0001     // interblock connection���ڲ���������>Tetrex bc 6 BndCond_6
#define BC_SOLIDSURFACE 0x0002 // solid surface �̱�----------------->Tetrex bc 3 Wall
#define BC_SYMMETRY 0x0003     // symmetry	�Գ�---------------------->Tetrex bc 4 Symmetry
#define BC_FARFIELD 0x0004     // farfield	Զ��---------------------->Tetrex bc 7 BndCond_7
#define BC_INFLOW 0x0005       // inflow	����---------------------->Tetrex bc 1 Inflow
#define BC_OUTFLOW 0x0006      // outflow  ����---------------------->Tetrex bc 2 Outflow
#define BC_I_POLE 000071       // i pole	I����
#define BC_J_POLE 000072       // j pole	J����
#define BC_K_POLE 000073       // k pole	K����
#define BC_GENERIC1 0x0008     // generic #1	������Ҫ����----------->Tetrex bc 8 BndCond_8
#define BC_GENERIC2 0x0009     // generic #2	������Ҫ����
#define BC_GENERIC3 0x000A     // generic #3	������Ҫ����
#define BC_OTHERS 0x0010       // �������Σ������������ж��塣

#define N_BCTYPE 9
const int nType1 = 4;

const unsigned int nTetrexBC[N_BCTYPE] =
    {
	BC_NONE,
	BC_INTERBLK,
	BC_SOLIDSURFACE,
	BC_SYMMETRY,
	BC_FARFIELD,
	BC_INFLOW,
	BC_OUTFLOW,
	0,
	BC_GENERIC1};

const char strTetrexBC[N_BCTYPE][12] =
    {
	"Unspecified",
	"BndCond_6",
	"Wall",
	"Symmetry",
	"BndCond_7",
	"Inflow",
	"Outflow",
	"POLEnotUsed",
	"BndCond_8"};

void bcSolidSurf(const double *cv, const double *cv0, const double *norm, double *ret);
void bcvisSolidSurf(const double *cv, const double *cv0, para *ppa, double *ret);

void bcFarField(
    const double *uh, const double *uh0, const double *norm, para *ppa, double *ret);

// �ԳƱ߽�����
void bcSymmetry(
    const double *uh, const double *uh0, const double *norm, double *ret);
// �ԳƱ߽�����
void bcSymmetry2(
    const double *uh, const double *uh0, const double *norm, double *ret);

#define GetElement(e0) ((unsigned int)(((unsigned int)(-(e0))) >> 4))
#define GetBCType(e0) ((unsigned int)(((unsigned int)(-(e0))) & 0xF))
#define SetBCType(e0, bctype) (-((int)((((unsigned int)(e0)) << 4) + ((unsigned int)(bctype)))))

#endif
