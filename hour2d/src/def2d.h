#ifndef _DEF_2D_H_
#define _DEF_2D_H_

#include <stdio.h>
//#include <conio.h>

#define PRINTERROR(str)      \
	do                   \
	{                    \
		printf(str); \
		getchar();   \
	} while (0)
#define PRINTERRORN(str)     \
	do                   \
	{                    \
		printf(str); \
		return;      \
	} while (0)
#define PRINTERROR0(str)     \
	do                   \
	{                    \
		printf(str); \
		getchar();   \
		return 0;    \
	} while (0)
#define PRINTERROR1(str)     \
	do                   \
	{                    \
		printf(str); \
		getchar();   \
		return 1;    \
	} while (0)
#define PRINTERRORE(str)     \
	do                   \
	{                    \
		printf(str); \
		getchar();   \
		exit(1);     \
	} while (0)

#include "quadrature.h"

//#define MAX_THREADS		1

// typedef double real;

#include "fluiddef.h"

#define TOL 0.000000000000000001 // ��ֹ���㣬����Ƚϵȴ�ʹ��

//#define TRI			(3)
//#define RECT			(4)
#define QUAD (4)

#ifdef TRI
#define ELEMTYPE TRI
#endif
#ifdef RECT
#define ELEMTYPE RECT
#endif
#ifdef QUAD
#define ELEMTYPE QUAD
#endif

#define NAuxVar (1) // Ϊ�������������������羫ȷ��֮���
#define NDIM (2)    // 2ά

#define MAXDGORD (3)  // 0-4
#define MAXDGDOF (10) // 2ά ���p4�׶���ʽ ((MAXDGORD+1)(MAXDGORD+2))/2
#define MAXFVORD (3)  // 0-4
#define MAXFVDOF (10) // 2ά ���p4�׶���ʽ ((MAXFVORD+1)(MAXFVORD+2))/2

#define IMPLICIT
#define MAXIMDOF 1

#define MAXSEGQDPTS (5) // max segment quadrature points
#define MAXQDORDTRI (6)
#define MAXQDPTSTRI (12)

#define ADVECTI 12
#define BURGERS 23
#define EULNSEQ 34
#define EQUATIO EULNSEQ

#if EQUATIO == ADVECTI
#define NEQ (1) // ��������
#elif EQUATIO == BURGERS
#define NEQ (1) // ��������
#elif EQUATIO == EULNSEQ
#define NEQ (4) //(2+NDIM) //euler����
#define CONSERV 123
#define PRIMITI 234
#define VARIABLE CONSERV
#define AUXVARIABLE CONSERV
#endif

#if EQUATIO == EULNSEQ
const char VAR_CV_STR[NEQ][4] = {"rho", "ru", "rv", "e"};
const char VAR_CV_STR2[NEQ][4] = {"r", "ru", "rv", "e"};
const char VAR_PV_STR[NEQ][4] = {"rho", "u", "v", "pre"};
#else
const char VAR_CV_STR[NEQ][4] = {"u"};
const char VAR_CV_STR2[NEQ][4] = {"u"};
const char VAR_PV_STR[NEQ][4] = {"u"};
#endif

#define SILENCE
//#define	Low_Storage

#ifdef Low_Storage
#define NTM 2 // ʹ�õʹ洢ʱ���ƽ�����ʱ 2
#else
#define NTM 2 // ʹ�÷ǵʹ洢ʱ���ƽ�����ʱ
#endif

#define GLOBTIME 0
#define LOCALEU 1
#define LOCALNS 2

#define UNSTEADY 0
#define STEADY 1
//#define STEADY2		2

typedef struct var2dPp
{
	double uh[NTM][MAXDGDOF][NEQ]; /* r,ru,rv,e  or r,u,v,p coes of polynomials */
	double rhs[MAXDGDOF][NEQ];     /* RHS of each element */

	double wh[MAXFVDOF][NEQ];
	double grad[NDIM][MAXFVDOF][NEQ];

	// double CFL;		/* for time step computation */
	double dt;

	unsigned int pdg; /* DG polynomial order 0-MAXDGORD */
	unsigned int pfv; /* FV reconstruction poly order 0-MAXFVORD */

	// unsigned int isDiscon[NEQ];
	double isDiscon[NEQ];
} uvar;

/*
	element cell face
	side edge line
	node point vertex
*/

const double XiofTri[4] = {0.0, 1.0, 0.0, 0.0};
const double EtofTri[4] = {0.0, 0.0, 1.0, 0.0};
// const unsigned int ndofsdtri[4]={ 0, 1, 2, 0 };// ��Ԫ�� ��s ��Ӧ��[s] [s+1]
// const unsigned int nd2sdtri[4] ={ 1, 2, 0 };
// const unsigned int nxtndtri[4] ={ 1, 2, 0 };
const unsigned int NxtSdTri[4] = {1, 2, 0, 1};
const unsigned int PreSdTri[4] = {2, 0, 1, 2};

const double XiofRect[8] = {-1.0, 1.0, 1.0, -1.0, -1.0};
const double EtofRect[8] = {-1.0, -1.0, 1.0, 1.0, -1.0};
// const unsigned int ndofsdrect[8] = { 0, 1, 2, 3, 0 };// ��Ԫ�� ��s ��Ӧ��[s] [s+1]
// const unsigned int nd2sdrect[4]  = { 1, 2, 3, 0 };
// const unsigned int nxtndrect[4]  = { 1, 2, 3, 0 };
const unsigned int NxtSdRect[4] = {1, 2, 3, 0};
const unsigned int PreSdRect[4] = {3, 0, 1, 2};

typedef struct elem2d
{
	int Node[ELEMTYPE];
	unsigned int Side[ELEMTYPE];
	// int nnodenum;
	// int Elem[ELEMTYPE];

	// unsigned int material;	// reserved

	// double xv, yv, xin, yin, R, r, Det;

	// unsigned int eleType;		/* for 2D 3 or 4 */
	double xyc[NDIM];
	double volr; /* ������ ����� ������ ��� */

	// #ifdef RECT
	//	double LocdGlb[0][1], LocdGlb[1][1];		/* xix, ety */
	// #else
	//	double LocdGlb[0][NDIM], LocdGlb[1][NDIM];  /* xix,xiy, etx,ety */
	// #endif

	double LocdGlb[NDIM][NDIM]; /* xidx, xidy, etdx, etdy */
	double dh;
	double dxy[NDIM];
} elem;

class basic_element;

#define LEFT__ELEM 1
#define RIGHT_ELEM 0

typedef struct side2d
{
	unsigned int Node[2];
	int ofElem[2]; // 1: left elem; 0: right element  �����߽����� 0 < 1

	unsigned int lsdidx[2]; // local side index of 1: left elem; 0: right elem

	// int bcType;
	double dl;	   // length
	double norm[NDIM]; //

	double etaf;	    // stabilization parameter for BR2
	unsigned int nqdpt; // number of quadrature points

	double flux[5][NEQ]; // double (*flux)[NEQ];MAXSDQDPTS
} side;

typedef struct node2d
{
	double xy[NDIM];

	// int bcType;

	unsigned int nNSide; // ����������Ŀ
	unsigned int nNElem; // ��������Ԫ��Ŀ

	unsigned int *ofSide; // ���������б�
	unsigned int *ofElem; // ��������Ԫ�б�
	double *ofElemWht;    // ��������ԪȨֵ
} node;

typedef struct grid2d
{
	unsigned int nele;
	unsigned int nsid;
	unsigned int nnod;
	unsigned int nbosid;
	unsigned int ninsid;
	unsigned int ncvbd;

	elem *pele;
	side *psid;
	node *pnod;
	uvar *pvar;
	unsigned int *pbosid;
	unsigned int *pinsid;
	// cvbd2d *pcvbd;
} grid;

typedef struct param2d
{
	int ivis;      // 0: Euler		1: Navier-Stokes
	int fluxtype;  // Euler fluxes 2: Roe 3:
	int pDG;       // polyord;		DG // 0-3
	int pFV;       // polyrecons;	FV //
	int isteady;   // 0: unsteady (global time step) 1: steady
	int itmax;     // number of time step
	int wffre;     // frequence for the solution to be saved
	int wsfre;     // frequence for the solution to be printed
	int dtfre;     // frequence for the dt to be determined
	int iconti;    // 0: start with uniform solution 1: from INIT_
	int resi;      // order of magnitude for the residual to be reduced (for steady problem)
	int CFLIncStp; // cfl linearly increase from step 1 -> this
	// runge-kutta
	int iRuKu;
	int bcTtype;	 // Temperature B.C. 1: isothermal wall 2: adiabatic wall
	int ImResSmoStp; // Implicit Residual Smoothing Steps

	double ImResSmoCoe; // Implicit Residual Smoothing Coefficient
	int ChgMaxStp;
	double M_TVB;	    // TVB_shu M_TVB
	double ChgMaxReduc; //
	double ChgMaxRatio; // rho pre change max fraction per step
	double Re;	    // Reynolds number
	double Ma;	    // inflow mach number
	double RatioP;	    // ratio Pout/Pin
	double Tinfd;	    // inflow temperature for Sutherland laws (dimensional)
	double Twald;	    // if isothermal walls, wall temperature (dimensional)
	double AoA;	    // angle of attack
	double CFL;	    // CFL number
	double tmax;	    // maximum physical time for run (for unsteady problem)

	char gfile[256]; // grid file name

	double a;
	double b;

#if EQUATIO == BURGERS
	// double (*ndcoordt)[NDIM];// burgers��ȷ���һ��������������������������
#endif

	double rhoinf, uinf, vinf, Preinf, Tinf, einf;
	double cosAoA, sinAoA;
	double Ttinf, Pretinf, Sinfn1, Rhoq2d2n1;

	double Twal;	// if isothermal walls, wall temperature (nondimensional)
	double c4suth;	// 110.4/Tinfd
	double R;	// P=rho R T
	double Tcoe;	// ==gamm1/R, T=(gamm1/R*(e/rho-0.5*q*q))
	double KapdMiu; // miu*KapdMiu*Tx KapdMiu = 1.0/Pr/gamm1/ma/ma;

	// double rhoout, uout, vout, Preout, Tout, eout;
	// double xkinf, xeinf;
} para;

typedef struct
{

} cvbd2d;

#endif
