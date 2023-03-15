#ifndef _MISCE_H_
#define _MISCE_H_

#include "def2d.h"

#include <math.h>

// ndof stand for poly dof degree of freedom cij * bi
inline void GetVar(unsigned int ndof, const double (*coe)[NEQ],
		   const double *basis, double *ret)
{
	// for (unsigned int j = 0; j < NEQ; j++)
	// {
	// 	double tret = 0.0;

	// 	for (unsigned int i = 0; i < ndof; i++)
	// 		tret += coe[i][j] * basis[i];

	// 	ret[j] = tret;
	// }
	double tret[NEQ] = {0.0};
	for (int i = 0; i < ndof; i++)
		for (int j = 0; j < NEQ; j++)
		{
			tret[j] += coe[i][j] * basis[i];
			// ret[j] = tret;
		}
	for (int j = 0; j < NEQ; j++)
	{
		ret[j] = tret[j];
		// ret[j] = tret;
	}
}

void GetValuede2d_aban(unsigned int nord, const double (*pcoe)[NEQ],
		       unsigned int twh, double xi, double et, double *ret);

double CalDt(unsigned int irk, double t, unsigned int nele, elem *pele, uvar *pvar, para *ppa);

int CalPre(const double *cv, const double *cv_0, double &pre);

// transfer conser vars to euler flux 2d
inline int CV2FluxI(const double *cv, const double *cv0, double (*fluxI)[NEQ])
{
	double pre;
	int ret = CalPre(cv, cv0, pre);

	const double *pcv = cv;
	if (ret == 1)
		pcv = cv0;

	double u = pcv[1] / pcv[0];
	double v = pcv[2] / pcv[0];
	double e = pcv[NEQ - 1];

	double ep = e + pre;

	fluxI[0][0] = pcv[1]; // rho u
	fluxI[0][1] = fluxI[0][0] * u + pre;
	fluxI[0][2] = fluxI[0][0] * v;
	fluxI[0][NEQ - 1] = ep * u;

	fluxI[1][0] = pcv[2];	   // rho v
	fluxI[1][1] = fluxI[0][2]; // fluxI[1][0]*u;
	fluxI[1][2] = fluxI[1][0] * v + pre;
	fluxI[1][NEQ - 1] = ep * v;

	return ret;
}

// transfer primitive vars to conservation vars
// ʹ��NEQ��Ϊ���޸ĳ�����ά��һά�ȷ���
inline void PV2CV(const double *pv, double *cv)
{
	cv[0] = pv[0];	       // rho
	cv[1] = pv[0] * pv[1]; // rho u
	cv[2] = pv[0] * pv[2]; // rho v
	cv[NEQ - 1] = pv[NEQ - 1] * gam1r + 0.5 * pv[0] * (pv[1] * pv[1] + pv[2] * pv[2]);
	// pv[NEQ-1]*gam1r + 0.5*pv[0]*(pv[1]*pv[1]+pv[2]*pv[2]);
}

// transfer conservation vars to primitive vars
inline int CV2PV(const double *cv, const double *cv0, double *pv)
{
	double pre;
	int ret = CalPre(cv, cv0, pre);

	const double *pcv = cv;
	if (ret == 1)
		pcv = cv0;

	pv[0] = pcv[0];
	pv[1] = pcv[1] / pcv[0];
	pv[2] = pcv[2] / pcv[0];
	pv[NEQ - 1] = pre;

	return ret;
}

inline void PV2FluxI(const double *pv, double (*fluxI)[NEQ])
{
	double rho = pv[0];
	double u = pv[1];
	double v = pv[2];
	double pre = pv[NEQ - 1];
	double e = pre * gam1r + 0.5 * rho * (u * u + v * v); // pre/gamm1 + 0.5*rho*(u*u+v*v);
	double ep = e + pre;

	fluxI[0][0] = rho * u; // rho u
	fluxI[0][1] = fluxI[0][0] * u + pre;
	fluxI[0][2] = fluxI[0][0] * v;
	fluxI[0][NEQ - 1] = ep * u;

	fluxI[1][0] = rho * v;	   // rho v
	fluxI[1][1] = fluxI[0][2]; // fluxI[1][0]*u;
	fluxI[1][2] = fluxI[1][0] * v + pre;
	fluxI[1][NEQ - 1] = ep * v;
}

inline double CalMu(double Tdimles, double c4suth)
{
	// currently for Couette flow
	return 1.0;
	// mix layer
	return Tdimles;

	double dnum = 1.0 + c4suth;	// numerator����
	double dden = Tdimles + c4suth; // denominator��ĸ
	return Tdimles * sqrt(Tdimles) * dnum / dden;
}

// use ux, uy, vx, vy, miu
inline void CalTau(double (*gradpv)[NEQ], double miu, double (*tau)[NDIM])
{
	double miu2d3 = miu * 0.66666666666666666666666666666667;
	tau[0][0] = (2 * gradpv[0][1] - gradpv[1][2]) * miu2d3;
	tau[0][1] = (gradpv[1][1] + gradpv[0][2]) * miu;

	tau[1][0] = tau[0][1]; // (    gradpv[1][1]+gradpv[0][2]) * miu;
	tau[1][1] = (2 * gradpv[1][2] - gradpv[0][1]) * miu2d3;
}

// ��֪�غ����cv���䵼��cvx cvy����ԭʼ����rho,u,v,T��x,yƫ����
inline void CV2GradPV(
    const double *cv, double (*gradcv)[NEQ], double Tcoe, double (*gradpv)[NEQ])
{
	double r = cv[0];
	double rr1 = 1.0 / r;
	double u = cv[1] * rr1;
	double v = cv[2] * rr1;
	double E = cv[3] * rr1;

	gradpv[0][0] = gradcv[0][0];
	gradpv[0][1] = (gradcv[0][1] - gradcv[0][0] * u) * rr1;
	gradpv[0][2] = (gradcv[0][2] - gradcv[0][0] * v) * rr1;
	gradpv[0][NEQ - 1] = (gradcv[0][NEQ - 1] - gradcv[0][0] * E) * rr1 - u * gradpv[0][1] - v * gradpv[0][2];
	gradpv[0][NEQ - 1] *= Tcoe; //(gamm1) * gamma *Ma*Ma;

	gradpv[1][0] = gradcv[1][0];
	gradpv[1][1] = (gradcv[1][1] - gradcv[1][0] * u) * rr1;
	gradpv[1][2] = (gradcv[1][2] - gradcv[1][0] * v) * rr1;
	gradpv[1][NEQ - 1] = (gradcv[1][NEQ - 1] - gradcv[1][0] * E) * rr1 - u * gradpv[1][1] - v * gradpv[1][2];
	gradpv[1][NEQ - 1] *= Tcoe; //(gamm1) * gamma *Ma*Ma;
}

// use u, v, tau, Tx, Ty, miu ,kappa
inline void PV2FluxV(
    double *pv, double (*gradpv)[NEQ], double miu, double kap, double (*fluxV)[NEQ])
{
	double u = pv[1];
	double v = pv[2];
	double tau[NDIM][NDIM];
	CalTau(gradpv, miu, tau);

	fluxV[0][0] = 0.0;
	fluxV[0][1] = tau[0][0];
	fluxV[0][2] = tau[0][1];
	fluxV[0][NEQ - 1] = u * fluxV[0][1] + v * fluxV[0][2] + kap * gradpv[0][NEQ - 1];

	fluxV[1][0] = 0.0;
	fluxV[1][1] = tau[1][0];
	fluxV[1][2] = tau[1][1];
	fluxV[1][NEQ - 1] = u * fluxV[1][1] + v * fluxV[1][2] + kap * gradpv[1][NEQ - 1];
}

inline void CalFluxV(
    const double *uh, const double *uh0, double (*graduh)[NEQ], para *ppa,
    double (*fluxV)[NEQ])
{
#if EQUATIO == EULNSEQ
// ����ʹ���غ����ʱ
#if VARIABLE == CONSERV
	double pv[NEQ];
	int ret = CV2PV(uh, uh0, pv);
	const double *cv = uh;
	if (ret == 1)
		cv = uh0;
#if AUXVARIABLE == CONSERV // ��������ʹ���غ����ʱ
	double gradpv[NDIM][NEQ];
	CV2GradPV(cv, graduh, ppa->Tcoe, gradpv);
#elif AUXVARIABLE == PRIMITI // ��������ʹ��ԭʼ����ʱ
	double(&gradpv)[NDIM][NEQ] = graduh;
#endif
// ����ʹ��ԭʼ����ʱ����ʱ��������ʹ��ԭʼ����
#elif VARIABLE == PRIMITI
	const double *pv = uh;
	double(&gradpv)[NDIM][NEQ] = graduh;
#endif

	double T = pv[NEQ - 1] / (pv[0] * ppa->R);
	double miu = ppa->Re * CalMu(T, ppa->c4suth);

	PV2FluxV(pv, gradpv, miu, miu * ppa->KapdMiu, fluxV);
#endif
}

// euler flux at norm dirction!
inline void CV2FluxIN2d(const double *cv, const double *cv0, double *norm, double *FluxIN)
{
#if EQUATIO == EULNSEQ
	double e;
#if VARIABLE == CONSERV
	double pv[NEQ];
	int ret = CV2PV(cv, cv0, pv);
	e = cv[NEQ - 1];
	if (ret == 1)
		e = cv0[NEQ - 1];
#elif VARIABLE == PRIMITI
	const double *pv = cv;
	e = pv[NEQ - 1] * gam1r + 0.5 * pv[0] * (pv[1] * pv[1] + pv[2] * pv[2]); // pv[NEQ-1]/gamm1 + 0.5*pv[0]*(pv[1]*pv[1]+pv[2]*pv[2]);
#endif
	double un = norm[0] * pv[1] + norm[1] * pv[2];
	FluxIN[0] = pv[0] * un;
	FluxIN[1] = pv[1] * FluxIN[0] + norm[0] * pv[NEQ - 1];
	FluxIN[2] = pv[2] * FluxIN[0] + norm[1] * pv[NEQ - 1];
	FluxIN[NEQ - 1] = (pv[NEQ - 1] + e) * un;
#endif
}

inline void Glb2LocTri(
    double &xi, double &et, const double *xy, const double *xy0,
    double xix, double xiy, double etx, double ety)
{
	xi = xix * (xy[0] - xy0[0]) + xiy * (xy[1] - xy0[1]);
	et = etx * (xy[0] - xy0[0]) + ety * (xy[1] - xy0[1]);
}

inline void Loc2GlbTri(
    double &x, double &y, double xi, double et,
    const double *ndx, const double *ndy)
{
	x = (ndx[1] - ndx[0]) * xi + (ndx[2] - ndx[0]) * et + ndx[0];
	y = (ndy[1] - ndy[0]) * xi + (ndy[2] - ndy[0]) * et + ndy[0];
}

inline void Glb2LocRect(
    double &xi, double &et, const double *xy, const double *xyc,
    double xix, double xiy, double etx, double ety)
{
	xi = xix * (xy[0] - xyc[0]) + xiy * (xy[1] - xyc[1]);
	et = etx * (xy[0] - xyc[0]) + ety * (xy[1] - xyc[1]);
}

inline void Loc2GlbRect(
    double &x, double &y, double xi, double et,
    const double *ndx, const double *ndy)
{
	x = 0.5 * ((ndx[1] - ndx[0]) * xi + (ndx[3] - ndx[0]) * et + ndx[1] + ndx[3]);
	y = 0.5 * ((ndy[1] - ndy[0]) * xi + (ndy[3] - ndy[0]) * et + ndy[1] + ndy[3]);
}

inline void SetRecoOrd(uvar *pvar, unsigned int e, unsigned int nfvord)
{
	pvar[e].pfv = nfvord;
}

inline void SetElemOrd(uvar *pvar, unsigned int e, unsigned int ndgord)
{
	pvar[e].pdg = ndgord;
}

inline void SetSideQdpts(side *psid, unsigned int s, unsigned int nqdpts)
{
	psid[s].nqdpt = nqdpts;
}

const int OrdofDOF2d[36] =
    {
	0,
	1, 1,
	2, 2, 2,
	3, 3, 3, 3,
	4, 4, 4, 4, 4,
	5, 5, 5, 5, 5, 5,
	6, 6, 6, 6, 6, 6, 6,
	7, 7, 7, 7, 7, 7, 7, 7};

const int DOFofOrd2d[16] =
    {
	1, 3, 6, 10,
	15, 21, 28, 36,
	45, 55, 66, 78,
	91, 105, 120, 136};

#endif