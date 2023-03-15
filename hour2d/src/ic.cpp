// initial conditions
// 2011 07 28 liming

#include <math.h>
#include <stdio.h>
//#include <conio.h>
#include <memory.h>

#include "ic.h"
#include "misc.h"
#include "solio.h"
#include "poly1d.h"
#include "poly2d.h"
#include "quadrature.h"

#include "m_math_def.h"

#include "omp.h"

#if EQUATIO == ADVECTI
double uICdt(double x, double y, double t, para *ppa)
{
	// return -0.5*PI/20.0*cos( PI*(x+y)/20.0 );
	return -0.5 * PI * cos(PI * (x + y));
}
double uICdt2(double x, double y, double t, para *ppa)
{
	// return -0.5*PI*PI/20.0/20.0*sin( PI*(x+y)/20.0 );
	return -0.5 * PI * PI * sin(PI * (x + y));
}
// exact f(x-at,y-bt)
double uIC(double x, double y, double t, para *ppa)
{
	// return x*x*x*y;
	if (x > 0.0)
		return 1.0 - exp(-750.0 * y / sqrt(x));
	else
		return 1.0;

	return sin(PI * x) * cos(PI2 * y);

	// return 0.25 + 0.5 * sin( PI*(x+y)/20.0 );
	return 0.75 + 0.5 * sin(PI * (x + y));

	if ((x + y) > -1.23)
		return 0.25 + 0.5 * sin(PI * (x + y));
	else
		return 1.25 + 0.5 * sin(PI * (x + y));

	return x * x * y;
}
#elif EQUATIO == BURGERS
double uIC(double x, double y, double t, para *ppa)
{
	// return x*x*x*y;
	// return 0.25 + 0.5 * sin( PI*(x+y)/20.0 );
	return 0.75 + 0.5 * sin(PI * (x + y));

	if ((x + y) > -1.23)
		return 0.25 + 0.5 * sin(PI * (x + y));
	else
		return 1.25 + 0.5 * sin(PI * (x + y));

	return x * x * y;
	return sin(1.0 * PI * x) * cos(3.0 * PI * y);
}
double uICdx(double x, double y, double t, para *ppa)
{
	return 0.5 * PI * cos(PI * (x + y));
}
double uICdy(double x, double y, double t, para *ppa)
{
	return 0.5 * PI * cos(PI * (x + y));
}
#elif EQUATIO == EULNSEQ
double mICdx(double x, double y, double t, para *ppa)
{
	// Couette flow
	return 0.0;
	// return 0.5*PI*cos( PI*(x+y) );
	// return ( 1.0 + 0.2*sin(PI*x)*cos(PI2*y) );
}
double mICdy(double x, double y, double t, para *ppa)
{
	// Couette flow
	return 0.5;
	// return PI*cos(PI*y);
	// return 0.5*PI*cos( PI*(x+y) );
	// return ( 1.0 + 0.2*sin(PI*x)*cos(PI2*y) );
}

double PIC(double x, double y, double t, para *ppa)
{
	//	return RIC( x, y, t, ppa ) * TIC( x, y, t, ppa ) / gamma / ppa->Ma / ppa->Ma;

	return ppa->Preinf;
}

double TIC(double x, double y, double t, para *ppa)
{
	double u_inf = 0.5;
	double u_Ninf = 0.25;
	double deltaw = 1.0;
	double u = 0.5 * (u_inf + u_Ninf + (u_inf - u_Ninf) * tanh(2.0 / deltaw * y));
	// double u = mIC( x, y, t, ppa ) / RIC( x, y, t, ppa );
	double ma2 = ppa->Ma * ppa->Ma;
	// u = (u-0.5*(u1+u2))/ur;
	// return ( 1.0 + 0.5*gamm1*ma2*sqrt(Pr)*(1.0-u*u) );
	// return ( 1.0 + 0.5*gamm1*ma2*sqrt(Pr)*(u_inf-u)*(u-u_Ninf) );
	///	return ( 1.0 + 0.18*(u_inf-u)*(u-u_Ninf) );

	// Couette flow
	///	return ( 0.8 + 0.025*y+ 0.125*y*(2.0-y)*ppa->Tcoe );

	// real ri = 0.5;
	// real ro = 1.0;
	double Ti = 1;

	double r = sqrt(x * x + y * y);
	double tmp = gamm1 * Pr * ppa->Ma * ppa->Ma * 8 / 9.0; // ma=0.1 0.00256, ma=0.3  0.02304
	double tmp1 = 1.3068528194400546905827678785418;       // 2 + log(0.5)

	return (Ti + tmp * (tmp1 - 0.5 / (r * r) - log(r)));
}

double RIC(double x, double y, double t, para *ppa)
{
	//	return PIC( x, y, t, ppa ) / TIC( x, y, t, ppa ) / ppa->R;

	return ppa->rhoinf;
	// return PIC( x, y, t, ppa )/(ppa->R*TIC( x, y, t, ppa ));

	// return 1.0 + 0.2*x*y*x*y;
	// return ( 1.0 + 0.2*sin(PI*x) );	return ( 1.0 + 0.2*x );
	// if( x<0.1 ) return 1.0+0.2*sin(PI*x)*cos(PI2*y);
	// else return 2.0+0.2*sin(PI*x)*cos(PI2*y);
	// return ( 1.0 + 0.2*sin(PI*x)*cos(PI2*y) );
	//  ������
	double xx, yy;
	xx = x;
	while (xx > 10.0) // ���ܲ��ᷢ��
		xx -= 20.0;
	while (xx < -10.0) // ���ܲ��ᷢ��
		xx += 20.0;
	yy = y;
	while (yy > 10.0) // ���ܲ��ᷢ��
		yy -= 20.0;
	while (yy < -10.0) // ���ܲ��ᷢ��
		yy += 20.0;

	double r2 = (xx - 0.0) * (xx - 0.0) + (yy - 0.0) * (yy - 0.0);
	double T = 1.0 - gamm1 / gammaa * 25.0 / 8.0 / PI / PI * pow(2.7182818284590452353602874713527, 1 - r2);
	return pow(T, gam1r);
	// return (40+x);
	// return (x*x+100)/400;
}

double mIC(double x, double y, double t, para *ppa)
{
	// real omega = 1;
	// real ri = 0.5;
	// real ro = 1.0;
	// real C = omega / 3; //(1/ri/ri - 1/ro/ro);

	//	double r = sqrt( x*x + y*y );
	//	double u_theta = 2*(1-r*r) / (3*r);

	//	double tx = -y/r;

	//	return RIC( x, y, t, ppa ) * u_theta * tx;
	// mix layer
	//	double u_inf, u_Ninf;
	//	double deltaw;

	//	u_inf  = 0.5;
	//	u_Ninf = 0.25;
	//	deltaw = 1.0;
	//	return RIC( x, y, t, ppa ) * 0.5*(u_inf+u_Ninf + (u_inf-u_Ninf)*tanh(2.0/deltaw*y));
	// double mach
	// double sq3 = 1.7320508075688772935274463415059;
	// double x0 = 1.0/6.0 + 20.0/sq3 * t;
	// double rho = RIC( x, y, t, ppa );
	// if( x<(x0+y/sq3) )
	//{
	//	return rho*8.25*0.86602540378443864676372317075294;
	//}
	// else
	//	return 0.0;
	// Couette flow
	// return 0.5*y*RIC( x, y, t, ppa );
	return ppa->rhoinf * ppa->uinf;

	// return sin(PI*y);
	// return RIC( x, y, t, ppa );
	//  ������
	double xx, yy;
	xx = x;
	while (xx > 10.0) // ���ܲ��ᷢ��
		xx -= 20.0;
	while (xx < -10.0) // ���ܲ��ᷢ��
		xx += 20.0;
	yy = y;
	while (yy > 10.0) // ���ܲ��ᷢ��
		yy -= 20.0;
	while (yy < -10.0) // ���ܲ��ᷢ��
		yy += 20.0;

	double r2 = (xx - 0.0) * (xx - 0.0) + (yy - 0.0) * (yy - 0.0);
	return RIC(xx, yy, t, ppa) * (0.0 + 5.0 / PI2 * pow(2.7182818284590452353602874713527, 0.5 * (1 - r2)) * (0.0 - yy));
}

double nIC(double x, double y, double t, para *ppa)
{
	// real omega = 1;
	// real ri = 0.5;
	// real ro = 1.0;
	// real C = omega / 3; //(1/ri/ri - 1/ro/ro);

	//	double r = sqrt( x*x + y*y );
	//	double u_theta = 2*(1-r*r) / (3*r);

	// real theta = atan2( y, x );
	//	double ty = x/r;

	//	return RIC( x, y, t, ppa ) * u_theta * ty;
	// mix layer
	//	return 0.0;

	// double mach
	// double sq3 = 1.7320508075688772935274463415059;
	// double x0 = 1.0/6.0 + 20.0/sq3 * t;
	// double rho = RIC( x, y, t, ppa );
	// if( x<(x0+y/sq3) )
	//{
	//	return -rho*8.25*0.5;
	//}
	// else
	//	return 0.0;
	// Couette flow
	// return 0.0;
	return ppa->rhoinf * ppa->vinf;

	// ������
	double xx, yy;
	xx = x;
	while (xx > 10.0) // ���ܲ��ᷢ��
		xx -= 20.0;
	while (xx < -10.0) // ���ܲ��ᷢ��
		xx += 20.0;
	yy = y;
	while (yy > 10.0) // ���ܲ��ᷢ��
		yy -= 20.0;
	while (yy < -10.0) // ���ܲ��ᷢ��
		yy += 20.0;

	double r2 = (xx - 0.0) * (xx - 0.0) + (yy - 0.0) * (yy - 0.0);
	return RIC(xx, yy, t, ppa) * (0.0 + 5.0 / PI2 * pow(2.7182818284590452353602874713527, 0.5 * (1 - r2)) * (xx - 0.0));
}

double eIC(double x, double y, double t, para *ppa)
{
	//	double rho = RIC( x, y, t, ppa );
	//	double Pre = PIC( x, y, t, ppa );
	//	double ru  = mIC( x, y, t, ppa );
	//	double rv  = nIC( x, y, t, ppa );
	//	return Pre*gam1r + 0.5*(ru*ru+rv*rv)/rho;

	// mix layer
	//	double rho = RIC( x, y, t, ppa );
	//	double Pre = PIC( x, y, t, ppa );
	//	double ru  = mIC( x, y, t, ppa );
	//	double rv  = nIC( x, y, t, ppa );
	//	return Pre*gam1r + 0.5*(ru*ru+rv*rv)/rho;

	// double mach
	// double r  = RIC ( x, y, t, ppa );
	// double ru = mIC( x, y, t, ppa );
	// double rv = nIC( x, y, t, ppa );
	// double Pre= PIC ( x, y, t, ppa );
	// return Pre*gam1r + 0.5*(ru*ru+rv*rv)/r;

	// Couette flow
	/*double r  = RIC ( x, y, t, ppa );
	double ru = mIC( x, y, t, ppa );
	double rv = nIC( x, y, t, ppa );
	double Pre= PIC ( x, y, t, ppa );
	return Pre*gam1r + 0.5*(ru*ru+rv*rv)/r;*/

	// return 0.5*y*y/4.0 + ppa->Preinf*gam1r;
	return ppa->einf;

	// return ( 1.0*gam1r + 0.5*RIC( x, y, t, ppa ) );
	//  ������
	double xx, yy;
	xx = x;
	while (xx > 10.0) // ���ܲ��ᷢ��
		xx -= 20.0;
	while (xx < -10.0) // ���ܲ��ᷢ��
		xx += 20.0;
	yy = y;
	while (yy > 10.0) // ���ܲ��ᷢ��
		yy -= 20.0;
	while (yy < -10.0) // ���ܲ��ᷢ��
		yy += 20.0;

	double rho = RIC(xx, yy, t, ppa);
	double Pre = pow(rho, gammaa);
	double ru = mIC(xx, yy, t, ppa);
	double rv = nIC(xx, yy, t, ppa);
	return Pre * gam1r + 0.5 * (ru * ru + rv * rv) / rho;
}
#endif

// initial condition
void Init(
    unsigned int irk, int &it0, double &t0, double &cputm0,
    unsigned int nele, unsigned int nsid, unsigned int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa,
    double *res10, double *res20, double *resx0)
{
	printf("i.c.\n");

	if (ppa->iconti == 1)
	{
		ReadData(irk, it0, t0, cputm0, res10, res20, resx0, nele, pvar, ppa);
		return;
	}

	// physical solution initialization
	it0 = 0;
	t0 = 0.0;
	cputm0 = 0.0;
	memset(res10, 0, sizeof(double) * NEQ);
	memset(res20, 0, sizeof(double) * NEQ);
	memset(resx0, 0, sizeof(double) * NEQ);

// ʹ����������exactic
#pragma omp parallel for
	for (int e = 0; e < (int)nele; e++)
	{
		unsigned int nType = ELEMTYPE;
		int *inode = pele[e].Node;
#ifdef QUAD
		if (inode[3] == inode[0])
			nType = 3;
#endif
		double ndx[ELEMTYPE], ndy[ELEMTYPE];
		for (unsigned int i = 0; i < ELEMTYPE; i++)
		{
			ndx[i] = pnod[inode[i]].xy[0];
			ndy[i] = pnod[inode[i]].xy[1];
		}

		unsigned int pfv = ppa->pFV;
		unsigned int pdg = ppa->pDG;
		unsigned int ndgdof = DOFofOrd2d[pdg];
		pvar[e].pfv = pfv;
		pvar[e].pdg = pdg;

		for (unsigned int i = 0; i < ndgdof; i++)
			for (unsigned int j = 0; j < NEQ; j++)
				pvar[e].uh[irk][i][j] = 0.0;

		unsigned int nqdord = 6; // pfv+pdg;//pfv*2;//pfv+1;//MAXEMQDORD;
		if (nqdord == 0)
			nqdord = 1;
		unsigned int nqdpt;
		double(*EmQdCoe)[4], (*EmQdBas)[MAXFVDOF];
		if (nType == 3) // tri
		{
			nqdpt = SymQdTriDim[nqdord];
			EmQdCoe = (double(*)[4])SymQdTriCoe[nqdord];
			EmQdBas = (double(*)[MAXFVDOF])QdPtBasTri[nqdord];
		}
		else // rect
		{
			nqdpt = QdNOrd2NPt1d[nqdord];
			EmQdCoe = (double(*)[4])QdPtCoeRect[nqdpt];
			EmQdBas = (double(*)[MAXFVDOF])QdPtBasRect[nqdpt];
			// nqdpt = nqdord/2+1;
			nqdpt *= nqdpt;
		}

		for (unsigned int qd = 0; qd < nqdpt; qd++)
		{
			double xi = EmQdCoe[qd][0];
			double et = EmQdCoe[qd][1];
			double coe = EmQdCoe[qd][3];

			double xx, yy;
			if (nType == 3)
				Loc2GlbTri(xx, yy, xi, et, ndx, ndy);
			else
				Loc2GlbRect(xx, yy, xi, et, ndx, ndy);

			double tmcv[NEQ];
#if EQUATIO == ADVECTI
			for (unsigned int j = 0; j < NEQ; j++)
				tmcv[j] = uIC(xx, yy, t0, ppa) * coe;
#elif EQUATIO == BURGERS
			for (unsigned int j = 0; j < NEQ; j++)
				tmcv[j] = uIC(xx, yy, t0, ppa) * coe;
#elif EQUATIO == EULNSEQ
			tmcv[0] = RIC(xx, yy, t0, ppa) * coe;
			tmcv[1] = mIC(xx, yy, t0, ppa) * coe;
			tmcv[2] = nIC(xx, yy, t0, ppa) * coe;
			tmcv[NEQ - 1] = eIC(xx, yy, t0, ppa) * coe;
#endif
			for (unsigned int i = 0; i < ndgdof; i++)
				for (unsigned int j = 0; j < NEQ; j++)
					pvar[e].uh[irk][i][j] += tmcv[j] * EmQdBas[qd][i];
		} // end of elem gauss quad
	}	  // end of elem
}
