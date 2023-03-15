// bc for 2d bc2d.cpp
// 2011 11 14 rewrite liming

#include "bc2d.h"
#include "misc.h"
#include "def2d.h"
#include "flux2d.h"
#include "quadrature.h"

#include <math.h>
//#include <conio.h>
#include <stdio.h>
#include <memory.h>



// 远场黎曼不变量边界条件
void bcFarField( 
	const double* uh, const double* uh0, const double* norm, para *ppa, double* ret)
{
#if		VARIABLE==CONSERV
	const double *cv = uh;
	double pre;
	if( CalPre( uh, uh0, pre ) )
		cv = uh0;

	double rho = cv[0];
	double u   = cv[1]/rho;
	double v   = cv[2]/rho;
#elif	VARIABLE==PRIMITI
	double rho = uh[0];
	double u   = uh[1];
	double v   = uh[2];
	double pre = uh[NEQ-1];
#endif

	if(rho<RMIN) { rho = RMIN; printf("far rho<0.0\n" ); getchar(); }
	if(rho>RMAX) { rho = RMAX; printf("far rho>max\n" ); getchar(); }
	if(pre<PMIN) { pre = PMIN; printf("far pre<0.0\n" ); getchar(); }
	if(pre>PMAX) { pre = PMAX; printf("far pre>max\n" ); getchar(); }

	double nx = norm[0];
	double ny = norm[1];

	double un = u*nx + v*ny;
	double a  = sqrt( gammaa*pre/rho );
	double rhoinf = ppa->rhoinf;
	double uinf   = ppa->uinf;
	double vinf   = ppa->vinf;
	double preinf = ppa->Preinf;
	double uninf = uinf*nx + vinf*ny;
	double ainf  = sqrt( gammaa*preinf/rhoinf );

	double rhob, ub, vb, preb;		// boundary
    double Sb;
	if( un>=a )	// super out
	{
		rhob = rho;
		ub   = u;
		vb   = v;
		preb = pre;
	}
	else if( uninf<=-ainf )	// super in
	//else if( 0.5*(un+uninf) <= -0.5*(a+ainf) )	// super in
	{
		rhob = ppa->rhoinf;
		ub   = ppa->uinf;
		vb   = ppa->vinf;
		preb = ppa->Preinf;
	}
	else	//sub
	{
		double unb = 0.5*(un+uninf) + (a-ainf)*gam1r;
		double ab  = 0.5*(a + ainf) + 0.25*gamm1*(un-uninf);


		if( unb>=0.0 )	// sub out 
		{
			ub = u + nx*(unb-un);
			vb = v + ny*(unb-un);
			Sb = pre * pow( rho, -gammaa );
		}
		else	// sub in
		{
			ub = uinf + nx*(unb-uninf);
			vb = vinf + ny*(unb-uninf);
			Sb = preinf * pow( rhoinf, -gammaa );
		}

		rhob = pow( ab*ab/gammaa/Sb, gam1r );
		preb = rhob*ab*ab/gammaa;

		if(rhob<RMIN) { rhob = RMIN; printf("far rho<0.0\n" ); getchar(); }
		if(rhob>RMAX) { rhob = RMAX; printf("far rho>max\n" ); getchar(); }
		if(preb<PMIN) { preb = PMIN; printf("far pre<0.0\n" ); getchar(); }
		if(preb>PMAX) { preb = PMAX; printf("far pre>max\n" ); getchar(); }
	}

#if		VARIABLE==CONSERV
	ret[0] = rhob;
	ret[1] = rhob*ub;
	ret[2] = rhob*vb;
	double temp1 = preb*gam1r;
	double temp2 = 0.5*rhob*(ub*ub+vb*vb);
	ret[3] = temp1 + temp2;
	//ret[3] = preb*gam1r + 0.5*rhob*(ub*ub+vb*vb);
	/*if(it == 1 && flag == 0){
		printf( " %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e\n",rhob, ub, vb, preb, temp1, temp2, ret[0],ret[1],ret[2],ret[3]);
	}*/
#elif	VARIABLE==PRIMITI
	ret[0] = rhob;
	ret[1] = ub;
	ret[2] = vb;
	ret[3] = preb;
#endif

}


// 欧拉壁面边界条件
void bcSolidSurf( const double* pVar, const double *pVar0, const double *norm, double *ret )
{
#if		VARIABLE==CONSERV
	const double*pcv = pVar;
	double pre;
	if( CalPre( pVar, pVar0, pre )==1 )
		pcv = pVar0;

	double rho = pcv[0];
	double u   = pcv[1]/rho;
	double v   = pcv[2]/rho;
#elif	VARIABLE==PRIMITI
	double rho = pVar[0];
	double u   = pVar[1];
	double v   = pVar[2];
	double pre = pVar[3];
#endif

	if(rho<RMIN) { rho = RMIN; printf("far rho<0.0\n" ); getchar(); }
	if(rho>RMAX) { rho = RMAX; printf("far rho>max\n" ); getchar(); }
	if(pre<PMIN) { pre = PMIN; printf("far pre<0.0\n" ); getchar(); }
	if(pre>PMAX) { pre = PMAX; printf("far pre>max\n" ); getchar(); }

	double nx = norm[0];
	double ny = norm[1];

	double un = u*nx + v*ny;
	double ub = u - un*nx;
	double vb = v - un*ny;

	double rhob = rho;
	double preb = pre;
	
#if		VARIABLE==CONSERV
	ret[0] = rhob;
	ret[1] = rhob*ub;
	ret[2] = rhob*vb;
	ret[3] = preb*gam1r+0.5*rhob*(ub*ub+vb*vb);
#elif	VARIABLE==PRIMITI
	ret[0] = rhob;
	ret[1] = ub;
	ret[2] = vb;
	ret[3] = preb;
#endif
}


void bcvisSolidSurf( const double* cv, const double* cv0, para *ppa, double* presu )
{
	if( ppa->ivis!=1 )
	{
		PRINTERRORN( "Not vis solid surface!\n" );
	}

	if( ppa->bcTtype==ISOTWALL )
	{
		// Couette flow
		double pre;
		CalPre( cv, cv0, pre );
		//presu[0] = pre/(ppa->R*0.8);
		//presu[0] = 1.0;
		presu[3] = pre*gam1r;
		presu[0] = gamm1*cv[3]/(ppa->R*ppa->Twal);
		presu[1] = 0.0;
		presu[2] = 0.0;
		//presu[3] = cv[3];
	}
	else if( ppa->bcTtype == ADIAWALL )
	{
		presu[0] = cv[0];
		presu[1] = 0.0;
		presu[2] = 0.0;
		presu[3] = cv[3];
	}
	else
	{
		PRINTERRORN( "temperature bc condition error!\n" );
	}
}


// 对称边界条件
void bcSymmetry( 
	const double *uh, const double *uh0, const double *norm, double *ret )
{
	ret[0] = uh[0];
	ret[3] = uh[3];

	double un = uh[1]*norm[0] + uh[2]*norm[1];
	ret[1] = uh[1] - un*norm[0];
	ret[2] = uh[2] - un*norm[1];
}

// 对称边界条件
void bcSymmetry2( 
	const double *uh, const double *uh0, const double *norm, double *ret )
{
	ret[0] = uh[0];
	ret[3] = uh[3];

	double un2 = 2.0*(uh[1]*norm[0] + uh[2]*norm[1]);
	ret[1] = uh[1] - un2*norm[0];
	ret[2] = uh[2] - un2*norm[1];
}
