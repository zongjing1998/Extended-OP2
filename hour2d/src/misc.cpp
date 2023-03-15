// misc
#include <math.h>
#include <stdio.h>
//#include <conio.h>
#include <stdlib.h>

#include "misc.h"
#include "def2d.h"
#include "poly2d.h"

#include "m_math_def.h"

void GetValuede2d_aban( unsigned int nord, const double (*pcoe)[NEQ], 
				   unsigned int twh, double xi, double et, double* presult )
{
	double tmbas[MAXFVDOF];
	unsigned int nnof = DOFofOrd2d[nord];

	for( unsigned int i=0; i<nnof; i++ )
		tmbas[i] = BasisPolyTriDe( i, xi, et, twh );
		
	for( int j=0; j<NEQ; j++ )
	{
		presult[j] = pcoe[0][j] * tmbas[0];

		for( unsigned int i=1; i<nnof; i++ )
			presult[j] += pcoe[i][j] * tmbas[i];
	}
}


// return 1 意味着出了"问题" 使用cv0来计算，但这也可能出负，
// 一旦出负需要返回然后修改时间步长
int CalPre( const double * cv, const double * cv_0, double & pre )
{
	double rho = cv[0];

	if( (rho<RMIN) || (rho>RMAX) )
	{
		rho = cv_0[0];
		if(rho<RMIN) { rho = RMIN; printf("rho<min,reset\n" ); getchar(); }
		if(rho>RMAX) { rho = RMAX; printf("rho>max,reset\n" ); getchar(); }

		// calculate pressure
		pre = gamm1 * (cv_0[NEQ-1] - 0.5*(cv_0[1]*cv_0[1]+cv_0[2]*cv_0[2])/rho);

		if(pre<PMIN) { pre = PMIN; printf("pre<min,reset\n" ); getchar(); }
		if(pre>PMAX) { pre = PMAX; printf("pre>max,reset\n" ); getchar(); }

		return 1;
	}

	// calculate pressure
	pre = gamm1 * (cv[NEQ-1] - 0.5*(cv[1]*cv[1]+cv[2]*cv[2])/rho);
	
	if( (pre<PMIN)||(pre>PMAX) )
	{
		rho = cv_0[0];
		if(rho<RMIN) { rho = RMIN; printf("rho<min,reset\n" ); getchar(); }
		if(rho>RMAX) { rho = RMAX; printf("rho>max,reset\n" ); getchar(); }

		// calculate pressure
		pre = gamm1 * (cv_0[NEQ-1] - 0.5*(cv_0[1]*cv_0[1]+cv_0[2]*cv_0[2])/rho);

		if(pre<PMIN) { pre = PMIN; printf("pre<min,reset\n" ); getchar(); }
		if(pre>PMAX) { pre = PMAX; printf("pre>max,reset\n" ); getchar(); }

		return 1;
	}

	return 0;
}


double CalDt( unsigned int irk, double t, unsigned int nele, elem *pele, uvar *pvar, para *ppa )
{
	//printf( "caldt\n" );
	if( ppa->CFL<0.0 ) return ( -ppa->CFL );

	double dtmin = 9e+99;

	//#pragma omp parallel for
	for( int e=0; e<(int)nele; e++ )
	{
#if		EQUATIO==ADVECTI
		pvar[e].dt = ppa->CFL / (2.0*pvar[e].pdg+1.0) * pele[e].dh / 
			( m_abs(ppa->a) + m_abs(ppa->b) + TOL );
#elif	EQUATIO==BURGERS
		double U = m_abs( pvar[e].uh[irk][0][0] ) + TOL;

		pvar[e].dt = ppa->CFL / (2.0*pvar[e].pdg+1.0) * pele[e].dh / U / 
			( m_abs(ppa->a) + m_abs(ppa->b) + TOL );
#elif	EQUATIO==EULNSEQ
		double rho = pvar[e].uh[irk][0][0];

	#if		VARIABLE==CONSERV
		// 计算单元守恒量平均值的”压力“
		double pre;
		CalPre( pvar[e].uh[irk][0], pvar[e].uh[irk][0], pre );
		double u = pvar[e].uh[irk][0][1]/rho;
		double v = pvar[e].uh[irk][0][2]/rho;
		//double uh = sqrt( u*u + v*v );

	#elif	VARIABLE==PRIMITI
		double pre = pvar[e].uh[irk][0][NEQ-1];
		double U = sqrt( 
			pvar[e].uh[irk][0][1]*pvar[e].uh[irk][0][1] + 
			pvar[e].uh[irk][0][2]*pvar[e].uh[irk][0][2] );
	#endif

		double drey = 0.0;
		if( ppa->ivis==1 )
		{
			drey = 2.0 * gammaa / Pr;
			double T = pre / ( ppa->R * rho );
			double retot = ppa->Re*CalMu( T, ppa->c4suth );
			//pvar[e].retot = pvar[e].relam;// + pvar[e].returb;
			drey *= retot / rho;
		}

		double C = sqrt( gammaa*pre/rho );	// sound speed

		pvar[e].dt = ppa->CFL * pele[e].dxy[0]*pele[e].dxy[1]/
			( (m_abs(u)+C)*pele[e].dxy[1]+ (m_abs(v)+C)*pele[e].dxy[0]);

#endif
		dtmin = m_min( dtmin, pvar[e].dt );
	}

	return dtmin;
}

