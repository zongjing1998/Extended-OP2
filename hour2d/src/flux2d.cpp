// numerical flux

#include <memory.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
//#include <conio.h>
#include <math.h>

#include "misc.h"
#include "def2d.h"
#include "flux2d.h"

#include "m_math_def.h"


#define Entropy_Delta 0.1
// qL  : left primitive variables rho,u,v,p
// qR  : right primitive variables rho,u,v,p
// diss: dissipation term of numerical flux 
void FluxRoe( const double* qL, const double* qR, 
	          double hL, double hR, const double* norm, double* diss )
{
	double Ua, ra, ua, va, ha, aa, tm, tm1;
	double lm1, lm2, lm4;// lm3;
	double dr, du, dv, dp, dU;
	double nx, ny;
	double aa2, tp1, tp2, tp4;

	nx = norm[0];
	ny = norm[1];

	tm = sqrt( qR[0]/qL[0] );
	tm1= 1.0 + tm;
	ra = tm * qL[0];				//sqrt( qR[0] * qL[0] );
	ua = ( qL[1] + qR[1]*tm )/tm1;	//(1.0+tm);
	va = ( qL[2] + qR[2]*tm )/tm1;	//(1.0+tm);
	ha = ( hL    + hR   *tm )/tm1;	//(1.0+tm);
	aa2= gamm1*( ha - 0.5*(ua*ua+va*va) );
	aa = sqrt( aa2 );
	Ua = ua*nx + va*ny;

	lm1 = m_abs(Ua-aa);
	if( lm1 < Entropy_Delta ) // EntropyModi
		lm1 = (lm1*lm1+Entropy_Delta*Entropy_Delta)/(2.0*Entropy_Delta);

	lm2 = m_abs(Ua);
	if( lm2 < Entropy_Delta ) // EntropyModi
		lm2 = (lm2*lm2+Entropy_Delta*Entropy_Delta)/(2.0*Entropy_Delta);
	//lm3 = lm2;

	lm4 = m_abs(Ua+aa);
	if( lm4 < Entropy_Delta ) // EntropyModi
		lm4 = (lm4*lm4+Entropy_Delta*Entropy_Delta)/(2.0*Entropy_Delta);
		
	dr = qR[0] - qL[0];
	du = qR[1] - qL[1];
	dv = qR[2] - qL[2];
	dp = qR[3] - qL[3];
	dU = du*nx + dv*ny;

	tp1 = (dp - ra*aa*dU)*0.5/aa2;
	tp2 =  dr - dp/aa2;
	tp4 = (dp + ra*aa*dU)*0.5/aa2;

	tp1*= lm1;
	tp2*= lm2;
	tp4*= lm4;

	diss[0] = 
		tp1 + 
		tp2 + 
		tp4;
	diss[1] = 
		tp1*(ua-nx*aa) + 
		tp2*ua + lm2*ra*(du-nx*dU) +  
		tp4*(ua+nx*aa);
	diss[2] = 
		tp1*(va-ny*aa) + 
		tp2*va + lm2*ra*(dv-ny*dU) +  
		tp4*(va+ny*aa);
	diss[3] = 
		tp1*(ha-Ua*aa) + 
		tp2*0.5*(ua*ua+va*va) + lm2*ra*(ua*du+va*dv-Ua*dU) +  
		tp4*(ha+Ua*aa);
}
#undef delta

// qL, left primitive variables rho, u, v, h, p
// qR, right primitive variables rho, u, v, h, p
// cvL, left primitive variables rho, ru, rv, h, p
// cvR, right primitive variables rho, ru, rv, h, p
// flu numerical flux
void FluxRoeCu( double* cvL, double* cvR, double* norm, double* flux )
{
	double ra, ua, va, ha, aa, tm, tm1;//Ua, 
	double lm1, lm2, lm3, lm4;
	//double dr, du, dv, dp, dU;
	double a1, a2, a3, a4;
	//double fl[NEQ], fr[NEQ], hL, hR;
	double du1, du2, du3, du4;

	//double nx = norm[0];
	//double ny = norm[1];

	tm = sqrt( cvR[0] / cvL[0] );
	tm1= 1.0 + tm;
	ra = tm * cvL[0];//sqrt( cvR[0] * cvL[0] );
	ua = ( cvL[1] / cvL[0] + cvR[1] / cvR[0] * tm ) / tm1;//( 1.0 + tm );
	va = ( cvL[2] / cvL[0] + cvR[2] / cvR[0] * tm ) / tm1;//( 1.0 + tm );
	ha = ( cvL[NEQ-1] + cvR[NEQ-1] * tm ) / tm1;//( 1.0 + tm );
	aa = sqrt( gamm1 * ( ha - ( ua * ua + va * va ) / 2.0 ) );

	lm1 = m_abs(ua-aa);
	lm2 = m_abs(ua);
	lm3 = lm2;
	lm4 = m_abs(ua+aa);
		
	du1 = cvR[0] - cvL[0];
	du2 = cvR[1] - cvL[1];
	du3 = cvR[2] - cvL[2];
	du4 = cvR[0]*cvR[NEQ-1]-cvR[NEQ] - cvL[0]*cvL[NEQ-1]+cvL[NEQ];

	a3 = gamm1 / aa / aa * ( 
		 (ha-ua*ua-va*va)*du1 + ua*du2 + va*du3 - du4 );
	a2 = du3 - va*du1;
	a4 = 0.5*(du1 - a3 +(du2-ua*du1)/aa);
	a1 = du1 - a3 - a4;

	flux[0] = 0.5*( cvL[1] + cvR[1] - 
		a1*lm1 - 
		a3*lm3 - 
		a4*lm4 );
	flux[1] = 0.5*( cvL[1]*cvL[1]/cvL[0]+cvL[NEQ] + cvR[1]*cvR[1]/cvR[0]+cvR[NEQ] - 
		a1*lm1*(ua-aa) - 
		a3*lm3*ua - 
		a4*lm4*(ua+aa) );
	flux[2] = 0.5*( cvL[1]*cvL[2]/cvL[0] + cvR[1]*cvR[2]/cvR[0] - 
		a1*lm1*va - 
		a2*lm2 - 
		a3*lm3*va - 
		a4*lm4*va );
	flux[NEQ-1] = 0.5*( cvL[1]*cvL[NEQ-1]/cvL[0] + cvR[1]*cvR[NEQ-1]/cvR[0] - 
		a1*lm1*(ha-ua*aa) - 
		a2*lm2*va - 
		a3*lm3*0.5*(ua*ua+va*va) - 
		a4*lm4*(ha+ua*aa) );
}

// cvL, left conservation variables rho, ru, rv, e
// cvR, right conservation variables rho, ru, rv, e
// flux numerical flux
// norm nx, ny, ( nz )
// local Lax-Friedrichs flux
// The Lax¨CFriedrichs (LF) flux and the local LF (LLF) flux
void FluxLLF( const double* cvL, const double* cvL0, 
			  const double* cvR, const double* cvR0, const double* norm, double* flux )
{
#if		EQUATIO==ADVECTI
	double outnorm = norm[0] + norm[1];
	for( unsigned int j=0; j<NEQ; j++ )
	{
	if( outnorm>TOL )
		flux[j] = outnorm*cvL[j];
	else if( outnorm<-TOL )
		flux[j] = outnorm*cvR[j];
	else
		flux[j] = 0.0;
	}
#elif	EQUATIO==BURGERS
	double outnorm = norm[0] + norm[1];
	for( unsigned int j = 0; j < NEQ; j++ )
	{
	double mlam;
	mlam = m_max( m_abs(outnorm*cvL[j]), m_abs(outnorm*cvR[j]) );
	flux[j] = 0.5*( 0.5*outnorm*cvL[j]*cvL[j]+0.5*outnorm*cvR[j]*cvR[j] - mlam*(cvR[j]-cvL[j]) );	
	}
#elif	EQUATIO==EULNSEQ
	double preL, preR;
	double vnL, vnR;
	double aL, aR;
	double mlambda;

	const double* pcvL = cvL;
	if( CalPre( cvL, cvL0, preL ) )
		pcvL = cvL0;

	const double* pcvR = cvR;
	if( CalPre( cvR, cvR0, preR ) )
		pcvR = cvR0;

	aL = sqrt( gammaa * preL / pcvL[0] );
	aR = sqrt( gammaa * preR / pcvR[0] );
	vnL = ( norm[0] * pcvL[1] + norm[1] * pcvL[2] ) / pcvL[0];
	vnR = ( norm[0] * pcvR[1] + norm[1] * pcvR[2] ) / pcvR[0];
	mlambda = m_max( m_abs( vnL ) + aL, m_abs( vnR ) + aR );
	//double a = sqrt( gamma*m_max(preL,preR)/m_min(pcvL[0],pcvR[0]) );
	//mlambda = m_max( m_abs( vnL ), m_abs( vnR ) ) + a;

	double prelr = preL+preR;
	flux[0] = 0.5*( pcvL[0]*vnL               + pcvR[0]*vnR       - mlambda*(pcvR[0]-pcvL[0]) );
	flux[1] = 0.5*( pcvL[1]*vnL+prelr*norm[0] + pcvR[1]*vnR       - mlambda*(pcvR[1]-pcvL[1]) );
	flux[2] = 0.5*( pcvL[2]*vnL+prelr*norm[1] + pcvR[2]*vnR       - mlambda*(pcvR[2]-pcvL[2]) );
	flux[3] = 0.5*((pcvL[3]+preL)*vnL         +(pcvR[3]+preR)*vnR - mlambda*(pcvR[3]-pcvL[3]) );
#endif
}

// cvL, left conservation variables rho, ru, rv, e
// cvR, right conservation variables rho, ru, rv, e
// fh numerical fh
// norm nx, ny, ( nz )
void FluxSel( 
	unsigned int sel, 
	const double *cvL, const double* cvL0, 
	const double *cvR, const double* cvR0, const double *norm, double *fh )
{
	switch( sel )
	{
	case FLUXNND:
		printf( "fluxnnd under developing!\n" );
		getchar();
        break;

	case FLUXROE:

#if		EQUATIO==EULNSEQ
		double qL[NEQ], qR[NEQ];// r,u,v,p
		double hL, hR;
		double diss[NEQ];

		const double* pcvL;
		pcvL = cvL;
		if( CV2PV( cvL, cvL0, qL )==1 )
			pcvL = cvL0;

		const double* pcvR;
		pcvR = cvR;
		if( CV2PV( cvR, cvR0, qR )==1 )
			pcvR = cvR0;

		// h = (e+p)/rho
		hL = ( pcvL[NEQ-1] + qL[NEQ-1] ) / pcvL[0];
		hR = ( pcvR[NEQ-1] + qR[NEQ-1] ) / pcvR[0];
		FluxRoe( qL, qR, hL, hR, norm, diss );
		
		double UnL;
		UnL = qL[1]*norm[0] + qL[2]*norm[1];
		double UnR;
		UnR = qR[1]*norm[0] + qR[2]*norm[1];

		double preLR;
		preLR = qL[NEQ-1] + qR[NEQ-1];
		fh[0]      = 0.5*( pcvL[0]*UnL + pcvR[0]*UnR                 - diss[0] );
		fh[1]      = 0.5*( pcvL[1]*UnL + pcvR[1]*UnR + preLR*norm[0] - diss[1] );
		fh[2]      = 0.5*( pcvL[2]*UnL + pcvR[2]*UnR + preLR*norm[1] - diss[2] );
		fh[NEQ-1] = 0.5*( (pcvL[NEQ-1]+qL[NEQ-1])*UnL + 
			(pcvR[NEQ-1]+qR[NEQ-1])*UnR - diss[NEQ-1] );
#else
		printf( "fluxroe is for euler/ns only!\n" );
		getchar();
#endif
        break;

	case FLUXLLF:
		FluxLLF( cvL, cvL0, cvR, cvR0, norm, fh );
        break;

	case FLUXKIN:
		printf( "fluxkin under developing!\n" );
		getchar();
        break;

	case FLUXOSH:
		printf( "fluxosh under developing!\n" );
		getchar();
        break;

	default:
		printf( "undefined flux type!\n" );
		getchar();
        // for test!
		double cvavg[NEQ];
		for( unsigned int j=0; j<NEQ; j++ )
			cvavg[j] = 0.5*( cvL[j] + cvR[j] );

		double pre;
		CalPre( cvavg, cvavg, pre );

		fh[0] = cvavg[1] * norm[0] + cvavg[2] * norm[1];
		double Un;
		Un = fh[0] / cvavg[0];
		fh[1] = cvavg[1] * Un + pre * norm[0];
		fh[2] = cvavg[2] * Un + pre * norm[1];
		fh[3] =(cvavg[3]+pre) * Un;
	}
}
