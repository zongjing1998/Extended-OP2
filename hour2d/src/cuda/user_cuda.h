#include "../solio.h"
#include "../recons.h"
#include "../poly2d.h"
#include "../gridparam.h"
#include "../ic.h"
#include "../quadrature.h"
#include "../gridio.h"
#include "../rhs.h"
#include "../fluiddef.h"
#include "../bc2d.h"
#include "../misc.h"
#include "../flux2d.h"
#include "../rcm.h"
#include "../solparam.h"
#include "../m_math_def.h"
#include "../limiter.h"
#include "../postproc.h"
#include "../def2d.h"
#include "../poly1d.h"
#ifndef Entropy_Delta
#define Entropy_Delta 0.1
#endif
__device__ void CV2FluxIN2d_gpu(const double *cv, const double *cv0, double *norm, double *FluxIN);
__device__ void PV2FluxI_gpu(const double *pv, double (*fluxI)[NEQ]);
__device__ void FluxRoe_gpu( const double* qL, const double* qR,
	          double hL, double hR, const double* norm, double* diss );
__device__ void bcSymmetry_gpu(
	const double *uh, const double *uh0, const double *norm, double *ret );
__device__ double CalMu_gpu(double Tdimles, double c4suth);
__device__ void bcvisSolidSurf_gpu( const double* cv, const double* cv0, para *ppa, double* presu );
__device__ void CalFluxV_gpu(
    const double *uh, const double *uh0, double (*graduh)[NEQ], para *ppa,
    double (*fluxV)[NEQ]);
__device__ void CV2GradPV_gpu(
    const double *cv, double (*gradcv)[NEQ], double Tcoe, double (*gradpv)[NEQ]);
__device__ void PV2FluxV_gpu(
    double *pv, double (*gradpv)[NEQ], double miu, double kap, double (*fluxV)[NEQ]);
__device__ int CV2FluxI_gpu(const double *cv, const double *cv0, double (*fluxI)[NEQ]);
__device__ int CV2PV_gpu(const double *cv, const double *cv0, double *pv);
__device__ void FluxLLF_gpu( const double* cvL, const double* cvL0,
			  const double* cvR, const double* cvR0, const double* norm, double* flux );
__device__ void bcSymmetry2_gpu(
	const double *uh, const double *uh0, const double *norm, double *ret );
__device__ void bcSolidSurf_gpu( const double* pVar, const double *pVar0, const double *norm, double *ret );
__device__ void CalTau_gpu(double (*gradpv)[NEQ], double miu, double (*tau)[NDIM]);
__device__ void bcFarField_gpu(
	const double* uh, const double* uh0, const double* norm, para *ppa, double* ret);
__device__ void GetVar_gpu(unsigned int ndof, const double (*coe)[NEQ],
		   const double *basis, double *ret);
__device__ int CalPre_gpu( const double * cv, const double * cv_0, double & pre );
__device__ void FluxSel_gpu(
	unsigned int sel,
	const double *cvL, const double* cvL0,
	const double *cvR, const double* cvR0, const double *norm, double *fh );
__device__ void FluxSel_gpu(
	unsigned int sel,
	const double *cvL, const double* cvL0,
	const double *cvR, const double* cvR0, const double *norm, double *fh )
{
	switch( sel )
	{
	case FLUXNND:
		printf( "fluxnnd under developing!\n" );
		/*getchar();*/
        break;

	case FLUXROE:

#if		EQUATIO==EULNSEQ
		double qL[NEQ], qR[NEQ];
		double hL, hR;
		double diss[NEQ];

		const double* pcvL;
		pcvL = cvL;
		if( CV2PV_gpu( cvL, cvL0, qL )==1 )
			pcvL = cvL0;

		const double* pcvR;
		pcvR = cvR;
		if( CV2PV_gpu( cvR, cvR0, qR )==1 )
			pcvR = cvR0;

		hL = ( pcvL[NEQ-1] + qL[NEQ-1] ) / pcvL[0];
		hR = ( pcvR[NEQ-1] + qR[NEQ-1] ) / pcvR[0];
		FluxRoe_gpu( qL, qR, hL, hR, norm, diss );

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
		/*getchar();*/
#endif
        break;

	case FLUXLLF:
		FluxLLF_gpu( cvL, cvL0, cvR, cvR0, norm, fh );
        break;

	case FLUXKIN:
		printf( "fluxkin under developing!\n" );
		/*getchar();*/
        break;

	case FLUXOSH:
		printf( "fluxosh under developing!\n" );
		/*getchar();*/
        break;

	default:
		printf( "undefined flux type!\n" );
		/*getchar();*/

		double cvavg[NEQ];
		for( unsigned int j=0; j<NEQ; j++ )
			cvavg[j] = 0.5*( cvL[j] + cvR[j] );

		double pre;
		CalPre_gpu( cvavg, cvavg, pre );

		fh[0] = cvavg[1] * norm[0] + cvavg[2] * norm[1];
		double Un;
		Un = fh[0] / cvavg[0];
		fh[1] = cvavg[1] * Un + pre * norm[0];
		fh[2] = cvavg[2] * Un + pre * norm[1];
		fh[3] =(cvavg[3]+pre) * Un;
	}

}
__device__ int CalPre_gpu( const double * cv, const double * cv_0, double & pre )
{
	double rho = cv[0];

	if( (rho<RMIN) || (rho>RMAX) )
	{
		rho = cv_0[0];
		if(rho<RMIN) { rho = RMIN; printf("rho<min,reset\n" ); /*getchar();*/ }
		if(rho>RMAX) { rho = RMAX; printf("rho>max,reset\n" ); /*getchar();*/ }

		pre = gamm1 * (cv_0[NEQ-1] - 0.5*(cv_0[1]*cv_0[1]+cv_0[2]*cv_0[2])/rho);

		if(pre<PMIN) { pre = PMIN; printf("pre<min,reset\n" ); /*getchar();*/ }
		if(pre>PMAX) { pre = PMAX; printf("pre>max,reset\n" ); /*getchar();*/ }

		return 1;
	}

	pre = gamm1 * (cv[NEQ-1] - 0.5*(cv[1]*cv[1]+cv[2]*cv[2])/rho);

	if( (pre<PMIN)||(pre>PMAX) )
	{
		rho = cv_0[0];
		if(rho<RMIN) { rho = RMIN; printf("rho<min,reset\n" ); /*getchar();*/ }
		if(rho>RMAX) { rho = RMAX; printf("rho>max,reset\n" ); /*getchar();*/ }

		pre = gamm1 * (cv_0[NEQ-1] - 0.5*(cv_0[1]*cv_0[1]+cv_0[2]*cv_0[2])/rho);

		if(pre<PMIN) { pre = PMIN; printf("pre<min,reset\n" ); /*getchar();*/ }
		if(pre>PMAX) { pre = PMAX; printf("pre>max,reset\n" ); /*getchar();*/ }

		return 1;
	}

	return 0;

}
__device__ void GetVar_gpu(unsigned int ndof, const double (*coe)[NEQ],
		   const double *basis, double *ret)
{







	double tret[NEQ] = {0.0};
	for (int i = 0; i < ndof; i++)
		for (int j = 0; j < NEQ; j++)
		{
			tret[j] += coe[i][j] * basis[i];

		}
	for (int j = 0; j < NEQ; j++)
	{
		ret[j] = tret[j];

	}

}
__device__ void bcFarField_gpu(
	const double* uh, const double* uh0, const double* norm, para *ppa, double* ret)
{
#if		VARIABLE==CONSERV
	const double *cv = uh;
	double pre;
	if( CalPre_gpu( uh, uh0, pre ) )
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

	if(rho<RMIN) { rho = RMIN; printf("far rho<0.0\n" ); /*getchar();*/ }
	if(rho>RMAX) { rho = RMAX; printf("far rho>max\n" ); /*getchar();*/ }
	if(pre<PMIN) { pre = PMIN; printf("far pre<0.0\n" ); /*getchar();*/ }
	if(pre>PMAX) { pre = PMAX; printf("far pre>max\n" ); /*getchar();*/ }

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

	double rhob, ub, vb, preb;
    double Sb;
	if( un>=a )
	{
		rhob = rho;
		ub   = u;
		vb   = v;
		preb = pre;
	}
	else if( uninf<=-ainf )

	{
		rhob = ppa->rhoinf;
		ub   = ppa->uinf;
		vb   = ppa->vinf;
		preb = ppa->Preinf;
	}
	else
	{
		double unb = 0.5*(un+uninf) + (a-ainf)*gam1r;
		double ab  = 0.5*(a + ainf) + 0.25*gamm1*(un-uninf);

		if( unb>=0.0 )
		{
			ub = u + nx*(unb-un);
			vb = v + ny*(unb-un);
			Sb = pre * pow( rho, -gammaa );
		}
		else
		{
			ub = uinf + nx*(unb-uninf);
			vb = vinf + ny*(unb-uninf);
			Sb = preinf * pow( rhoinf, -gammaa );
		}

		rhob = pow( ab*ab/gammaa/Sb, gam1r );
		preb = rhob*ab*ab/gammaa;

		if(rhob<RMIN) { rhob = RMIN; printf("far rho<0.0\n" ); /*getchar();*/ }
		if(rhob>RMAX) { rhob = RMAX; printf("far rho>max\n" ); /*getchar();*/ }
		if(preb<PMIN) { preb = PMIN; printf("far pre<0.0\n" ); /*getchar();*/ }
		if(preb>PMAX) { preb = PMAX; printf("far pre>max\n" ); /*getchar();*/ }
	}

#if		VARIABLE==CONSERV
	ret[0] = rhob;
	ret[1] = rhob*ub;
	ret[2] = rhob*vb;
	double temp1 = preb*gam1r;
	double temp2 = 0.5*rhob*(ub*ub+vb*vb);
	ret[3] = temp1 + temp2;


#elif	VARIABLE==PRIMITI
	ret[0] = rhob;
	ret[1] = ub;
	ret[2] = vb;
	ret[3] = preb;
#endif


}
__device__ void CalTau_gpu(double (*gradpv)[NEQ], double miu, double (*tau)[NDIM])
{
	double miu2d3 = miu * 0.66666666666666666666666666666667;
	tau[0][0] = (2 * gradpv[0][1] - gradpv[1][2]) * miu2d3;
	tau[0][1] = (gradpv[1][1] + gradpv[0][2]) * miu;

	tau[1][0] = tau[0][1];
	tau[1][1] = (2 * gradpv[1][2] - gradpv[0][1]) * miu2d3;

}
__device__ void bcSolidSurf_gpu( const double* pVar, const double *pVar0, const double *norm, double *ret )
{
#if		VARIABLE==CONSERV
	const double*pcv = pVar;
	double pre;
	if( CalPre_gpu( pVar, pVar0, pre )==1 )
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

	if(rho<RMIN) { rho = RMIN; printf("far rho<0.0\n" ); /*getchar();*/ }
	if(rho>RMAX) { rho = RMAX; printf("far rho>max\n" ); /*getchar();*/ }
	if(pre<PMIN) { pre = PMIN; printf("far pre<0.0\n" ); /*getchar();*/ }
	if(pre>PMAX) { pre = PMAX; printf("far pre>max\n" ); /*getchar();*/ }

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
__device__ void bcSymmetry2_gpu(
	const double *uh, const double *uh0, const double *norm, double *ret )
{
	ret[0] = uh[0];
	ret[3] = uh[3];

	double un2 = 2.0*(uh[1]*norm[0] + uh[2]*norm[1]);
	ret[1] = uh[1] - un2*norm[0];
	ret[2] = uh[2] - un2*norm[1];

}
__device__ void FluxLLF_gpu( const double* cvL, const double* cvL0,
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
	if( CalPre_gpu( cvL, cvL0, preL ) )
		pcvL = cvL0;

	const double* pcvR = cvR;
	if( CalPre_gpu( cvR, cvR0, preR ) )
		pcvR = cvR0;

	aL = sqrt( gammaa * preL / pcvL[0] );
	aR = sqrt( gammaa * preR / pcvR[0] );
	vnL = ( norm[0] * pcvL[1] + norm[1] * pcvL[2] ) / pcvL[0];
	vnR = ( norm[0] * pcvR[1] + norm[1] * pcvR[2] ) / pcvR[0];
	mlambda = m_max( m_abs( vnL ) + aL, m_abs( vnR ) + aR );



	double prelr = preL+preR;
	flux[0] = 0.5*( pcvL[0]*vnL               + pcvR[0]*vnR       - mlambda*(pcvR[0]-pcvL[0]) );
	flux[1] = 0.5*( pcvL[1]*vnL+prelr*norm[0] + pcvR[1]*vnR       - mlambda*(pcvR[1]-pcvL[1]) );
	flux[2] = 0.5*( pcvL[2]*vnL+prelr*norm[1] + pcvR[2]*vnR       - mlambda*(pcvR[2]-pcvL[2]) );
	flux[3] = 0.5*((pcvL[3]+preL)*vnL         +(pcvR[3]+preR)*vnR - mlambda*(pcvR[3]-pcvL[3]) );
#endif

}
__device__ int CV2PV_gpu(const double *cv, const double *cv0, double *pv)
{
	double pre;
	int ret = CalPre_gpu(cv, cv0, pre);

	const double *pcv = cv;
	if (ret == 1)
		pcv = cv0;

	pv[0] = pcv[0];
	pv[1] = pcv[1] / pcv[0];
	pv[2] = pcv[2] / pcv[0];
	pv[NEQ - 1] = pre;

	return ret;

}
__device__ int CV2FluxI_gpu(const double *cv, const double *cv0, double (*fluxI)[NEQ])
{
	double pre;
	int ret = CalPre_gpu(cv, cv0, pre);

	const double *pcv = cv;
	if (ret == 1)
		pcv = cv0;

	double u = pcv[1] / pcv[0];
	double v = pcv[2] / pcv[0];
	double e = pcv[NEQ - 1];

	double ep = e + pre;

	fluxI[0][0] = pcv[1];
	fluxI[0][1] = fluxI[0][0] * u + pre;
	fluxI[0][2] = fluxI[0][0] * v;
	fluxI[0][NEQ - 1] = ep * u;

	fluxI[1][0] = pcv[2];
	fluxI[1][1] = fluxI[0][2];
	fluxI[1][2] = fluxI[1][0] * v + pre;
	fluxI[1][NEQ - 1] = ep * v;

	return ret;

}
__device__ void PV2FluxV_gpu(
    double *pv, double (*gradpv)[NEQ], double miu, double kap, double (*fluxV)[NEQ])
{
	double u = pv[1];
	double v = pv[2];
	double tau[NDIM][NDIM];
	CalTau_gpu(gradpv, miu, tau);

	fluxV[0][0] = 0.0;
	fluxV[0][1] = tau[0][0];
	fluxV[0][2] = tau[0][1];
	fluxV[0][NEQ - 1] = u * fluxV[0][1] + v * fluxV[0][2] + kap * gradpv[0][NEQ - 1];

	fluxV[1][0] = 0.0;
	fluxV[1][1] = tau[1][0];
	fluxV[1][2] = tau[1][1];
	fluxV[1][NEQ - 1] = u * fluxV[1][1] + v * fluxV[1][2] + kap * gradpv[1][NEQ - 1];

}
__device__ void CV2GradPV_gpu(
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
	gradpv[0][NEQ - 1] *= Tcoe;

	gradpv[1][0] = gradcv[1][0];
	gradpv[1][1] = (gradcv[1][1] - gradcv[1][0] * u) * rr1;
	gradpv[1][2] = (gradcv[1][2] - gradcv[1][0] * v) * rr1;
	gradpv[1][NEQ - 1] = (gradcv[1][NEQ - 1] - gradcv[1][0] * E) * rr1 - u * gradpv[1][1] - v * gradpv[1][2];
	gradpv[1][NEQ - 1] *= Tcoe;

}
__device__ void CalFluxV_gpu(
    const double *uh, const double *uh0, double (*graduh)[NEQ], para *ppa,
    double (*fluxV)[NEQ])
{
#if EQUATIO == EULNSEQ

#if VARIABLE == CONSERV
	double pv[NEQ];
	int ret = CV2PV_gpu(uh, uh0, pv);
	const double *cv = uh;
	if (ret == 1)
		cv = uh0;
#if AUXVARIABLE == CONSERV
	double gradpv[NDIM][NEQ];
	CV2GradPV_gpu(cv, graduh, ppa->Tcoe, gradpv);
#elif AUXVARIABLE == PRIMITI
	double(&gradpv)[NDIM][NEQ] = graduh;
#endif

#elif VARIABLE == PRIMITI
	const double *pv = uh;
	double(&gradpv)[NDIM][NEQ] = graduh;
#endif

	double T = pv[NEQ - 1] / (pv[0] * ppa->R);
	double miu = ppa->Re * CalMu_gpu(T, ppa->c4suth);

	PV2FluxV_gpu(pv, gradpv, miu, miu * ppa->KapdMiu, fluxV);
#endif

}
__device__ void bcvisSolidSurf_gpu( const double* cv, const double* cv0, para *ppa, double* presu )
{
	if( ppa->ivis!=1 )
	{
		PRINTERRORN( "Not vis solid surface!\n" );
	}

	if( ppa->bcTtype==ISOTWALL )
	{

		double pre;
		CalPre_gpu( cv, cv0, pre );


		presu[3] = pre*gam1r;
		presu[0] = gamm1*cv[3]/(ppa->R*ppa->Twal);
		presu[1] = 0.0;
		presu[2] = 0.0;

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
__device__ double CalMu_gpu(double Tdimles, double c4suth)
{

	return 1.0;

	return Tdimles;

	double dnum = 1.0 + c4suth;
	double dden = Tdimles + c4suth;
	return Tdimles * sqrt(Tdimles) * dnum / dden;

}
__device__ void bcSymmetry_gpu(
	const double *uh, const double *uh0, const double *norm, double *ret )
{
	ret[0] = uh[0];
	ret[3] = uh[3];

	double un = uh[1]*norm[0] + uh[2]*norm[1];
	ret[1] = uh[1] - un*norm[0];
	ret[2] = uh[2] - un*norm[1];

}
__device__ void FluxRoe_gpu( const double* qL, const double* qR,
	          double hL, double hR, const double* norm, double* diss )
{
	double Ua, ra, ua, va, ha, aa, tm, tm1;
	double lm1, lm2, lm4;
	double dr, du, dv, dp, dU;
	double nx, ny;
	double aa2, tp1, tp2, tp4;

	nx = norm[0];
	ny = norm[1];

	tm = sqrt( qR[0]/qL[0] );
	tm1= 1.0 + tm;
	ra = tm * qL[0];
	ua = ( qL[1] + qR[1]*tm )/tm1;
	va = ( qL[2] + qR[2]*tm )/tm1;
	ha = ( hL    + hR   *tm )/tm1;
	aa2= gamm1*( ha - 0.5*(ua*ua+va*va) );
	aa = sqrt( aa2 );
	Ua = ua*nx + va*ny;

	lm1 = m_abs(Ua-aa);
	if( lm1 < Entropy_Delta )
		lm1 = (lm1*lm1+Entropy_Delta*Entropy_Delta)/(2.0*Entropy_Delta);

	lm2 = m_abs(Ua);
	if( lm2 < Entropy_Delta )
		lm2 = (lm2*lm2+Entropy_Delta*Entropy_Delta)/(2.0*Entropy_Delta);


	lm4 = m_abs(Ua+aa);
	if( lm4 < Entropy_Delta )
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
__device__ void PV2FluxI_gpu(const double *pv, double (*fluxI)[NEQ])
{
	double rho = pv[0];
	double u = pv[1];
	double v = pv[2];
	double pre = pv[NEQ - 1];
	double e = pre * gam1r + 0.5 * rho * (u * u + v * v);
	double ep = e + pre;

	fluxI[0][0] = rho * u;
	fluxI[0][1] = fluxI[0][0] * u + pre;
	fluxI[0][2] = fluxI[0][0] * v;
	fluxI[0][NEQ - 1] = ep * u;

	fluxI[1][0] = rho * v;
	fluxI[1][1] = fluxI[0][2];
	fluxI[1][2] = fluxI[1][0] * v + pre;
	fluxI[1][NEQ - 1] = ep * v;

}
__device__ void CV2FluxIN2d_gpu(const double *cv, const double *cv0, double *norm, double *FluxIN)
{
#if EQUATIO == EULNSEQ
	double e;
#if VARIABLE == CONSERV
	double pv[NEQ];
	int ret = CV2PV_gpu(cv, cv0, pv);
	e = cv[NEQ - 1];
	if (ret == 1)
		e = cv0[NEQ - 1];
#elif VARIABLE == PRIMITI
	const double *pv = cv;
	e = pv[NEQ - 1] * gam1r + 0.5 * pv[0] * (pv[1] * pv[1] + pv[2] * pv[2]);
#endif
	double un = norm[0] * pv[1] + norm[1] * pv[2];
	FluxIN[0] = pv[0] * un;
	FluxIN[1] = pv[1] * FluxIN[0] + norm[0] * pv[NEQ - 1];
	FluxIN[2] = pv[2] * FluxIN[0] + norm[1] * pv[NEQ - 1];
	FluxIN[NEQ - 1] = (pv[NEQ - 1] + e) * un;
#endif

}
