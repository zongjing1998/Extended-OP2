// param.cpp
// 2011 07 28 liming
// read parameters
// process parameters
// write screen

#include <math.h>
//#include <conio.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "bc2d.h"
#include "flux2d.h"
#include "solparam.h"

#include "m_math_def.h"

void ReadParam( para *ppa, const char* strParaFile )
{
	FILE *sparam = NULL;
       	//errno_t err;
	
	// read parameter file (will fail if file does not exist)
	/*if( (err = fopen_s( &sparam, strParaFile, "r" )) !=0 )
	{
		printf( "The file %s was not opened\n", strParaFile );
		exit(1);
	}*/
        sparam = fopen(strParaFile, "r");
        if(sparam == NULL){
                printf( "The file %s was not opened\n", strParaFile );
                exit(1);
        }
	char stemp[256];
	fgets( stemp, 256, sparam );
	ppa->ivis = atoi( stemp );
	fgets( stemp, 256, sparam );
	ppa->fluxtype = atoi( stemp );
	fgets( stemp, 256, sparam );
	ppa->pDG = atoi( stemp );
	ppa->pFV = ppa->pDG;

	fgets( stemp, 256, sparam );
	ppa->isteady = atoi( stemp );
	fgets( stemp, 256, sparam );
	ppa->itmax = atoi( stemp );
	fgets( stemp, 256, sparam );
	ppa->wffre = atoi( stemp );
	fgets( stemp, 256, sparam );
	ppa->wsfre = atoi( stemp );
	fgets( stemp, 256, sparam );
	ppa->dtfre = atoi( stemp );
	fgets( stemp, 256, sparam );
	ppa->iconti = atoi( stemp );
	fgets( stemp, 256, sparam );
	ppa->resi = atoi( stemp );
	fgets( stemp, 256, sparam );
	ppa->CFLIncStp = atoi( stemp );
    fgets( stemp, 256, sparam );
    ppa->iRuKu = atoi( stemp );
	fgets( stemp, 256, sparam );
	ppa->bcTtype = atoi( stemp );
    fgets( stemp, 256, sparam );
    ppa->ImResSmoStp = atoi( stemp );

    fgets( stemp, 256, sparam );
    ppa->ImResSmoCoe = atof( stemp );

    fgets( stemp, 256, sparam );
    ppa->ChgMaxStp = atoi( stemp );

	fgets( stemp, 256, sparam );
	ppa->M_TVB = atof( stemp );

	fgets( stemp, 256, sparam );
	ppa->ChgMaxReduc = atof( stemp );

	fgets( stemp, 256, sparam );
	ppa->ChgMaxRatio = atof( stemp );

	fgets( stemp, 256, sparam );
	ppa->Re = atof( stemp );
	fgets( stemp, 256, sparam );
	ppa->Ma = atof( stemp );
	fgets( stemp, 256, sparam );
	ppa->RatioP = atof( stemp );
	fgets( stemp, 256, sparam );
	ppa->Tinfd = atof( stemp );
	fgets( stemp, 256, sparam );
	ppa->Twald = atof( stemp );
	fgets( stemp, 256, sparam );
	ppa->AoA = atof( stemp );
	fgets( stemp, 256, sparam );
	ppa->CFL = atof( stemp );
	fgets( stemp, 256, sparam );
	ppa->tmax = atof( stemp );

	fgets( stemp, 256, sparam );	
	int pos = strcspn( stemp, " \t\n" );
	stemp[pos] = 0;
	strcpy( ppa->gfile, stemp );

	fclose( sparam );

	// configuration
	printf( "*********************************************\n" );
	printf( "        NSNNDGFV : RELEASE 0.99, 2012        \n" );
	printf( "*********************************************\n" );

#if		EQUATIO==ADVECTI
	printf( "advective equation\n" );

#elif	EQUATIO==BURGERS
	//ppa->ndcoordt = NULL;
	printf( "Burgers equation\n" );

#elif	EQUATIO==EULNSEQ
	if( ppa->ivis==0 ) printf( "Euler eqs\n" );
	else printf( "Navier-Stokes eqs\n" );
	printf( "*********************************************\n" );

	// temperature B.C.
	if( ppa->ivis==1 )
	{
		if( ppa->bcTtype == ISOTWALL ) printf( "isothermal wall\n" );
		if( ppa->bcTtype == ADIAWALL ) printf( "adiabatic  wall\n" );

		printf( "reynolds number %e\n", ppa->Re );
		ppa->Re = 1.0/ppa->Re;
	}

	//printf( "Froude number %f\n", ppa->Fr );
	// uniform free stream conditions // external flow around a body
	printf( "free stream mach number %f\n", ppa->Ma );
	printf( "angle of attack %f\n", ppa->AoA );

	double Ma2 = ppa->Ma*ppa->Ma;
	ppa->rhoinf = 1.0;
	ppa->Preinf = 1.0 / (gammaa*Ma2);
	ppa->Tinf   = 1.0;
	ppa->einf   = 0.5 + ppa->Preinf*gam1r;

	ppa->AoA	= ppa->AoA * PI/180.0;
	ppa->cosAoA	= cos( ppa->AoA );
	ppa->sinAoA	= sin( ppa->AoA );
	ppa->uinf	= ppa->cosAoA;
	ppa->vinf	= ppa->sinAoA;

	ppa->Ttinf	   = 1.0 + 0.5*gamm1 * Ma2;
	ppa->Pretinf   = pow( ppa->Ttinf, gammaa*gam1r ) * ppa->Preinf;
	ppa->Sinfn1	   = pow( ppa->rhoinf, gammaa ) / ppa->Preinf;
	ppa->Rhoq2d2n1 = 2.0; ///(ppa->rhoinf);//*(u*u+v*v);

	ppa->Twal	= ppa->Twald / ppa->Tinfd;
	ppa->c4suth = 110.4 / ppa->Tinfd;
	ppa->R		= 1.0 / (gammaa*Ma2);
	ppa->Tcoe	= gamm1/ppa->R;
	ppa->KapdMiu= 1.0 / (Pr*gamm1*Ma2);
#endif
	printf( "*********************************************\n" );

	if( ppa->isteady==UNSTEADY ) printf( "unsteady problems\n" );
	else						 printf( "steady   problems\n" );

	if( ppa->CFL<0.0 )
	{
		printf( "global time step, time step is %f\n", -ppa->CFL );
	}
	else
	{
		if( ppa->isteady==STEADY ) printf( "local  time step, " );
		else					   printf( "minmum time step, " );

		printf( "CFL number is %f\n", ppa->CFL );
	}

	if( ppa->isteady==UNSTEADY ) ppa->ImResSmoStp = 0;
	if( ppa->ImResSmoStp!=0 )
	{
		printf( "Implicit Residual Smoothing Steps %d\n",  ppa->ImResSmoStp );
		printf( "Implicit Residual Smoothing Coefficient %.2f\n",  ppa->ImResSmoCoe );
	}

	if( ppa->ChgMaxRatio>0.0 )
	{
		printf( "pressure change steps %d\n",  ppa->ChgMaxStp );
		printf( "rho pre change limit %.2f\n",  ppa->ChgMaxRatio );
	}

	if( ppa->iconti==1 ) printf( "read ini data\n" );

	// DG poly order pDG, FV poly order pFV 
	if( ppa->pDG<0 ) ppa->pDG = 0;
	if( ppa->pDG>MAXDGORD ) ppa->pDG = MAXDGORD;
	if( ppa->pFV<0 ) ppa->pFV = 0;
	if( ppa->pFV>MAXFVORD ) ppa->pFV = MAXFVORD;

	if( ppa->pFV<ppa->pDG ) ppa->pDG = ppa->pFV;
	printf( "DG%d/FV%d\n", ppa->pDG,ppa->pFV );


	if( ppa->M_TVB!=0.0 ) printf( "M_TVB of shu %f\n", ppa->M_TVB );

	if( ppa->fluxtype==FLUXNND ) printf( "NND's scheme\n" );
	if( ppa->fluxtype==FLUXROE ) printf( "Roe's scheme\n" );
	if( ppa->fluxtype==FLUXLLF ) printf( "local Lax-Friedrichs' scheme\n" );
	if( ppa->fluxtype==FLUXKIN ) printf( "Kinetic's scheme\n" );
	if( ppa->fluxtype==FLUXOSH ) printf( "Osher's scheme\n" );


	printf( "*********************************************\n" );
}
