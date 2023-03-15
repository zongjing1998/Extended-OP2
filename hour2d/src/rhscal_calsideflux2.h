inline void rhscal_calsideflux2(const uvar *pvar0, const uvar *pvar1, const elem *pele0, const elem *pele1, side *pbsesid, side *psid)
{
		int bceR = (*psid).ofElem[0];
		 int  nbctype = GetBCType( bceR );
//		int eL  = (*psid).ofElem[1];
		int tsL = (*psid).lsdidx[1];
		 int  pfvL    = (*pvar0).pfv;
		 int  nfvdofL = DOFofOrd2d[pfvL];
		 int   nTypeL = ELEMTYPE;
#ifdef QUAD
		 const int * inodeL = (*pele0).Node;
		if( inodeL[3]==inodeL[0] ) nTypeL = 3;
#endif
		double norm[NDIM];
		for(  int  nd=0; nd<NDIM; nd++ )
			norm[nd] = (*psid).norm[nd];

		 int  nqdpt = (*psid).nqdpt;
		double (*SdQdPtBasL)[MAXFVDOF];
		if( nTypeL==3 ) SdQdPtBasL = (double (*)[MAXFVDOF]) SdQdPtBasTri [tsL+3][nqdpt];
		else            SdQdPtBasL = (double (*)[MAXFVDOF]) SdQdPtBasRect[tsL+4][nqdpt];
		switch( nbctype )
		{
		case BC_GENERIC1:
		case BC_SOLIDSURFACE:
			for(  int  qd = 0; qd < nqdpt; ++qd ) // each gauss quad pts
			{
				double tcvL[NEQ];
				GetVar( nfvdofL, (*pvar0).wh, SdQdPtBasL[qd], tcvL );
				double tcvR[NEQ];
				//bcSolidSurf3( tcvL, (*pvar0).wh[0], norm, tcvR );
				FluxSel( param.fluxtype, tcvL, (*pvar0).wh[0], 
					tcvR, (*pvar0).wh[0], norm, (*psid).flux[qd] );
				if( param.ivis==0 )
					bcSolidSurf( tcvL, (*pvar0).wh[0], norm, tcvR );
				else
				{
					bcvisSolidSurf ( tcvL, (*pvar0).wh[0], &param, tcvR );
				}
				CV2FluxIN2d( tcvR, tcvR, norm, (*psid).flux[qd] );
				if( param.ivis!=1 ) continue;
				double gradcvL[NDIM][NEQ];
				for(  int  i=0; i<NDIM; i++ )
					GetVar( nfvdofL, (*pvar0).grad[i], SdQdPtBasL[qd], gradcvL[i] );
				double fluxV[NDIM][NEQ];
				CalFluxV( tcvR, tcvR, gradcvL, &param, fluxV );
				double fluxVN[NEQ];
				for(  int  j=0; j<NEQ; j++ )
					fluxVN[j] = fluxV[0][j]*norm[0]+fluxV[1][j]*norm[1];
#if	EQUATIO==EULNSEQ
				if( param.bcTtype==ADIAWALL )
					fluxVN[3] = 0.0;
#endif
				for(  int  j=0; j<NEQ; j++ )
					(*psid).flux[qd][j] -= fluxVN[j];
			}
			break;

		case BC_FARFIELD:
			for(  int  qd=0; qd<nqdpt; qd++ )
			{
				double tcvL[NEQ];
				GetVar( nfvdofL, (*pvar0).wh, SdQdPtBasL[qd], tcvL );
				double tcvR[NEQ];
				bcFarField( tcvL, (*pvar0).wh[0], norm, &param, tcvR );
				CV2FluxIN2d( tcvR, tcvR, norm, (*psid).flux[qd] );	
				if( param.ivis!=1 ) continue;
				double gradcvL[NDIM][NEQ];
				for(  int  i=0; i<NDIM; i++ )
					GetVar( nfvdofL, (*pvar0).grad[i], SdQdPtBasL[qd], gradcvL[i] );
				double fluxV[NDIM][NEQ];
				CalFluxV( tcvR, tcvR, gradcvL, &param, fluxV );
				for(  int  j=0; j<NEQ; j++ )
					(*psid).flux[qd][j] -= fluxV[0][j]*norm[0]+fluxV[1][j]*norm[1];
			}
			break;
		case BC_SYMMETRY:
			for(  int  qd=0; qd<nqdpt; qd++ )
			{
				double tcvL[NEQ];
				GetVar( nfvdofL, (*pvar0).wh, SdQdPtBasL[qd], tcvL );
				double tcvR[NEQ];
				bcSymmetry( tcvL, (*pvar0).wh[0], norm, tcvR );
				CV2FluxIN2d( tcvR, tcvR, norm, (*psid).flux[qd] );
				//bcSymmetry2( tcvL, (*pvar0).wh[0], norm, tcvR );

				if( param.ivis!=1 ) continue;
				double gradcvL[NDIM][NEQ];
				for(  int  i=0; i<NDIM; i++ )
					GetVar( nfvdofL, (*pvar0).grad[i], SdQdPtBasL[qd], gradcvL[i] );
				bcSymmetry2( tcvL, (*pvar0).wh[0], norm, tcvR );
				double gradcvR[NDIM][NEQ];
				double nx = norm[0];
				double ny = norm[1];
				for(  int  j=0; j<NEQ; j++ )
				{
					double gradcvRn =  nx*gradcvR[0][j] + ny*gradcvR[1][j];
					double gradcvRt = -ny*gradcvR[0][j] + nx*gradcvR[1][j];
					gradcvR[0][j] = -nx*gradcvRn - ny*gradcvRt;
					gradcvR[1][j] = -ny*gradcvRn + nx*gradcvRt;
				}
				double fluxV[NDIM][NEQ];
				CalFluxV( tcvL, (*pvar0).wh[0], gradcvL, &param, fluxV );

				for(  int  j=0; j<NEQ; j++ )
					(*psid).flux[qd][j] -= 0.5*(fluxV[0][j]*norm[0]+fluxV[1][j]*norm[1]);

				double tcvR0[NEQ];
				bcSymmetry2( (*pvar0).wh[0], (*pvar0).wh[0], norm, tcvR0 );
				CalFluxV( tcvR, tcvR0, gradcvR, &param, fluxV );
				for(  int  j=0; j<NEQ; j++ )
					(*psid).flux[qd][j] -= 0.5*(fluxV[0][j]*norm[0]+fluxV[1][j]*norm[1]);
			}
			break;
		case BC_INTERBLK:
//			 int  eR = GetElement( bceR );
			int tsR = (*psid).lsdidx[0];	
			 int  pfvR    = (*pvar1).pfv;
			 int  nfvdofR = DOFofOrd2d[pfvR];
			 int   nTypeR = ELEMTYPE;
#ifdef QUAD
			 const int * inodeR = (*pele1).Node;
			if( inodeR[3]==inodeR[0] ) nTypeR = 3;
#endif
			double (*SdQdPtBasR)[MAXFVDOF];
			if( nTypeR==3 ) SdQdPtBasR = (double (*)[MAXFVDOF]) SdQdPtBasTri [tsR][nqdpt];
			else            SdQdPtBasR = (double (*)[MAXFVDOF]) SdQdPtBasRect[tsR][nqdpt];
			for(  int  qd=0; qd<nqdpt; qd++ )
			{
				double tcvL[NEQ];
				GetVar( nfvdofL, (*pvar0).wh, SdQdPtBasL[qd], tcvL );
				double tcvR[NEQ];
				GetVar( nfvdofR, (*pvar1).wh, SdQdPtBasR[qd], tcvR );
double* tnbflux = (*pbsesid).flux[nqdpt-1-qd];
#if		EQUATIO!=EULNSEQ
				double normmodi[NDIM];
				normmodi[0] = norm[0]*param.a;
				normmodi[1] = norm[1]*param.b;
				FluxSel( param->fluxtype, tcvL, (*pvar0).wh[0], 
					tcvR, (*pvar1).wh[0], normmodi, (*psid).flux[qd] );
#elif	EQUATIO==EULNSEQ
				FluxSel( param.fluxtype, tcvL, (*pvar0).wh[0], 
					tcvR, (*pvar1).wh[0], norm, (*psid).flux[qd] );
#endif
				for(  int  j=0; j<NEQ; j++ )
					tnbflux[j] = -(*psid).flux[qd][j];
				if( param.ivis!=1 ) continue;
				double gradcvL[NDIM][NEQ], gradcvR[NDIM][NEQ];
				for(  int  i=0; i<NDIM; i++ )
				{
					GetVar( nfvdofL, (*pvar0).grad[i], SdQdPtBasL[qd], gradcvL[i] );
					GetVar( nfvdofR, (*pvar1).grad[i], SdQdPtBasR[qd], gradcvR[i] );
				}
				double fluxV[NDIM][NEQ];
				CalFluxV( tcvL, (*pvar0).wh[0], gradcvL, &param, fluxV );
				for(  int  j=0; j<NEQ; j++ )
					(*psid).flux[qd][j] -= 0.5*(fluxV[0][j]*norm[0]+fluxV[1][j]*norm[1]);

				CalFluxV( tcvR, (*pvar1).wh[0], gradcvR, &param, fluxV );
				for(  int  j=0; j<NEQ; j++ )
					(*psid).flux[qd][j] -= 0.5*(fluxV[0][j]*norm[0]+fluxV[1][j]*norm[1]);

				for(  int  j=0; j<NEQ; j++ )
					tnbflux[j] = -(*psid).flux[qd][j];
			}// end of interblk bc
			break;

		}// end of nbctype
}
