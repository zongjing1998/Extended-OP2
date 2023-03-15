inline void rhscal_calsideflux(const uvar *pvar0, const uvar *pvar1, const elem *pele0, const elem *pele1, side *psid)
{
//		int eL  = (*psid).ofElem[1];
		int tsL = (*psid).lsdidx[1];	
//		int eR  = (*psid).ofElem[0];
		int tsR = (*psid).lsdidx[0];	
		 int  pfvL    = (*pvar0).pfv;
		 int  nfvdofL = DOFofOrd2d[pfvL];
		 int  pfvR    = (*pvar1).pfv;
		 int  nfvdofR = DOFofOrd2d[pfvR];
		 int   nTypeL = ELEMTYPE;
		 int   nTypeR = ELEMTYPE;
#ifdef QUAD
		 const int * inodeL = (*pele0).Node;
		if( inodeL[3]==inodeL[0] ) nTypeL = 3;
		 const int * inodeR = (*pele1).Node;
		if( inodeR[3]==inodeR[0] ) nTypeR = 3;
#endif
		double norm[NDIM];
		for(  int  nd=0; nd<NDIM; nd++ )
			norm[nd] = (*psid).norm[nd];
		 int  nqdpt = (*psid).nqdpt;

		double (*SdQdPtBasL)[MAXFVDOF], (*SdQdPtBasR)[MAXFVDOF];
		if( nTypeL==3 ) SdQdPtBasL = (double (*)[MAXFVDOF]) SdQdPtBasTri [tsL+3][nqdpt];
		else            SdQdPtBasL = (double (*)[MAXFVDOF]) SdQdPtBasRect[tsL+4][nqdpt];
		if( nTypeR==3 ) SdQdPtBasR = (double (*)[MAXFVDOF]) SdQdPtBasTri [tsR  ][nqdpt];
		else            SdQdPtBasR = (double (*)[MAXFVDOF]) SdQdPtBasRect[tsR  ][nqdpt];
		// each gauss quad pts
		for(  int  qd=0; qd<nqdpt; qd++ )
		{
			double tcvL[NEQ];
			GetVar( nfvdofL, (*pvar0).wh, SdQdPtBasL[qd], tcvL );
			double tcvR[NEQ];
			GetVar( nfvdofR, (*pvar1).wh, SdQdPtBasR[qd], tcvR );
#if	EQUATIO!=EULNSEQ
			double normmodi[NDIM];
			normmodi[0] = norm[0]*param.a;
			normmodi[1] = norm[1]*param.b;
			FluxSel( param.fluxtype, tcvL, (*pvar0).wh[0], 
				tcvR, (*pvar1).wh[0], normmodi, (*psid).flux[qd] );
#else
			FluxSel( param.fluxtype, tcvL, (*pvar0).wh[0], 
				tcvR, (*pvar1).wh[0], norm,     (*psid).flux[qd] );
#endif
			if( param.ivis!=1 )
				continue;

			double gradcvL[NDIM][NEQ], gradcvR[NDIM][NEQ];
			for(  int  i=0; i<NDIM; i++ )
			{
				GetVar( nfvdofL, (*pvar0).grad[i], SdQdPtBasL[qd], gradcvL[i] );
				//GetVar( nfvdofL, (*pvar0).grad[1], SdQdPtBasL[qd], gradcvL[1] );
				GetVar( nfvdofR, (*pvar1).grad[i], SdQdPtBasR[qd], gradcvR[i] );
			}

			double fluxVL[NDIM][NEQ], fluxVR[NDIM][NEQ];
			CalFluxV( tcvL, (*pvar0).wh[0], gradcvL, &param, fluxVL );
			CalFluxV( tcvR, (*pvar1).wh[0], gradcvR, &param, fluxVR );

			for(  int  j=0; j<NEQ; j++ )
				(*psid).flux[qd][j] -= 0.5*( 
				(fluxVL[0][j]+fluxVR[0][j])*norm[0] + 
				(fluxVL[1][j]+fluxVR[1][j])*norm[1] );
		}// end of qd
}
