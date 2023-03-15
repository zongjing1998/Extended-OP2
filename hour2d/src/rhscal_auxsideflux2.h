inline void rhscal_auxsideflux2(side *psid, side *pbsesid, const uvar *pvar0, const uvar *pvar1, const elem *pele0, const elem *pele1)
{	
		int bceR = (*psid).ofElem[0];
		if( bceR >= 0 )
		{
			printf( "calauxsideflux bc error!\n" );
	//		getchar();    exit(1);
		}
		 int  nbctype = GetBCType( bceR );
//		int eL  = (*psid).ofElem[1];
		int tsL = (*psid).lsdidx[1];	
		 int  nfvordL = (*pvar0).pfv;
		 int  nfvdofL = DOFofOrd2d[nfvordL];
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
		// each gauss quad pts
		for(  int  qd = 0; qd < nqdpt; ++qd )
		{
			double tcvL[NEQ];
			GetVar( nfvdofL, (*pvar0).wh, SdQdPtBasL[qd], tcvL );
			double tcvR[NEQ];
			if( ( nbctype==BC_SOLIDSURFACE )||( nbctype==BC_GENERIC1 ) )
			{
				bcvisSolidSurf( tcvL, (*pvar0).wh[0], &param, tcvR );
			}
			else if( nbctype==BC_FARFIELD )
			{
				bcFarField( tcvL, (*pvar0).wh[0], norm, &param, tcvR );
			}
			else if( nbctype==BC_SYMMETRY )
			{
				bcSymmetry( tcvL, (*pvar0).wh[0], norm, tcvR );
			}
			else if( nbctype==BC_INTERBLK )
			{
//				 int  eR = GetElement( bceR );
				int tsR = (*psid).lsdidx[0];	
				 int  nfvordR = (*pvar1).pfv;
				 int  nfvdofR = DOFofOrd2d[nfvordR];
				 int  qdposR = qd;
				 int   nTypeR = ELEMTYPE;
#ifdef QUAD
				 const int * inodeR = (*pele1).Node;
				if( inodeR[3]==inodeR[0] ) nTypeR = 3;
#endif
				if( nTypeR==3 ) 
					GetVar( nfvdofR, (*pvar1).wh, SdQdPtBasTri [tsR][nqdpt][qdposR], tcvR );
				else // rect
					GetVar( nfvdofR, (*pvar1).wh, SdQdPtBasRect[tsR][nqdpt][qdposR], tcvR );

				double* tnbflux = (*pbsesid).flux[nqdpt-1-qd];
#if		AUXVARIABLE==CONSERV
				for(  int  j=0; j<NEQ; j++ )
				{
					(*psid).flux[qd][j] = 0.5*(tcvL[j] + tcvR[j]);
					tnbflux[j] = (*psid).flux[qd][j];
				}
#elif	AUXVARIABLE==PRIMITI
				double tpvL[NEQ];
				CV2PV( tcvL, (*pvar0).wh[0], tpvL );
				double tpvR[NEQ];
				CV2PV( tcvR, (*pvar1).wh[0], tpvR );
				tpvL[NEQ-1] /= tpvL[0]*param.R;
				tpvR[NEQ-1] /= tpvR[0]*param.R;
				for(  int  j=0; j<NEQ; j++ )
				{
					(*psid).flux[qd][j] = 0.5*(tpvL[j] + tpvR[j]);
					tnbflux[j] = (*psid).flux[qd][j];
				}
#endif
				continue;// next quad pts
			}

#if		AUXVARIABLE==CONSERV
			for(  int  j=0; j<NEQ; j++ )
			{
				(*psid).flux[qd][j] = tcvR[j];
			}
#elif	AUXVARIABLE==PRIMITI
			double tpvR[NEQ];
			CV2PV( tcvR, tcvR, tpvR );
			tpvR[NEQ-1] /= tpvR[0]*param.R;
			for(  int  j=0; j<NEQ; j++ )
			{
				(*psid).flux[qd][j] = tpvR[j];
			}
#endif
			continue;// next quad pts
		}// end of no pts
}
