inline void rhscal_auxsideflux(const uvar *pvar0, const uvar *pvar1, const elem *pele0, const elem *pele1, side *psid)
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
		 int  nqdpt = (*psid).nqdpt;
		double (*SdQdPtBasL)[MAXFVDOF], (*SdQdPtBasR)[MAXFVDOF];
		if( nTypeL==3 ) SdQdPtBasL = (double (*)[MAXFVDOF]) SdQdPtBasTri [tsL+3][nqdpt];
		else            SdQdPtBasL = (double (*)[MAXFVDOF]) SdQdPtBasRect[tsL+4][nqdpt];
		if( nTypeR==3 ) SdQdPtBasR = (double (*)[MAXFVDOF]) SdQdPtBasTri [tsR  ][nqdpt];
		else            SdQdPtBasR = (double (*)[MAXFVDOF]) SdQdPtBasRect[tsR  ][nqdpt];

		// each gauss quad pts
		for(  int  qd = 0; qd < nqdpt; ++qd )
		{
			double tcvL[NEQ];
			GetVar( nfvdofL, (*pvar0).wh, SdQdPtBasL[qd], tcvL );
			double tcvR[NEQ];
			GetVar( nfvdofR, (*pvar1).wh, SdQdPtBasR[qd], tcvR );

#if		AUXVARIABLE==CONSERV
			for(  int  j=0; j<NEQ; j++ )
			{
				(*psid).flux[qd][j] = 0.5*(tcvL[j] + tcvR[j]);
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
			}
#endif
		}// end of qd
}
