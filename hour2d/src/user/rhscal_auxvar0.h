inline void rhscal_auxvar(const side **psid, const int *aele, const elem *pele, uvar *pvar)
{
	double xix = (*pele).LocdGlb[0][0];
	double xiy = (*pele).LocdGlb[0][1];
	double etx = (*pele).LocdGlb[1][0];
	double ety = (*pele).LocdGlb[1][1];

	int nfvord = (*pvar).pfv;
	int nfvdof = DOFofOrd2d[nfvord];
	//		memset( (*pvar).grad, 0, sizeof(double) *NDIM*MAXFVDOF*NEQ );

	for (int i = 0; i < NDIM; i++)
		for (int j = 0; j < MAXFVDOF; j++)
			for (int k = 0; k < NEQ; k++)
				(*pvar).grad[i][j][k] = 0;

	// get number of gauss points of element gauss quad
	int nqdord = nfvord + nfvord;
	if (nqdord < 1)
		nqdord = 1;

	int nType = ELEMTYPE;
#ifdef QUAD
	const int *inode = (*pele).Node;
	if (inode[3] == inode[0])
		nType = 3;
#endif
	int nqdpt;
	double(*EmQdCoe)[4], (*EmQdBas)[MAXFVDOF], (*EmQdBasDxi)[MAXFVDOF], (*EmQdBasDet)[MAXFVDOF];

	if (nType == 3) // tri
	{
		nqdpt = SymQdTriDim[nqdord];
		EmQdCoe = (double(*)[4])SymQdTriCoe[nqdord];
		EmQdBas = (double(*)[MAXFVDOF])QdPtBasTri[nqdord];
		EmQdBasDxi = (double(*)[MAXFVDOF])QdPtBasDxiTri[nqdord];
		EmQdBasDet = (double(*)[MAXFVDOF])QdPtBasDetTri[nqdord];
	}
	else // rect
	{
		// nqdpt = (nqdord+1)/2;
		nqdpt = QdNOrd2NPt1d[nqdord];
		EmQdCoe = (double(*)[4])QdPtCoeRect[nqdpt];
		EmQdBas = (double(*)[MAXFVDOF])QdPtBasRect[nqdpt];
		EmQdBasDxi = (double(*)[MAXFVDOF])QdPtBasDxiRect[nqdpt];
		EmQdBasDet = (double(*)[MAXFVDOF])QdPtBasDetRect[nqdpt];
		nqdpt *= nqdpt;
	}

	// zongjing
	side psid_temp[4];
	// double flux_temp[4][5][4];
	for (int i = 0; i < 4; i++)
	{
		psid_temp[i].Node[0] = psid[i][0].Node[0];
		psid_temp[i].Node[1] = psid[i][0].Node[1];
		psid_temp[i].ofElem[0] = psid[i][0].ofElem[0];
		psid_temp[i].ofElem[1] = psid[i][0].ofElem[1];
		psid_temp[i].lsdidx[0] = psid[i][0].lsdidx[0];
		psid_temp[i].lsdidx[1] = psid[i][0].lsdidx[1];
		psid_temp[i].dl = psid[i][0].dl;
		for (int jj = 0; jj < NDIM; jj++)
			psid_temp[i].norm[jj] = psid[i][0].norm[jj];
		psid_temp[i].etaf = psid[i][0].etaf;
		psid_temp[i].nqdpt = psid[i][0].nqdpt;
		for (int np = 0; np < 5; np++)
		{
			for (int j = 0; j < NEQ; j++)
			{
				// flux_temp[i][np][j]=psid[i][0].flux[np][j];
				psid_temp[i].flux[np][j] = psid[i][0].flux[np][j];
			}
		}
	}
	double grad_temp[NDIM][MAXFVDOF][NEQ];
	for (int i = 0; i < NDIM; i++)
	{
		for (int j = 0; j < MAXFVDOF; j++)
		{
			for (int k = 0; k < NEQ; k++)
			{
				grad_temp[i][j][k] = (*pvar).grad[i][j][k];
			}
		}
	}
	// zongjing

	for (int qd = 0; qd < nqdpt; qd++) // rhs_CalAuxVar_1
	{
		double coe = EmQdCoe[qd][3];

		double tcv[NEQ];
		GetVar(nfvdof, (*pvar).wh, EmQdBas[qd], tcv);

#if AUXVARIABLE == PRIMITI
		double tpv[NEQ];
		CV2PV(tcv, (*pvar).wh[0], tpv);
		tpv[NEQ - 1] /= tpv[0] * (param.R);
		for (int j = 0; j < NEQ; j++)
			tcv[j] = tpv[j];
#endif
		for (int j = 0; j < NEQ; j++)
			tcv[j] *= coe;

		for (int k = 0; k < NDIM; k++)
			for (int i = 0; i < nfvdof; i++)
			{
				double basxy[NDIM];
				basxy[0] = EmQdBasDxi[qd][i] * xix + EmQdBasDet[qd][i] * etx;
				basxy[1] = EmQdBasDxi[qd][i] * xiy + EmQdBasDet[qd][i] * ety;

				for (int j = 0; j < NEQ; j++)
					grad_temp[k][i][j] -= tcv[j] * basxy[k]; // zongjing
										 // (*pvar).grad[k][i][j] -= tcv[j] * basxy[k];
			} // end of nfvdof
	}		  // element gauss quad

	for (int ts = 0; ts < nType; ts++)
	{
		//			 int  s = (*pele).Side[ts];

		double norm[NDIM];
		for (int nd = 0; nd < NDIM; nd++)
			norm[nd] = psid[ts][0].norm[nd];

		int nqdpt = psid[ts][0].nqdpt;
		double dldds = -psid[ts][0].dl * (*pele).volr;
		int sign = 1;
		int baspos = 0;

		if (((int)*aele) == psid[ts][0].ofElem[1])
		{
			dldds = -dldds;
			sign = -1;
			baspos = nqdpt - 1;
		}

		// each gauss points
		for (int np = 0; np < nqdpt; np++, baspos += sign) // rhs_CalAuxVar_2
		{
			double *tmbas;
			if (nType == 3)
				tmbas = SdQdPtBasTri[ts][nqdpt][baspos];
			else
				tmbas = SdQdPtBasRect[ts][nqdpt][baspos];
			double quadcoe = GauQd1dCoe[nqdpt][np][2] * dldds;
			for (int k = 0; k < NDIM; k++)
				for (int i = 0; i < nfvdof; i++)
				{
					double temp = quadcoe * tmbas[i];
					for (int j = 0; j < NEQ; j++)
					{
						double temp2 = temp * psid[ts][0].flux[np][j];
						//				for(  int  k=0; k<NDIM ; k++ )
						// (*pvar).grad[k][i][j] += temp2 * norm[k];//zongjing
						grad_temp[k][i][j] += temp2 * norm[k];
					} // each variables
				}
		} // each quad points
	}	  // end ts each side of an elem

	// zongjing
	for (int i = 0; i < NDIM; i++)
	{
		for (int j = 0; j < MAXFVDOF; j++)
		{
			for (int k = 0; k < NEQ; k++)
			{
				(*pvar).grad[i][j][k] = grad_temp[i][j][k];
			}
		}
	}
}
