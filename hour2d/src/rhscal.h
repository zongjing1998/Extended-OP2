inline void rhscal(const side **psid, const int *aele, const elem *pele, double *rhs, uvar *pvar)
{
	int pdg = (*pvar).pdg;
	int ndgdof = DOFofOrd2d[pdg];
	// clear rhs
	//		memset( (*pvar).rhs, 0, sizeof(double) *ndgdof*NEQ );

	for (int j = 0; j < ndgdof; j++)
	{
		for (int k = 0; k < 4; k++)
		{
			rhs[j * 4 + k] = 0;
		}
	}

	int nType = ELEMTYPE;
#ifdef QUAD
	const int *inode = (*pele).Node;
	if (inode[3] == inode[0])
		nType = 3;
#endif
	double(*SdQdBas)[MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVDOF];
	if (nType == 3)
		SdQdBas = (double(*)[MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVDOF]) SdQdPtBasTri;
	else
		SdQdBas = (double(*)[MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVDOF]) SdQdPtBasRect;
	for (int ts = 0; ts < nType; ts++)
	{
		//			 int  s =(*pele).Side[ts];
		int nqdpt = psid[ts][0].nqdpt;

		double dldds, (*Bas)[MAXFVDOF];
		if (((int)*aele) == psid[ts][0].ofElem[1])
		{
			dldds = -psid[ts][0].dl * (*pele).volr;
			Bas = (double(*)[MAXFVDOF])(&SdQdBas[ts + nType][nqdpt][0][0]);
		}
		else
		{
			dldds = psid[ts][0].dl * (*pele).volr;
			Bas = (double(*)[MAXFVDOF])(&SdQdBas[ts][nqdpt][0][0]);
		}

		for (int qd = 0; qd < nqdpt; qd++)
		{
			double coe = QdPtCoe1d[nqdpt][qd][1] * dldds;

			double fluxcoe[NEQ];
			for (int j = 0; j < NEQ; j++)
				fluxcoe[j] = psid[ts][0].flux[qd][j] * coe;

			for (int i = 0; i < ndgdof; i++)
				for (int j = 0; j < NEQ; j++)
				{
					rhs[i * 4 + j] += fluxcoe[j] * Bas[qd][i];
					//					(*pvar).rhs[i][j] += fluxcoe[j] * Bas[qd][i];
				}
		}
	}

	//		if( pdg<=0 ) goto end;

	int pfv = (*pvar).pfv;
	int nfvdof = DOFofOrd2d[pfv];
	int nqdord = pfv + pdg - 1;
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
	else
	{
		nqdpt = QdNOrd2NPt1d[nqdord];
		EmQdCoe = (double(*)[4])QdPtCoeRect[nqdpt];
		EmQdBas = (double(*)[MAXFVDOF])QdPtBasRect[nqdpt];
		EmQdBasDxi = (double(*)[MAXFVDOF])QdPtBasDxiRect[nqdpt];
		EmQdBasDet = (double(*)[MAXFVDOF])QdPtBasDetRect[nqdpt];
		nqdpt *= nqdpt;
	}
	double xix = (*pele).LocdGlb[0][0];
	double xiy = (*pele).LocdGlb[0][1];
	double etx = (*pele).LocdGlb[1][0];
	double ety = (*pele).LocdGlb[1][1];

	for (int qd = 0; qd < nqdpt; qd++)
	{

		double tcv[NEQ];
		GetVar(nfvdof, (*pvar).wh, EmQdBas[qd], tcv);
		double flux[NDIM][NEQ];
#if EQUATIO == ADVECTI
		for (int j = 0; j < NEQ; j++)
		{
			flux[0][j] = param.a * tcv[j];
			flux[1][j] = param.b * tcv[j];
		}
#elif EQUATIO == BURGERS
		for (int j = 0; j < NEQ; j++)
		{
			double tmp = 0.5 * tcv[j] * tcv[j];
			flux[0][j] = tmp * param.a;
			flux[1][j] = tmp * param.b;
		}
#elif EQUATIO == EULNSEQ
#if VARIABLE == CONSERV
		CV2FluxI(tcv, (*pvar).wh[0], flux);
#elif VARIABLE == PRIMITI
		PV2FluxI(tcv, flux);
#endif
#endif
		if (param.ivis == 1)
		{
			double graduh[NDIM][NEQ];
			for (int i = 0; i < NDIM; i++)
				GetVar(nfvdof, (*pvar).grad[i], EmQdBas[qd], graduh[i]);
			double fluxV[NDIM][NEQ];
			CalFluxV(tcv, (*pvar).wh[0], graduh, &param, fluxV);

			for (int i = 0; i < NDIM; i++)
				for (int j = 0; j < NEQ; j++)
					flux[i][j] -= fluxV[i][j];
		}

		double coe = EmQdCoe[qd][3];
		for (int i = 1; i < ndgdof; i++)
		{
			double basxcoe = (EmQdBasDxi[qd][i] * xix + EmQdBasDet[qd][i] * etx) * coe;
			double basycoe = (EmQdBasDxi[qd][i] * xiy + EmQdBasDet[qd][i] * ety) * coe;

			for (int j = 0; j < NEQ; j++)
			{
				rhs[i * 4 + j] += (flux[0][j] * basxcoe + flux[1][j] * basycoe);
				//(*pvar).rhs[i][j] += (flux[0][j]*basxcoe + flux[1][j]*basycoe);
			}
		}
	}

	// end:
}
