inline void recons_0wh(uvar *pvar)
{
		unsigned int pdg    = (*pvar).pdg;
		unsigned int pfv    = (*pvar).pfv;
		unsigned int ndgdof = DOFofOrd2d[pdg];
		unsigned int nfvdof = DOFofOrd2d[pfv];
		
		for( unsigned int i=0; i<ndgdof; i++ )
		for( unsigned int j=0; j<NEQ;    j++ )
			(*pvar).wh[i][j] = (*pvar).uh[0][i][j];

		for( unsigned int i=ndgdof; i<nfvdof; i++ )
		for( unsigned int j=0; j<NEQ;         j++ )
			(*pvar).wh[i][j] = 0.0;
}



