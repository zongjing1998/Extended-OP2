inline void hour_ele(const double *rhs, uvar *pvar)
{
	double tdt = (*pvar).dt;
	unsigned int ndgord = (*pvar).pdg;
	unsigned int ndgdof = DOFofOrd2d[ndgord];
	for( unsigned int i=0; i<ndgdof; i++ )
	for( unsigned int j=0; j<NEQ;    j++ )
	{
	 double tmprhs = tdt*rhs[i*4+j];
	 (*pvar).uh[1][i][j] = (*pvar).uh[0][i][j] + tmprhs;
	 (*pvar).uh[0][i][j] += 0.5*tmprhs;
	} 		
}

