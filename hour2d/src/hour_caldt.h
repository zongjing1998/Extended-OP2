inline void hour_caldt(const elem *pele, uvar *pvar, const double *CFL, double *dtmin)
{
#if		EQUATIO==ADVECTI
		(*pvar).dt = param.CFL / (2.0*(*pvar).pdg+1.0) * (*pele).dh / 
			( m_abs(param.a) + m_abs(param.b) + TOL );
#elif	EQUATIO==BURGERS
		double U = m_abs( (*pvar).uh[0][0][0] ) + TOL;

		(*pvar).dt = param.CFL / (2.0*(*pvar).pdg+1.0) * (*pele).dh / U / 
			( m_abs(param.a) + m_abs(param.b) + TOL );
#elif	EQUATIO==EULNSEQ
		double rho = (*pvar).uh[0][0][0];

	#if		VARIABLE==CONSERV

		double pre;
		CalPre( (*pvar).uh[0][0], (*pvar).uh[0][0], pre );
		double u = (*pvar).uh[0][0][1]/rho;
		double v = (*pvar).uh[0][0][2]/rho;


	#elif	VARIABLE==PRIMITI
		double pre = (*pvar).uh[0][0][NEQ-1];
		double U = sqrt( 
			(*pvar).uh[0][0][1]*(*pvar).uh[0][0][1] + 
			(*pvar).uh[0][0][2]*(*pvar).uh[0][0][2] );
	#endif

		double drey = 0.0;
		if( param.ivis==1 )
		{
			drey = 2.0 * gammaa / Pr;
			double T = pre / ( param.R * rho );
			double retot = param.Re*CalMu( T, param.c4suth );

			drey *= retot / rho;
		}

		double C = sqrt( gammaa*pre/rho );

		(*pvar).dt = *CFL * (*pele).dxy[0]*(*pele).dxy[1]/
			( (m_abs(u)+C)*(*pele).dxy[1]+ (m_abs(v)+C)*(*pele).dxy[0]);

#endif
*dtmin = m_min( 9e+99, (*pvar).dt );
}
