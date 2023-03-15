//
// auto-generated by op2.py
//

//user function
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

// host stub function
void op_par_loop_hour_caldt(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  //create aligned pointers for dats
  ALIGNED_elem const elem * __restrict__ ptr0 = (elem *) arg0.data;
  DECLARE_PTR_ALIGNED(ptr0,elem_ALIGN);
  ALIGNED_uvar       uvar * __restrict__ ptr1 = (uvar *) arg1.data;
  DECLARE_PTR_ALIGNED(ptr1,uvar_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(0);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  hour_caldt");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      double dat2[SIMD_VEC];
      for ( int i=0; i<SIMD_VEC; i++ ){
        dat2[i] = *((double*)arg2.data);
      }
      double dat3[SIMD_VEC];
      for ( int i=0; i<SIMD_VEC; i++ ){
        dat3[i] = INFINITY;
      }
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        hour_caldt(
          &(ptr0)[1 * (n+i)],
          &(ptr1)[1 * (n+i)],
          &dat2[i],
          &dat3[i]);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        *(double*)arg3.data = MIN(*(double*)arg3.data,dat3[i]);
      }
    }
    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      hour_caldt(
        &(ptr0)[1*n],
        &(ptr1)[1*n],
        (double*)arg2.data,
        (double*)arg3.data);
    }
  }

  // combine reduction data
  op_mpi_reduce(&arg3,(double*)arg3.data);
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[0].name      = name;
  OP_kernels[0].count    += 1;
  OP_kernels[0].time     += wall_t2 - wall_t1;
  OP_kernels[0].transfer += (float)set->size * arg0.size;
  OP_kernels[0].transfer += (float)set->size * arg1.size * 2.0f;
}
