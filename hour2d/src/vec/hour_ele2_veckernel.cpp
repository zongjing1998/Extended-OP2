//
// auto-generated by op2.py
//

//user function
inline void hour_ele2(const double *rhs, uvar *pvar)
{
	double tdt = (*pvar).dt;
	unsigned int ndgord = (*pvar).pdg;
	unsigned int ndgdof = DOFofOrd2d[ndgord];
	for( unsigned int i=0; i<ndgdof; i++ )
	for( unsigned int j=0; j<NEQ;    j++ )
	{
	 (*pvar).uh[0][i][j] += 0.5*tdt * rhs[i*4+j];
	}		
}


// host stub function
void op_par_loop_hour_ele2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;
  //create aligned pointers for dats
  ALIGNED_double const double * __restrict__ ptr0 = (double *) arg0.data;
  DECLARE_PTR_ALIGNED(ptr0,double_ALIGN);
  ALIGNED_uvar       uvar * __restrict__ ptr1 = (uvar *) arg1.data;
  DECLARE_PTR_ALIGNED(ptr1,uvar_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(3);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  hour_ele2");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        hour_ele2(
          &(ptr0)[40 * (n+i)],
          &(ptr1)[1 * (n+i)]);
      }
    }
    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      hour_ele2(
        &(ptr0)[40*n],
        &(ptr1)[1*n]);
    }
  }

  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[3].name      = name;
  OP_kernels[3].count    += 1;
  OP_kernels[3].time     += wall_t2 - wall_t1;
  OP_kernels[3].transfer += (float)set->size * arg0.size;
  OP_kernels[3].transfer += (float)set->size * arg1.size * 2.0f;
}
