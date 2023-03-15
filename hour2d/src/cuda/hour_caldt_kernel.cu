//
// auto-generated by op2.py
//

//user function
__device__ void hour_caldt_gpu( const elem *pele, uvar *pvar, const double *CFL, double *dtmin) {
#if		EQUATIO==ADVECTI
		(*pvar).dt = param_cuda.CFL / (2.0*(*pvar).pdg+1.0) * (*pele).dh /
			( m_abs(param_cuda.a) + m_abs(param_cuda.b) + TOL );
#elif	EQUATIO==BURGERS
		double U = m_abs( (*pvar).uh[0][0][0] ) + TOL;

		(*pvar).dt = param_cuda.CFL / (2.0*(*pvar).pdg+1.0) * (*pele).dh / U /
			( m_abs(param_cuda.a) + m_abs(param_cuda.b) + TOL );
#elif	EQUATIO==EULNSEQ
		double rho = (*pvar).uh[0][0][0];

	#if		VARIABLE==CONSERV

		double pre;
		CalPre_gpu( (*pvar).uh[0][0], (*pvar).uh[0][0], pre );
		double u = (*pvar).uh[0][0][1]/rho;
		double v = (*pvar).uh[0][0][2]/rho;

	#elif	VARIABLE==PRIMITI
		double pre = (*pvar).uh[0][0][NEQ-1];
		double U = sqrt(
			(*pvar).uh[0][0][1]*(*pvar).uh[0][0][1] +
			(*pvar).uh[0][0][2]*(*pvar).uh[0][0][2] );
	#endif

		double drey = 0.0;
		if( param_cuda.ivis==1 )
		{
			drey = 2.0 * gammaa / Pr;
			double T = pre / ( param_cuda.R * rho );
			double retot = param_cuda.Re*CalMu_gpu( T, param_cuda.c4suth );

			drey *= retot / rho;
		}

		double C = sqrt( gammaa*pre/rho );

		(*pvar).dt = *CFL * (*pele).dxy[0]*(*pele).dxy[1]/
			( (m_abs(u)+C)*(*pele).dxy[1]+ (m_abs(v)+C)*(*pele).dxy[0]);

#endif
*dtmin = m_min( 9e+99, (*pvar).dt );

}

// CUDA kernel function
__global__ void op_cuda_hour_caldt(
  const elem *__restrict arg0,
  uvar *arg1,
  const double *arg2,
  double *arg3,
  int   set_size ) {

  double arg3_l[1];
  for ( int d=0; d<1; d++ ){
    arg3_l[d]=arg3[d+blockIdx.x*1];
  }

  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    hour_caldt_gpu(arg0+n*1,
               arg1+n*1,
               arg2,
               arg3_l);
  }

  //global reductions

  for ( int d=0; d<1; d++ ){
    op_reduction<OP_MIN>(&arg3[d+blockIdx.x*1],arg3_l[d]);
  }
}


//host stub function
void op_par_loop_hour_caldt(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  double*arg2h = (double *)arg2.data;
  double*arg3h = (double *)arg3.data;
  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(0);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[0].name      = name;
  OP_kernels[0].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  hour_caldt");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    //transfer constants to GPU
    int consts_bytes = 0;
    consts_bytes += ROUND_UP(1*sizeof(double));
    reallocConstArrays(consts_bytes);
    consts_bytes = 0;
    arg2.data   = OP_consts_h + consts_bytes;
    arg2.data_d = OP_consts_d + consts_bytes;
    for ( int d=0; d<1; d++ ){
      ((double *)arg2.data)[d] = arg2h[d];
    }
    consts_bytes += ROUND_UP(1*sizeof(double));
    mvConstArraysToDevice(consts_bytes);

    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_0
      int nthread = OP_BLOCK_SIZE_0;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    //transfer global reduction data to GPU
    int maxblocks = nblocks;
    int reduct_bytes = 0;
    int reduct_size  = 0;
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    reduct_size   = MAX(reduct_size,sizeof(double));
    reallocReductArrays(reduct_bytes);
    reduct_bytes = 0;
    arg3.data   = OP_reduct_h + reduct_bytes;
    arg3.data_d = OP_reduct_d + reduct_bytes;
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        ((double *)arg3.data)[d+b*1] = arg3h[d];
      }
    }
    reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));
    mvReductArraysToDevice(reduct_bytes);

    int nshared = reduct_size*nthread;
    op_cuda_hour_caldt<<<nblocks,nthread,nshared>>>(
      (elem *) arg0.data_d,
      (uvar *) arg1.data_d,
      (double *) arg2.data_d,
      (double *) arg3.data_d,
      set->size );
    //transfer global reduction data back to CPU
    mvReductArraysToHost(reduct_bytes);
    for ( int b=0; b<maxblocks; b++ ){
      for ( int d=0; d<1; d++ ){
        arg3h[d] = MIN(arg3h[d],((double *)arg3.data)[d+b*1]);
      }
    }
    arg3.data = (char *)arg3h;
    op_mpi_reduce(&arg3,arg3h);
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[0].time     += wall_t2 - wall_t1;
  OP_kernels[0].transfer += (float)set->size * arg0.size;
  OP_kernels[0].transfer += (float)set->size * arg1.size * 2.0f;
}
