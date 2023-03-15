//
// auto-generated by op2.py
//

__constant__ int direct_hour_ele2_stride_OP2CONSTANT;
int direct_hour_ele2_stride_OP2HOST=-1;
//user function
__device__ void hour_ele2_gpu( const double *rhs, uvar *pvar) {
	double tdt = (*pvar).dt;
	unsigned int ndgord = (*pvar).pdg;
	unsigned int ndgdof = DOFofOrd2d_cuda[ndgord];
	for( unsigned int i=0; i<ndgdof; i++ )
	for( unsigned int j=0; j<NEQ;    j++ )
	{
	 (*pvar).uh[0][i][j] += 0.5*tdt * rhs[(i*4+j)*direct_hour_ele2_stride_OP2CONSTANT];
	}

}

// CUDA kernel function
__global__ void op_cuda_hour_ele2(
  const double *__restrict arg0,
  uvar *arg1,
  int   set_size ) {


  //process set elements
  for ( int n=threadIdx.x+blockIdx.x*blockDim.x; n<set_size; n+=blockDim.x*gridDim.x ){

    //user-supplied kernel call
    hour_ele2_gpu(arg0+n,
              arg1+n*1);
  }
}


//host stub function
void op_par_loop_hour_ele2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1){

  int nargs = 2;
  op_arg args[2];

  args[0] = arg0;
  args[1] = arg1;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(3);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[3].name      = name;
  OP_kernels[3].count    += 1;


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  hour_ele2");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    if ((OP_kernels[3].count==1) || (direct_hour_ele2_stride_OP2HOST != getSetSizeFromOpArg(&arg0))) {
      direct_hour_ele2_stride_OP2HOST = getSetSizeFromOpArg(&arg0);
      cudaMemcpyToSymbol(direct_hour_ele2_stride_OP2CONSTANT,&direct_hour_ele2_stride_OP2HOST,sizeof(int));
    }
    //set CUDA execution parameters
    #ifdef OP_BLOCK_SIZE_3
      int nthread = OP_BLOCK_SIZE_3;
    #else
      int nthread = OP_block_size;
    #endif

    int nblocks = 200;

    op_cuda_hour_ele2<<<nblocks,nthread>>>(
      (double *) arg0.data_d,
      (uvar *) arg1.data_d,
      set->size );
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[3].time     += wall_t2 - wall_t1;
  OP_kernels[3].transfer += (float)set->size * arg0.size;
  OP_kernels[3].transfer += (float)set->size * arg1.size * 2.0f;
}
