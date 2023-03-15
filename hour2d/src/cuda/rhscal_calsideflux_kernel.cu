//
// auto-generated by op2.py
//

//user function
__device__ void rhscal_calsideflux_gpu( const uvar *pvar0, const uvar *pvar1, const elem *pele0, const elem *pele1, side *psid) {

		int tsL = (*psid).lsdidx[1];

		int tsR = (*psid).lsdidx[0];
		 int  pfvL    = (*pvar0).pfv;
		 int  nfvdofL = DOFofOrd2d_cuda[pfvL];
		 int  pfvR    = (*pvar1).pfv;
		 int  nfvdofR = DOFofOrd2d_cuda[pfvR];
		 int   nTypeL = ELEMTYPE;
		 int   nTypeR = ELEMTYPE;
#ifdef QUAD
		 const int * inodeL = (*pele0).Node;
		if( inodeL[3]==inodeL[0] ) nTypeL = 3;
		 const int * inodeR = (*pele1).Node;
		if( inodeR[3]==inodeR[0] ) nTypeR = 3;
#endif
		double norm[NDIM];
		for(  int  nd=0; nd<NDIM; nd++ )
			norm[nd] = (*psid).norm[nd];
		 int  nqdpt = (*psid).nqdpt;

		double (*SdQdPtBasL)[MAXFVDOF], (*SdQdPtBasR)[MAXFVDOF];
		if( nTypeL==3 ) SdQdPtBasL = (double (*)[MAXFVDOF]) SdQdPtBasTri_cuda [tsL+3][nqdpt];
		else            SdQdPtBasL = (double (*)[MAXFVDOF]) SdQdPtBasRect_cuda[tsL+4][nqdpt];
		if( nTypeR==3 ) SdQdPtBasR = (double (*)[MAXFVDOF]) SdQdPtBasTri_cuda [tsR  ][nqdpt];
		else            SdQdPtBasR = (double (*)[MAXFVDOF]) SdQdPtBasRect_cuda[tsR  ][nqdpt];

		for(  int  qd=0; qd<nqdpt; qd++ )
		{
			double tcvL[NEQ];
			GetVar_gpu( nfvdofL, (*pvar0).wh, SdQdPtBasL[qd], tcvL );
			double tcvR[NEQ];
			GetVar_gpu( nfvdofR, (*pvar1).wh, SdQdPtBasR[qd], tcvR );
#if	EQUATIO!=EULNSEQ
			double normmodi[NDIM];
			normmodi[0] = norm[0]*param_cuda.a;
			normmodi[1] = norm[1]*param_cuda.b;
			FluxSel_gpu( param_cuda.fluxtype, tcvL, (*pvar0).wh[0],
				tcvR, (*pvar1).wh[0], normmodi, (*psid).flux[qd] );
#else
			FluxSel_gpu( param_cuda.fluxtype, tcvL, (*pvar0).wh[0],
				tcvR, (*pvar1).wh[0], norm,     (*psid).flux[qd] );
#endif
			if( param_cuda.ivis!=1 )
				continue;

			double gradcvL[NDIM][NEQ], gradcvR[NDIM][NEQ];
			for(  int  i=0; i<NDIM; i++ )
			{
				GetVar_gpu( nfvdofL, (*pvar0).grad[i], SdQdPtBasL[qd], gradcvL[i] );

				GetVar_gpu( nfvdofR, (*pvar1).grad[i], SdQdPtBasR[qd], gradcvR[i] );
			}

			double fluxVL[NDIM][NEQ], fluxVR[NDIM][NEQ];
			CalFluxV_gpu( tcvL, (*pvar0).wh[0], gradcvL, &param_cuda, fluxVL );
			CalFluxV_gpu( tcvR, (*pvar1).wh[0], gradcvR, &param_cuda, fluxVR );

			for(  int  j=0; j<NEQ; j++ )
				(*psid).flux[qd][j] -= 0.5*(
				(fluxVL[0][j]+fluxVR[0][j])*norm[0] +
				(fluxVL[1][j]+fluxVR[1][j])*norm[1] );
		}

}

// CUDA kernel function
__global__ void op_cuda_rhscal_calsideflux(
  const uvar *__restrict ind_arg0,
  const elem *__restrict ind_arg1,
  side *__restrict ind_arg2,
  const int *__restrict opDat0Map,
  const int *__restrict opDat4Map,
  int start,
  int end,
  int *col_reord,
  int   set_size) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid + start < end) {
    int n = col_reord[tid + start];
    //initialise local variables
    int map0idx;
    int map1idx;
    int map4idx;
    map0idx = opDat0Map[n + set_size * 0];
    map1idx = opDat0Map[n + set_size * 1];
    map4idx = opDat4Map[n + set_size * 0];

    //user-supplied kernel call
    rhscal_calsideflux_gpu(ind_arg0+map0idx*1,
                       ind_arg0+map1idx*1,
                       ind_arg1+map0idx*1,
                       ind_arg1+map1idx*1,
                       ind_arg2+map4idx*1);
  }
}


//host stub function
void op_par_loop_rhscal_calsideflux(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(8);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[8].name      = name;
  OP_kernels[8].count    += 1;


  int    ninds   = 3;
  int    inds[5] = {0,0,1,1,2};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: rhscal_calsideflux\n");
  }

  //get plan
  #ifdef OP_PART_SIZE_8
    int part_size = OP_PART_SIZE_8;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 2);
  if (set_size > 0) {

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);

    //execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==Plan->ncolors_core) {
        op_mpi_wait_all_grouped(nargs, args, 2);
      }
      #ifdef OP_BLOCK_SIZE_8
      int nthread = OP_BLOCK_SIZE_8;
      #else
      int nthread = OP_block_size;
      #endif

      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];
      int nblocks = (end - start - 1)/nthread + 1;
      op_cuda_rhscal_calsideflux<<<nblocks,nthread>>>(
      (uvar *)arg0.data_d,
      (elem *)arg2.data_d,
      (side *)arg4.data_d,
      arg0.map_data_d,
      arg4.map_data_d,
      start,
      end,
      Plan->col_reord,
      set->size+set->exec_size);

    }
    OP_kernels[8].transfer  += Plan->transfer;
    OP_kernels[8].transfer2 += Plan->transfer2;
  }
  op_mpi_set_dirtybit_cuda(nargs, args);
  cutilSafeCall(cudaDeviceSynchronize());
  //update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[8].time     += wall_t2 - wall_t1;
}
