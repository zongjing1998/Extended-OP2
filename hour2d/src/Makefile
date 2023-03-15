include make-common.inc

OP2_INC      = -I$(OP2_INSTALL_PATH)/c/include
OP2_LIB      = -L$(OP2_INSTALL_PATH)/c/lib
CC  = g++ -cpp -O2 -fopenmp -w
#CC  = icpc -cpp -O2 -qopenmp -w
CUDA_INC     = -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB     = -L$(CUDA_INSTALL_PATH)/lib64
NVCCFLAGS    = -O3 -arch=sm_70 -Xptxas=-v -Dlcm=ca  -use_fast_math -maxrregcount=256
CFLAGS += -DOP_USER_DATATYPES
all: clean hour2d_seq  hour2d_cuda hour2d_genseq hour2d_openmp
# all: clean hour2d_openmp

hour2d_cuda: hour2d_kernels_cu.o  hour2d_cuda1
hour2d_cuda1:
	$(CC) $(CFLAGS) $(OP2_INC) $(OP2_LIB) $(CUDA_INC) $(CUDA_LIB) \
         bc2d.cpp  flux2d.cpp  gridio.cpp  gridparam.cpp  ic.cpp  misc.cpp  poly1d.cpp  poly2d.cpp  postproc.cpp  rcm.cpp \
	 recons.cpp  rhs_op.cpp  solio.cpp  solparam.cpp \
	 hour2d_op.cpp hour2d_kernels_cu.o -lcudart -lop2_cuda -o hour2d_cuda

hour2d_kernels_cu.o:
		nvcc  $(NVCCFLAGS) $(OP2_INC) -c -o hour2d_kernels_cu.o cuda/hour2d_kernels.cu	
hour2d_openmp:
	$(CC) $(CFLAGS) $(OP2_INC) $(OP2_LIB) \
         bc2d.cpp  flux2d.cpp  gridio.cpp  gridparam.cpp  ic.cpp  misc.cpp  poly1d.cpp \
	 poly2d.cpp  postproc.cpp  rcm.cpp  recons.cpp  rhs_op.cpp  solio.cpp  solparam.cpp  hour2d_op.cpp \
	openmp/hour2d_kernels.cpp -lm -lop2_openmp -o hour2d_openmp

hour2d_genseq:
	$(CC) $(CFLAGS) $(OP2_INC) $(OP2_LIB) \
         bc2d.cpp  flux2d.cpp  gridio.cpp  gridparam.cpp  ic.cpp  misc.cpp  poly1d.cpp \
	 poly2d.cpp  postproc.cpp  rcm.cpp  recons.cpp  rhs_op.cpp  solio.cpp  solparam.cpp  hour2d.cpp \
	seq/hour2d_seqkernels.cpp -lm -lop2_seq -o hour2d_genseq

hour2d_seq:
	$(CC) $(CFLAGS) $(OP2_INC) $(OP2_LIB) \
         bc2d.cpp  flux2d.cpp  gridio.cpp  gridparam.cpp  ic.cpp  misc.cpp  poly1d.cpp \
	 poly2d.cpp  postproc.cpp  rcm.cpp  recons.cpp  rhs.cpp  solio.cpp  solparam.cpp  hour2d.cpp \
	 -lm -lop2_seq -o hour2d_seq

# hour2d_mpi_cuda:	hour2d_mpi_op.cpp hour2d_kernels_mpi_cu.o
# 	mpicc -O2 -fopenmp -w -DOP_USER_DATATYPES  $(OP2_INC) $(OP2_LIB) $(PARMETIS_INC) $(CUDA_INC) $(CUDA_LIB) \
# 	bc2d.cpp flux2d.cpp gridio.cpp gridparam.cpp ic.cpp misc.cpp poly1d.cpp poly2d.cpp postproc.cpp rcm.cpp \
# 	recons.cpp rhs_op.cpp solio.cpp solparam.cpp \
# 	hour2d_mpi_op.cpp hour2d_kernels_mpi_cu.o -lop2_mpi_cuda -lcudart -lm $(PARMETIS_LIB) -o hour2d_mpi_cuda

# hour2d_kernels_mpi_cu.o:
# 		nvcc  $(NVCCFLAGS) $(OP2_INC) -c -o hour2d_kernels_mpi_cu.o cuda/hour2d_mpi_kernels.cu	

# hour2d_mpi_cuda:	hour2d_mpi_op.cpp hour2d_kernels_mpi_cu.o Makefile
# 		  mpicc -O2 -fopenmp -w -DOP_USER_DATATYPES  hour2d_mpi_op.cpp hour2d_kernels_mpi_cu.o \
# 		  $(OP2_INC) $(PARMETIS_INC) $(PTSCOTCH_INC) \
# 		  $(OP2_LIB) -lop2_mpi_cuda $(PARMETIS_LIB) $(PTSCOTCH_LIB) \
# 		  $(CUDA_LIB) -lcudart -o hour2d_mpi_cuda -lm

# hour2d_kernels_mpi_cu.o:
# 		  nvcc  $(INC) $(NVCCFLAGS) $(OP2_INC) -I $(MPI_INSTALL_PATH)/include \
# 		  -I. -Icuda -c -o hour2d_kernels_mpi_cu.o cuda/hour2d_mpi_kernels.cu


clean:
	rm -f hour2d_seq hour2d_genseq hour2d_openmp hour2d_cuda *.o	
		

