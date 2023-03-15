//
// auto-generated by op2.py
//

#include "user_genseq.h"
// global constants
extern const int DOFofOrd2d[16];
extern const int SymQdTriDim[7];
extern int QdNOrd2NPt1d[10];
extern const double GauQd1dCoe[6][5][3];
extern double QdPtCoe1d[(MAXSEGQDPTS + 1)][MAXSEGQDPTS][2];
extern const double SymQdTriCoe[7][12][4];
extern double QdPtBasTri[7][12][10];
extern double QdPtCoeRect[6][25][4];
extern double QdPtBasRect[6][25][10];
extern double SdQdPtBasTri[6][6][5][10];
extern double SdQdPtBasRect[8][6][5][10];
extern double QdPtBasDxiTri[7][12][10];
extern double QdPtBasDetTri[7][12][10];
extern double QdPtBasDxiRect[6][25][10];
extern double QdPtBasDetRect[6][25][10];

// header
#include "op_lib_cpp.h"

// user kernel files
#include "hour_caldt_seqkernel.cpp"
#include "hour_ele_seqkernel.cpp"
#include "recons_1wh_seqkernel.cpp"
#include "hour_ele2_seqkernel.cpp"
#include "recons_0wh_seqkernel.cpp"
#include "rhscal_auxsideflux_seqkernel.cpp"
#include "rhscal_auxsideflux2_seqkernel.cpp"
#include "rhscal_auxvar_seqkernel.cpp"
#include "rhscal_calsideflux_seqkernel.cpp"
#include "rhscal_calsideflux2_seqkernel.cpp"
#include "rhscal_seqkernel.cpp"
