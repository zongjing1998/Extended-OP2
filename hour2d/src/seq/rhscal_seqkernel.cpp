//
// auto-generated by op2.py
//

//user function
#include "../rhscal.h"

// host stub function
void op_par_loop_rhscal(char const *name, op_set set,
  op_arg arg0,
  op_arg arg4,
  op_arg arg5,
  op_arg arg6,
  op_arg arg7){

  int nargs = 8;
  op_arg args[8];

  arg0.idx = 0;
  args[0] = arg0;
  for ( int v=1; v<4; v++ ){
    args[0 + v] = op_arg_dat(arg0.dat, v, arg0.map, 1, "side", OP_READ);
  }

  args[4] = arg4;
  args[5] = arg5;
  args[6] = arg6;
  args[7] = arg7;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(10);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: rhscal\n");
  }

  int set_size = op_mpi_halo_exchanges_grouped(set, nargs, args, 1);

  if (set_size > 0) {

    for ( int n=0; n<set_size; n++ ){
      if (n==set->core_size) {
        op_mpi_wait_all_grouped(nargs, args, 1);
      }
      int map0idx;
      int map1idx;
      int map2idx;
      int map3idx;
      map0idx = arg0.map_data[n * arg0.map->dim + 0];
      map1idx = arg0.map_data[n * arg0.map->dim + 1];
      map2idx = arg0.map_data[n * arg0.map->dim + 2];
      map3idx = arg0.map_data[n * arg0.map->dim + 3];

      const side* arg0_vec[] = {
         &((side*)arg0.data)[1 * map0idx],
         &((side*)arg0.data)[1 * map1idx],
         &((side*)arg0.data)[1 * map2idx],
         &((side*)arg0.data)[1 * map3idx]};

      rhscal(
        arg0_vec,
        &((int*)arg4.data)[1 * n],
        &((elem*)arg5.data)[1 * n],
        &((double*)arg6.data)[40 * n],
        &((uvar*)arg7.data)[1 * n]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[10].name      = name;
  OP_kernels[10].count    += 1;
  OP_kernels[10].time     += wall_t2 - wall_t1;
  OP_kernels[10].transfer += (float)set->size * arg0.size;
  OP_kernels[10].transfer += (float)set->size * arg4.size;
  OP_kernels[10].transfer += (float)set->size * arg5.size;
  OP_kernels[10].transfer += (float)set->size * arg6.size;
  OP_kernels[10].transfer += (float)set->size * arg7.size * 2.0f;
  OP_kernels[10].transfer += (float)set->size * arg0.map->dim * 4.0f;
}
