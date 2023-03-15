//
// auto-generated by op2.py
//

//user function
#include "../rhscal_auxsideflux2.h"

// host stub function
void op_par_loop_rhscal_auxsideflux2(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4,
  op_arg arg5){

  int nargs = 6;
  op_arg args[6];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  args[5] = arg5;

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(6);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: rhscal_auxsideflux2\n");
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
      map1idx = arg1.map_data[n * arg1.map->dim + 0];
      map2idx = arg2.map_data[n * arg2.map->dim + 0];
      map3idx = arg2.map_data[n * arg2.map->dim + 1];


      rhscal_auxsideflux2(
        &((side*)arg0.data)[1 * map0idx],
        &((side*)arg1.data)[1 * map1idx],
        &((uvar*)arg2.data)[1 * map2idx],
        &((uvar*)arg2.data)[1 * map3idx],
        &((elem*)arg4.data)[1 * map2idx],
        &((elem*)arg4.data)[1 * map3idx]);
    }
  }

  if (set_size == 0 || set_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[6].name      = name;
  OP_kernels[6].count    += 1;
  OP_kernels[6].time     += wall_t2 - wall_t1;
  OP_kernels[6].transfer += (float)set->size * arg0.size * 2.0f;
  OP_kernels[6].transfer += (float)set->size * arg2.size;
  OP_kernels[6].transfer += (float)set->size * arg4.size;
  OP_kernels[6].transfer += (float)set->size * arg0.map->dim * 4.0f;
  OP_kernels[6].transfer += (float)set->size * arg1.map->dim * 4.0f;
  OP_kernels[6].transfer += (float)set->size * arg2.map->dim * 4.0f;
}