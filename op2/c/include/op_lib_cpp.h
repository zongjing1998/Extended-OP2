#ifndef __OP_LIB_CPP_H
#define __OP_LIB_CPP_H

/*
 * include core definitions and C declarations of op2 user-level routines
 */

#include <op_lib_c.h>
#include <op_lib_core.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * run-time type-checking routines
 */

inline int type_error(const double *a, const char *type)
{
  (void)a;
  return (strcmp(type, "double") && strcmp(type, "double:soa"));
}
inline int type_error(const float *a, const char *type)
{
  (void)a;
  return (strcmp(type, "float") && strcmp(type, "float:soa"));
}
inline int type_error(const int *a, const char *type)
{
  (void)a;
  return (strcmp(type, "int") && strcmp(type, "int:soa"));
}
inline int type_error(const uint *a, const char *type)
{
  (void)a;
  return (strcmp(type, "uint") && strcmp(type, "uint:soa"));
}
inline int type_error(const ll *a, const char *type)
{
  (void)a;
  return (strcmp(type, "ll") && strcmp(type, "ll:soa"));
}
inline int type_error(const ull *a, const char *type)
{
  (void)a;
  return (strcmp(type, "ull") && strcmp(type, "ull:soa"));
}
inline int type_error(const bool *a, const char *type)
{
  (void)a;
  return (strcmp(type, "bool") && strcmp(type, "bool:soa"));
}

/*
 * add in user's datatypes
 */
#ifdef OP_USER_DATATYPES
#include "op_user_datatypes.h"
#endif

/*
 * zero constants
 */

#define ZERO_double 0.0;
#define ZERO_float 0.0f;
#define ZERO_int 0;
#define ZERO_uint 0;
#define ZERO_ll 0;
#define ZERO_ull 0;
#define ZERO_bool 0;

/*
 * external variables declared in op_lib_core.cpp
 */

extern int OP_diags, OP_part_size, OP_block_size, OP_gpu_direct;

extern int OP_set_index, OP_set_max, OP_map_index, OP_map_max, OP_dat_index,
    OP_plan_index, OP_plan_max, OP_kern_max;

extern op_set *OP_set_list;
extern op_map *OP_map_list;
extern Double_linked_list OP_dat_list;
extern op_kernel *OP_kernels;
extern double OP_plan_time;

op_dat op_decl_dat_char(op_set, int, char const *, int, char *, char const *);
op_dat op_decl_dat_temp_char(op_set, int, char const *, int, char const *);
int op_free_dat_temp_char(op_dat dat);

/* Implementation */

template <class T>
op_dat op_decl_dat(op_set set, int dim, char const *type, T *data,
                   char const *name)
{

  if (type_error(data, type))
  {
    printf("incorrect type specified for dataset \"%s\" \n", name);
    exit(1);
  }

  return op_decl_dat_char(set, dim, type, sizeof(T), (char *)data, name);
}

template <class T>
void op_decl_const2(char const *name, int dim, char const *type, T *data)
{
  if (type_error(data, type))
  {
    printf("incorrect type specified for constant \"%s\" \n", name);
    exit(1);
  }
  op_decl_const_char(dim, type, sizeof(T), (char *)data, name);
}
template <class T, class T1>
void op_decl_const2(char const *name, int dim, char const *type, T *data, T1 *data1)
{
  if (type_error(data, type))
  {
    printf("incorrect type specified for constant \"%s\" \n", name);
    exit(1);
  }
  op_decl_const_char(dim, type, sizeof(T1), (char *)data, name);
}

template <class T>
void op_decl_const(int dim, char const *type, T *data)
{
  (void)dim;
  if (type_error(data, type))
  {
    printf("incorrect type specified for constant in op_decl_const");
    exit(1);
  }
}
template <class T>
void op_decl_const(int dim, char const *type, T *data, int gbl_flag)
{
  (void)dim;
  if (type_error(data, type))
  {
    printf("incorrect type specified for constant in op_decl_const");
    exit(1);
  }
}

template <class T>
op_arg op_arg_gbl(T *data, int dim, char const *type, op_access acc)
{
  if (type_error(data, type))
    return op_arg_gbl_char((char *)data, dim, "error", sizeof(T), acc);
  else
    return op_arg_gbl_char((char *)data, dim, type, sizeof(T), acc);
}

template <class T>
op_arg op_opt_arg_gbl(int opt, T *data, int dim, char const *type,
                      op_access acc)
{
  if (type_error(data, type))
    return op_opt_arg_gbl_char(opt, (char *)data, dim, "error", sizeof(T), acc);
  else
    return op_opt_arg_gbl_char(opt, (char *)data, dim, type, sizeof(T), acc);
}

//
// temporary dats
//
template <class T>
op_dat op_decl_dat_temp(op_set set, int dim, char const *type, T *data,
                        char const *name)
{
  return op_decl_dat_temp_char(set, dim, type, sizeof(T), name);
}

inline int op_free_dat_temp(op_dat dat) { return op_free_dat_temp_char(dat); }

//
// fetch data
//
template <class T>
void op_fetch_data(op_dat dat, T *usr_ptr)
{
  op_fetch_data_char(dat, (char *)usr_ptr);
}
template <class T>
void op_fetch_data(op_dat dat, T *usr_ptr, char *temp_pinned)
{
  op_fetch_data_char2(dat, (char *)usr_ptr, temp_pinned);
}

template <class T>
void op_fetch_data_idx(op_dat dat, T *usr_ptr, int low, int high)
{
  op_fetch_data_idx_char(dat, (char *)usr_ptr, low, high);
}

//
// wrapper functions to handle MPI global reductions
//

inline void op_mpi_reduce(op_arg *args, float *data)
{
  op_mpi_reduce_float(args, data);
}

inline void op_mpi_reduce(op_arg *args, double *data)
{
  op_mpi_reduce_double(args, data);
}

inline void op_mpi_reduce(op_arg *args, int *data)
{
  op_mpi_reduce_int(args, data);
}

// needed as a dummy, "do nothing" routine for the non-mpi backends
template <class T>
void op_mpi_reduce(op_arg *args, T *data) {}

#endif /* __OP_LIB_CPP_H */
