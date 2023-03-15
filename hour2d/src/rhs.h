#ifndef _RHS_H_
#define _RHS_H_

#include "def2d.h"
#include "op_seq.h"
extern para param;
void CalRHS(
    op_set eles, op_set sids, op_set isds, op_set bsds,
    op_map peside_s, op_map pisid, op_map pnisd_eL, op_map pbsid, op_map pbsesid, op_map pnbsd_eL,
    op_dat p_pvar, op_dat p_pele, op_dat p_psid, op_dat p_a, op_dat p_rhs);

#endif