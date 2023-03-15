#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>

#include "def2d.h"
extern para param;

#include <iostream>
#include <time.h>

#include "ic.h"
#include "rhs.h"
#include "bc2d.h"
#include "misc.h"
#include "flux2d.h"
#include "poly1d.h"
#include "poly2d.h"
#include "quadrature.h"
#include "m_math_def.h"

#include "op_seq.h"

#include "rhscal.h"
#include "rhscal_auxvar.h"
#include "rhscal_auxsideflux.h"
#include "rhscal_auxsideflux2.h"
#include "rhscal_calsideflux.h"
#include "rhscal_calsideflux2.h"

void CalRHS(op_set eles, op_set sids, op_set isds, op_set bsds,
	    op_map peside_s, op_map pisid, op_map pnisd_eL, op_map pbsid, op_map pbsesid, op_map pnbsd_eL,
	    op_dat p_pvar, op_dat p_pele, op_dat p_psid, op_dat p_a, op_dat p_rhs)

{
	if (param.ivis == 1)
	{
		op_par_loop(rhscal_auxsideflux, "rhscal_auxsideflux", isds,
			    op_arg_dat(p_pvar, 0, pnisd_eL, 1, "uvar", OP_READ),
			    op_arg_dat(p_pvar, 1, pnisd_eL, 1, "uvar", OP_READ),
			    op_arg_dat(p_pele, 0, pnisd_eL, 1, "elem", OP_READ),
			    op_arg_dat(p_pele, 1, pnisd_eL, 1, "elem", OP_READ),
			    op_arg_dat(p_psid, 0, pisid, 1, "side", OP_RW));

		op_par_loop(rhscal_auxsideflux2, "rhscal_auxsideflux2", bsds,
			    op_arg_dat(p_psid, 0, pbsid, 1, "side", OP_RW),
			    op_arg_dat(p_psid, 0, pbsesid, 1, "side", OP_READ),
			    op_arg_dat(p_pvar, 0, pnbsd_eL, 1, "uvar", OP_READ),
			    op_arg_dat(p_pvar, 1, pnbsd_eL, 1, "uvar", OP_READ),
			    op_arg_dat(p_pele, 0, pnbsd_eL, 1, "elem", OP_READ),
			    op_arg_dat(p_pele, 1, pnbsd_eL, 1, "elem", OP_READ));

		op_par_loop(rhscal_auxvar, "rhscal_auxvar", eles,
			    op_arg_dat(p_psid, -4, peside_s, 1, "side", OP_READ),
			    op_arg_dat(p_a, -1, OP_ID, 1, "int", OP_READ),
			    op_arg_dat(p_pele, -1, OP_ID, 1, "elem", OP_READ),
			    op_arg_dat(p_pvar, -1, OP_ID, 1, "uvar", OP_RW));
	}

	op_par_loop(rhscal_calsideflux, "rhscal_calsideflux", isds,
		    op_arg_dat(p_pvar, 0, pnisd_eL, 1, "uvar", OP_READ),
		    op_arg_dat(p_pvar, 1, pnisd_eL, 1, "uvar", OP_READ),
		    op_arg_dat(p_pele, 0, pnisd_eL, 1, "elem", OP_READ),
		    op_arg_dat(p_pele, 1, pnisd_eL, 1, "elem", OP_READ),
		    op_arg_dat(p_psid, 0, pisid, 1, "side", OP_RW));

	op_par_loop(rhscal_calsideflux2, "rhscal_calsideflux2", bsds,
		    op_arg_dat(p_pvar, 0, pnbsd_eL, 1, "uvar", OP_READ),
		    op_arg_dat(p_pvar, 1, pnbsd_eL, 1, "uvar", OP_READ),
		    op_arg_dat(p_pele, 0, pnbsd_eL, 1, "elem", OP_READ),
		    op_arg_dat(p_pele, 1, pnbsd_eL, 1, "elem", OP_READ),
		    op_arg_dat(p_psid, 0, pbsesid, 1, "side", OP_READ),
		    op_arg_dat(p_psid, 0, pbsid, 1, "side", OP_RW));

	op_par_loop(rhscal, "rhscal", eles,
		    op_arg_dat(p_psid, -4, peside_s, 1, "side", OP_READ),
		    op_arg_dat(p_a, -1, OP_ID, 1, "int", OP_READ),
		    op_arg_dat(p_pele, -1, OP_ID, 1, "elem", OP_READ),
		    op_arg_dat(p_rhs, -1, OP_ID, 40, "double", OP_WRITE),
		    op_arg_dat(p_pvar, -1, OP_ID, 1, "uvar", OP_RW));
}
