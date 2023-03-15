
// hour2d.cpp : Defines the entry point for the console application.
// rewrite liming 2011 11 14

// #include <cuda.h>
// #include <cuda_runtime.h>
// #include <cuda_runtime_api.h>
#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>

#include "def2d.h"
para param;

#include "bc2d.h"
#include "ic.h"
#include "rhs.h"
#include "misc.h"
#include "solio.h"
#include "recons.h"
#include "gridio.h"
#include "poly1d.h"
#include "poly2d.h"
#include "postproc.h"
#include "solparam.h"
#include "gridparam.h"
#include "rcm.h"
#include "m_math_def.h"
#include "omp.h"

#include "op_seq.h"
#include "hour_ele.h"
#include "hour_ele2.h"
#include "hour_caldt.h"
#include "recons_1wh.h"
#include "recons_0wh.h"
#include "user/op_save_map.h"

int main(int argc, char *argv[])
{
// omp_set_dynamic( 0 );
// omp_set_num_threads( MAX_THREADS );
#pragma omp parallel
	{
#pragma omp single
		{
			printf("Number of CPU cores: %d\n", omp_get_num_procs());
			printf("the threadID: %d\n", omp_get_thread_num());
			printf("Number of Threads: %d\n", omp_get_num_threads());
		}
	}
	int nele = 0;
	int nsid = 0;
	int nnod = 0;
	int nbosid = 0;
	int ninsid = 0;
	int ncvbd = 0;

	elem *pele = NULL;
	side *psid = NULL;
	node *pnod = NULL;
	uvar *pvar = NULL;
	int *pbosid = NULL;
	int *pinsid = NULL;
	cvbd2d *pcvbd = NULL;

	// read the parameters
	ReadParam(&param, "../param2d.dat");
	int rcm = 2;
	if (argc > 1)
	{
		int iteNum = atoi(argv[1]);
		param.itmax = iteNum;
	}
	if (argc > 2)
	{
		char *stemp;
		stemp = argv[2];
		strcpy(param.gfile, stemp);
	}
	if (argc > 3)
	{
		rcm = atoi(argv[3]);
	}

	int *pbsd_inv;
	pbsd_inv = ReadGridTetrex(param.gfile, nele, nsid, nnod, pele, psid, pnod, pvar,
				  nbosid, ninsid, ncvbd, pbosid, pinsid, pcvbd, rcm);

	int colorNum;
	int **coloredSide = NULL;
	// coloredSide = colorEdge( nele, nsid, colorNum, pele, psid );

	CheckGrid(nele, nsid, nnod, pele, psid, pnod);

	CalBasArr1d(&param);
	CalBasArrTri();
	CalBasArrRect();

	CalSidePara(nele, nsid, nnod, pele, psid, pnod, pvar, &param);

	// calculate element parameter
	CalElemPara(nele, nsid, nnod, pele, psid, pnod, pvar, &param);

	printf("\nafter CalElemPara\n\n");

	int it, it0;
	double t, t0, cputm0;
	double res10[NEQ], res20[NEQ], resx0[NEQ];

	// initial condition
	Init(
	    0, it0, t0, cputm0, nele, nsid, nnod, pele, psid, pnod, pvar, &param,
	    res10, res20, resx0);
	it = it0;
	t = t0;

	int nModi = PosiCorreProc(
	    0, it, t, cputm0, nele, nsid, nnod, pele, psid, pnod, pvar, &param);

	if (nModi < 0)
	{
		printf("Init p0 negative: %d\n", nModi);
		getchar();
		return 1;
	}
	if (nModi > 0)
	{
		printf("Init pk negative: %d\n", nModi);
		getchar();
	} // return 1; }

	// saving the initial solution
	printf("\nsaving the initial solution?\n");
	// getchar();

	// the residual file
	FILE *sTxtRes = NULL;
	FILE *sDatRes = NULL;
	// Ö»ÔÚ¶¨³£Çé¿ö
	if (param.isteady != UNSTEADY)
	{
		const char *filemode;
		if (param.iconti == 0)
			filemode = "w";
		else
			filemode = "a+";
		sTxtRes = fopen("res.txt", filemode); // linux file open mode
		if (sTxtRes == NULL)
		{
			printf("The file 'res.txt' was not opened\n");
			exit(1);
		}
		sDatRes = fopen("res.dat", filemode);
		if (sDatRes == NULL)
		{
			printf("The file 'res.dat' was not opened\n");
			exit(1);
		}
		if (param.iconti == 0)
		{
			fprintf(sDatRes, "VARIABLES = it, t, cpu,");
			for (int j = 0; j < NEQ; j++)
				fprintf(sDatRes, " L<sub><math>5</math></sub>%s L<sub>2</sub>%s L<sub>1</sub>%s",
					VAR_CV_STR[j], VAR_CV_STR[j], VAR_CV_STR[j]);
			fprintf(sDatRes, "\n");
		}
	}
	if (param.iconti == 0)
	{
		FILE *stxt = fopen("preci.dat", "a");
		if (stxt != NULL)
		{
			fclose(stxt);
		}
	}

	DisDetect(
	    0, it, t, 0, nele, nsid, nnod, pele, psid, pnod, pvar, &param);

	Recons(
	    0, it, t, 0, nele, nsid, nnod, pele, psid, pnod, pvar, &param,
	    nbosid, ninsid, pbosid, pinsid, pcvbd);

	// basPreciAna( 0, it, t, cputm0, nele, nsid, nnod, pele, psid, pnod, pvar, &param );
	// getchar();
	// PreciAna(
	//	0, it, t, cputm0, nele, nsid, nnod, pele, psid, pnod, pvar, &param );
	// getchar();
	OutSol(
	    0, it, t, cputm0, nele, nsid, nnod, pele, psid, pnod, pvar, &param);

	// return 0;
	// getchar();
	double reslimi = pow(10.0, double(param.resi));
	/// int iRuKu = param.iRuKu;
	double tCFL = param.CFL;
	int glbdgdof = DOFofOrd2d[param.pDG];

	if ((param.isteady != UNSTEADY) && (tCFL > 0.0))
		if ((it0 + 1) <= param.CFLIncStp)
			param.CFL = (it0 + 1) * tCFL / param.CFLIncStp;

	double dt = CalDt(0, t, nele, pele, pvar, &param);
	// ·Ç¶¨³£ÎÊÌâ
	if (param.isteady == UNSTEADY)
		dt = m_min(dt, param.tmax - t);

	if ((param.isteady != STEADY) || (param.CFL < 0.0)) // È«³¡Ê¹ÓÃÍ³Ò»µÄÊ±¼ä²½³¤
		for (int e = 0; e < nele; e++)
			pvar[e].dt = dt;

	// start calculate
	printf("start at \nit= % 8d, t= % 14.7f, cpu= % 10.3f\n", it, t, cputm0);
	// getchar();
	clock_t cputmst = (clock_t)((double)clock() - cputm0 * CLOCKS_PER_SEC);

	double CFL = param.CFL;
	double *rhs = (double *)malloc(10 * 4 * nele * sizeof(double));
	memset(rhs, 0, 10 * 4 * nele * sizeof(double));

	int *aele = (int *)malloc(nele * sizeof(int));
	for (int times = 0; times < nele; times++)
		aele[times] = times;
	// Ó³Éä±í
	// ÄÚ±ß
	int *nisd_eL = (int *)malloc(2 * ninsid * sizeof(int));
	memset(nisd_eL, 0, 2 * ninsid * sizeof(int));
#pragma omp parallel for
	for (int idx = 0; idx < (int)ninsid; ++idx)
	{
		unsigned int s = pinsid[idx];
		nisd_eL[2 * idx] = psid[s].ofElem[1];
		nisd_eL[2 * idx + 1] = psid[s].ofElem[0];
	}
	// ±ß½ç
	int *nbsd_eL = (int *)malloc(2 * nbosid * sizeof(int));
	memset(nbsd_eL, 0, 2 * nbosid * sizeof(int));
	int *bsesid = (int *)malloc(nbosid * sizeof(int));
	memset(bsesid, 0, nbosid * sizeof(int));
#pragma omp parallel for
	for (int idx = 0; idx < (int)nbosid; ++idx)
	{
		unsigned int s = pbosid[idx];
		int bceR = psid[s].ofElem[0];
		nbsd_eL[2 * idx] = psid[s].ofElem[1];
		nbsd_eL[2 * idx + 1] = GetElement(bceR);

		unsigned int eR = GetElement(bceR);
		int tsR = psid[s].lsdidx[0];
		bsesid[idx] = pele[eR].Side[tsR];
	}
	// µ¥Ôª
	int *eside_s = (int *)malloc(nType1 * nele * sizeof(int));
	memset(eside_s, 0, nType1 * nele * sizeof(int));
#pragma omp parallel for
	for (int e = 0; e < (int)nele; e++)
	{
		for (unsigned int ts = 0; ts < nType1; ts++)
			eside_s[e * nType1 + ts] = pele[e].Side[ts];
	}

	op_init(argc, argv, 2);

	double cpu_t1, cpu_t2, wall_t1, wall_t2;

	op_set eles = op_decl_set(nele, "eles");
	op_set sids = op_decl_set(nsid, "sids");
	op_set isds = op_decl_set(ninsid, "isds");
	op_set bsds = op_decl_set(nbosid, "bsds");

	op_map pisid = op_decl_map(isds, sids, 1, pinsid, "pisid");
	op_map pnisd_eL = op_decl_map(isds, eles, 2, nisd_eL, "pnisd_eL");

	op_map pbsid = op_decl_map(bsds, sids, 1, pbosid, "pbsid");
	op_map pbsesid = op_decl_map(bsds, sids, 1, bsesid, "pbsesid");
	op_map pnbsd_eL = op_decl_map(bsds, eles, 2, nbsd_eL, "pnbsd_eL");

	op_map peside_s = op_decl_map(eles, sids, nType1, eside_s, "peside_s");

	// zongjing
	saveMap(peside_s);
	saveMap(pbsid);
	// zongjing

	op_dat p_a = op_decl_dat(eles, 1, "int", aele, "p_a");
	op_dat p_pele = op_decl_dat(eles, 1, "elem", pele, "p_pele");
	op_dat p_psid = op_decl_dat(sids, 1, "side", psid, "p_psid");
	op_dat p_pvar = op_decl_dat(eles, 1, "uvar", pvar, "p_pvar");
	op_dat p_rhs = op_decl_dat(eles, 40, "double", rhs, "p_rhs");
	double *pvar_grad = (double *)malloc(nele * 80 * sizeof(double));
	op_dat p_pvar_grad = op_decl_dat(eles, 2 * 10 * 4, "double", pvar_grad, "p_pvar_grad");

	op_decl_const(16, "int", DOFofOrd2d);
	op_decl_const(7, "int", SymQdTriDim);
	op_decl_const(10, "int", QdNOrd2NPt1d);
	op_decl_const(6 * 5 * 3, "double", GauQd1dCoe);
	op_decl_const((MAXSEGQDPTS + 1) * MAXSEGQDPTS * 2, "double", QdPtCoe1d);

	op_decl_const(1, "para", &param);

	op_decl_const(7 * 12 * 4, "double", SymQdTriCoe);
	op_decl_const(7 * 12 * 10, "double", QdPtBasTri);

	op_decl_const(6 * 25 * 4, "double", QdPtCoeRect);
	op_decl_const(6 * 25 * 10, "double", QdPtBasRect);

	op_decl_const(6 * 6 * 5 * 10, "double", SdQdPtBasTri);
	op_decl_const(8 * 6 * 5 * 10, "double", SdQdPtBasRect);

	op_decl_const(7 * 12 * 10, "double", QdPtBasDxiTri, 1);
	op_decl_const(7 * 12 * 10, "double", QdPtBasDetTri, 1);
	op_decl_const(6 * 25 * 10, "double", QdPtBasDxiRect, 1);
	op_decl_const(6 * 25 * 10, "double", QdPtBasDetRect, 1);
	op_diagnostic_output();
	op_timers(&cpu_t1, &wall_t1);

	double time_begin, time_end, pinned_begin, pinned_end;
	time_begin = omp_get_wtime();
	double kernel_time = 0, fetch_dat_time = 0;
	double *temp_pinned;
	// cudaHostAlloc((void **)&temp_pinned, p_rhs->size * p_rhs->set->size * sizeof(char), cudaHostAllocDefault);

	for (it = it0 + 1;; it++)
	{
		pinned_begin = omp_get_wtime();
		if ((it <= param.CFLIncStp) || ((it % param.dtfre) == 0))
		{
			op_par_loop(hour_caldt, "hour_caldt", eles,
				    op_arg_dat(p_pele, -1, OP_ID, 1, "elem", OP_READ),
				    op_arg_dat(p_pvar, -1, OP_ID, 1, "uvar", OP_RW),
				    op_arg_gbl(&CFL, 1, "double", OP_READ),
				    op_arg_gbl(&dt, 1, "double", OP_MIN));
		}

		CalRHS(eles, sids, isds, bsds, peside_s, pisid, pnisd_eL, pbsid,
		       pbsesid, pnbsd_eL, p_pvar, p_pele, p_psid, p_a, p_rhs, p_pvar_grad);

		op_par_loop(hour_ele, "hour_ele", eles,
			    op_arg_dat(p_rhs, -1, OP_ID, 40, "double", OP_READ),
			    op_arg_dat(p_pvar, -1, OP_ID, 1, "uvar", OP_RW));

		op_par_loop(recons_1wh, "recons_1wh", eles,
			    op_arg_dat(p_pvar, -1, OP_ID, 1, "uvar", OP_RW));

		CalRHS(eles, sids, isds, bsds, peside_s, pisid, pnisd_eL, pbsid,
		       pbsesid, pnbsd_eL, p_pvar, p_pele, p_psid, p_a, p_rhs, p_pvar_grad);

		op_par_loop(hour_ele2, "hour_ele2", eles,
			    op_arg_dat(p_rhs, -1, OP_ID, 40, "double", OP_READ),
			    op_arg_dat(p_pvar, -1, OP_ID, 1, "uvar", OP_INC));

		op_par_loop(recons_0wh, "recons_0wh", eles,
			    op_arg_dat(p_pvar, -1, OP_ID, 1, "uvar", OP_RW));
		pinned_end = omp_get_wtime();
		kernel_time += pinned_end - pinned_begin;
		pinned_begin = omp_get_wtime();
		// op_fetch_data(p_pvar, pvar);
		op_fetch_data(p_rhs, rhs);
		// op_fetch_data(p_rhs, rhs, (char *)temp_pinned);
		pinned_end = omp_get_wtime();
		fetch_dat_time += pinned_end - pinned_begin;

		//////////////////////////////////////////////////////////////////////////CalRes:
		t += dt;

		// RKÊ±¼äÍÆ½øºó²Ù×÷res£¬ÒªÇó´ËÊ±rhs»¹Ã»±»ÆÆ»µ(×îºóÒ»²½µÄ)
		int eresx[MAXDGDOF][NEQ];
		double res1[MAXDGDOF][NEQ], res2[MAXDGDOF][NEQ], resx[MAXDGDOF][NEQ];
		double L2rho;

		// Ö»ÔÚ¶¨³£Çé¿ö
		if (param.isteady != UNSTEADY)
		{
			// RKÊ±¼äÍÆ½øºó²Ù×÷res£¬ÒªÇó´ËÊ±rhs»¹Ã»±»ÆÆ»µ
			memset(eresx, 0, sizeof(int) * (NEQ * MAXDGDOF));
			memset(res1, 0, sizeof(double) * (NEQ * MAXDGDOF));
			memset(res2, 0, sizeof(double) * (NEQ * MAXDGDOF));
			memset(resx, 0, sizeof(double) * (NEQ * MAXDGDOF));
			for (int e = 0; e < nele; e++)
			{
				int ndgord = pvar[e].pdg;
				int ndgdof = DOFofOrd2d[ndgord];
				for (int i = 0; i < ndgdof; i++)
				{
					for (int j = 0; j < NEQ; j++)
					{
						//						double tmp = m_abs( pvar[e].rhs[i][j] );
						double tmp = m_abs(rhs[e * 40 + i * 4 + j]);
						// double tmp = m_abs(temp_pinned[e + i * 4 * nele + j * nele]);

						res1[i][j] += tmp;
						res2[i][j] += tmp * tmp;

						if (tmp > resx[i][j])
						{
							resx[i][j] = tmp;
							eresx[i][j] = e;
						}
					}
				}
			}
			for (int i = 0; i < glbdgdof; i++)
				for (int j = 0; j < NEQ; j++)
				{
					res1[i][j] /= nele;
					res2[i][j] = sqrt(res2[i][j] / nele);
				}
			// Ç°¶þÊ®²½Æ½¾ùÖµ it = 1--20
			if (it < 21)
				for (int j = 0; j < NEQ; j++)
				{
					res10[j] = (res10[j] * (it - 1) + res1[0][j]) / it;
					res20[j] = (res20[j] * (it - 1) + res2[0][j]) / it;
					resx0[j] = (resx0[j] * (it - 1) + resx[0][j]) / it;
				}

			L2rho = res2[0][0] / res20[0];

			if ((tCFL > 0.0) && ((it + 1) <= param.CFLIncStp))
				CFL = (it + 1) * tCFL / param.CFLIncStp;

			// printf( "main res end\n" );
		}

		// ²»ÊÇÃ¿²½¶¼Êä³ö¼ÆËãÐÅÏ¢
		if ((it < 241) || ((it % param.wsfre) == 0))
		{
			double elaps = (double)(clock() - cputmst) / CLOCKS_PER_SEC;
			for (int i = 0; i < glbdgdof; i++)
			{
				for (int j = 0; j < NEQ; j++)
				{
					// printf( "%d %d %d %.18e %.18e %.18e\n",it, i, j, resx[i][j], res1[i][j], res2[i][j] );
				}
			}
			// ÔÚ¶¨³£Çé¿ö
			if (param.isteady != UNSTEADY)
			{
				fprintf(sDatRes, "% 8d % 12.5f % 10.3f", it, t, elaps); // res.dat
				for (int j = 0; j < NEQ; j++)
					fprintf(sDatRes, " % 9.3e % 9.3e % 9.3e", resx[0][j], res2[0][j], res1[0][j]);
				fprintf(sDatRes, "\n");

				// res.txt
				fprintf(sTxtRes, "it=% 8d, t=% 10.5f, cpu=% 10.3f, L2=% 9.3e*% 9.3e, Linf=% 9.3e, "
						 "xc=% 9.3e, yc=% 9.3e",
					it, t, elaps, L2rho, res20[0], resx[0][0],
					pele[eresx[0][0]].xyc[0], pele[eresx[0][0]].xyc[1]);
				for (int i = 0; i < glbdgdof; i++)
					for (int j = 0; j < NEQ; j++)
						fprintf(sTxtRes, " % 9.3e % 9.3e % 9.3e", resx[i][j], res2[i][j], res1[i][j]);
				fprintf(sTxtRes, "\n");

				printf("it= % 8d, t= % 14.7f, cpu= % 10.3f, L2rhoDrop= % 9.3e\n",
				       it, t, elaps, L2rho);
			}
			else
				printf("it= % 8d, t= % 14.7f, cpu= % 10.3f\n", it, t, elaps);
		}

		if ((it % param.wffre) == 0)
		{
			double elaps = (double)(clock() - cputmst) / CLOCKS_PER_SEC; // elapsed_time

			printf("save data at \nit= % 8d, t= % 14.7f, cpu= % 10.3f\n", it, t, elaps);

			//			PreciAna( 0, it, t, elaps, nele, nsid, nnod, pele, psid, pnod, pvar, &param );
			//			OutSol(0, it, t, elaps, nele, nsid, nnod, pele, psid, pnod, pvar, &param );
			//			SaveData( 0, it, t, elaps, res10, res20, resx0, nele, pvar, &param );

			// getchar();
		}

		// testing for the end of the iteration loop
		if (it == param.itmax)
		{
			printf("main End of execution: maximal number of time steps : %d\n", it);
			break;
		}

		// Ö»ÔÚ¶¨³£Çé¿ö testing for the end of the res drop
		if (param.isteady != UNSTEADY)
		{
			if (L2rho < reslimi)
			{
				printf("main End of execution: L2rho drop : % 9.3e\n", L2rho);
				break;
			}
		}

		// Ö»ÔÚ·Ç¶¨³£Çé¿ö testing for the end of the physical time
		if (param.isteady == UNSTEADY)
		{
			// if( fabs( t - param.dtmax ) < 0.000005 )
			if (t > (param.tmax - 0.000005))
			{
				printf("main End of execution : maximal physical time: % 14.7f\n", t);
				break;
			}
		}
	} // end of it
	// clock_t print_st = clock();
	// double print_sec = (double)(print_st - ite_st) / CLOCKS_PER_SEC;

	// cudaFreeHost(temp_pinned);

	time_end = omp_get_wtime();
	printf("the cfd solver uses: %14.7f\n", time_end - time_begin);

	op_timers(&cpu_t2, &wall_t2);
	op_printf("Max total runtime = %f\n", wall_t2 - wall_t1);
	op_printf("kernel_time = %f\n", kernel_time);
	op_printf("fetch_dat_time = %f\n", fetch_dat_time);
	op_printf("cpu runtime = %f\n", wall_t2 - wall_t1 - kernel_time - fetch_dat_time);

	op_fetch_data(p_pvar, pvar);
	op_fetch_data(p_psid, psid);
	op_timing_output();
	op_printf("op2 exit \n\n");
	op_exit();

	int *nbosidtemp = (int *)calloc(nbosid, sizeof(int));
	// EndofMain it:
	if (true) //(it%param.wffre)!=0 ) // µÚÒ»²½Ê±Îª1 Êä³ö¼ÆËã½á¹û
	{
		if (rcm == 2)
		{
			for (int i = 0; i < nbosid; i++)
			{
				nbosidtemp[i] = pbsd_inv[i];
				// printf( "%d  %d\n", i+1, pbsd_inv[i]+1 );
			}
		}
		else
		{
			for (int i = 0; i < nbosid; i++)
			{
				nbosidtemp[i] = pbosid[i];
			}
		}
		double elaps = (double)(clock() - cputmst) / CLOCKS_PER_SEC; // elapsed_time
		printf("save data at \nit= % 8d, t= % 14.7f, cpu= % 10.3f\n", it, t, elaps);
		PreciAna(0, it, t, elaps, nele, nsid, nnod, pele, psid, pnod, pvar, &param);
		OutSol(0, it, t, elaps, nele, nsid, nnod, pele, psid, pnod, pvar, &param);
		OutWallSol(0, it, t, elaps, nbosid, nbosidtemp, pele, psid, pnod, pvar, &param);
		SaveData(0, it, t, elaps, res10, res20, resx0, nele, pvar, &param);
	}

	if (pele != NULL)
		free(pele);
	if (psid != NULL)
		free(psid);
	if (pnod != NULL)
		free(pnod);
	if (pvar != NULL)
		free(pvar);
	if (aele != NULL)
		free(aele);
	if (nisd_eL != NULL)
		free(nisd_eL);
	if (nbsd_eL != NULL)
		free(nbsd_eL);
	if (bsesid != NULL)
		free(bsesid);
	if (eside_s != NULL)
		free(eside_s);
	free(pbsd_inv);
	free(nbosidtemp);

	/* All  files are closed: */
	// int numclosed;
	// numclosed = fcloseall( );
	// printf( "\n%u files saved and closed\n", numclosed );
	// getchar();

	return 0;
}
