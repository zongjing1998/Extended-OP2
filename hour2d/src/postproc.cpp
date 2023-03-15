// post process
// 2011 07 28 liming

#include <math.h>
#include <stdio.h>
//#include <conio.h>
#include <string.h>
#include <stdlib.h>

#include "ic.h"
#include "bc2d.h"
#include "misc.h"
#include "def2d.h"
#include "poly1d.h"
#include "poly2d.h"
#include "limiter.h"
#include "postproc.h"
#include "quadrature.h"

#include "m_math_def.h"

#include "omp.h"
// for precision analysis
void PreciAna(
    unsigned int irk, int it, double t, double cputm,
    unsigned int nele, unsigned int nsid, unsigned int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa)
{
	FILE *stxt = NULL;
	// errno_t err;

	/*// append to the file
	if( ( err = fopen_s( &stxt, "preci.dat", "a" ) ) !=0 )
	{
		PRINTERRORN( "The file 'preci.dat' was not opened\n" );
	}*/
	stxt = fopen("preci.dat", "a");
	if (stxt == NULL)
	{
		PRINTERRORN("The file 'preci.dat' was not opened\n");
	}
#if EQUATIO == ADVECTI
	// fprintf( stxt, " advective eq\n" );
#elif EQUATIO == BURGERS
	// fprintf( stxt, " Burgers eq\n" );
#elif EQUATIO == EULNSEQ
	// fprintf( stxt, " Euler eqs\n" );
#endif

	fprintf(stxt, "% 8d % 12.5f % 10.3f", it, t, cputm);

	// double dt = m_abs( ppa->CFL );	//ע��˴�

	double avgL1 = 0.0;
	double avgL2 = 0.0;
	double ptLinf = 0.0;
	double tolvol = 0.0;

	for (unsigned int e = 0; e < nele; e++)
	{
		// tolu += pvar[e].wh[0][0]/pele[e].volr;
		tolvol += 1.0 / pele[e].volr;

		int *inode = pele[e].Node;
		unsigned int nType = ELEMTYPE;
#ifdef QUAD
		if (inode[3] == inode[0])
			nType = 3;
#endif
		double ndx[ELEMTYPE], ndy[ELEMTYPE];
		for (unsigned int i = 0; i < ELEMTYPE; i++)
		{
			ndx[i] = pnod[inode[i]].xy[0];
			ndy[i] = pnod[inode[i]].xy[1];
		}

		unsigned int pfv = pvar[e].pfv;
		unsigned int nfvdof = DOFofOrd2d[pfv];

		unsigned int nqdord = pfv * 2;
		if (nqdord == 0)
			nqdord = 1;

		double tmpL2 = 0.0;
		double tmpL1 = 0.0;

		unsigned int nqdpt;
		double(*EmQdCoe)[4], (*EmQdBas)[MAXFVDOF];
		if (nType == 3) // tri
		{
			nqdpt = SymQdTriDim[nqdord];
			EmQdCoe = (double(*)[4])SymQdTriCoe[nqdord];
			EmQdBas = (double(*)[MAXFVDOF])QdPtBasTri[nqdord];
		}
		else // rect
		{
			nqdpt = MAXSEGQDPTS; // QdNOrd2NPt1d[nqdord];
			// nqdpt = nqdord/2+1;
			EmQdCoe = (double(*)[4])QdPtCoeRect[nqdpt];
			EmQdBas = (double(*)[MAXFVDOF])QdPtBasRect[nqdpt];
			nqdpt *= nqdpt;
		}

		for (unsigned int qd = 0; qd < nqdpt; qd++)
		{
			double xi = EmQdCoe[qd][0];
			double et = EmQdCoe[qd][1];
			double coe = EmQdCoe[qd][3];

			double xx, yy;
			if (nType == 3)
				Loc2GlbTri(xx, yy, xi, et, ndx, ndy);
			else
				Loc2GlbRect(xx, yy, xi, et, ndx, ndy);

			double tcv[NEQ];
			GetVar(nfvdof, pvar[e].wh, EmQdBas[qd], tcv);

			double tmcv[NEQ];
			double diff;
#if EQUATIO == ADVECTI
			tmcv[0] = uIC(xx, yy, t, ppa);
			diff = m_abs(tcv[0] - tmcv[0]);
#elif EQUATIO == BURGERS
			tmcv[0] = Burg2d(xx, yy, t, nele, nsid, nnod, pele, psid, pnod, pvar, ppa);
			diff = m_abs(tcv[0] - tmcv[0]);
#elif EQUATIO == EULNSEQ

			//			double tmcvT = TIC( xx, yy, t, ppa );
			//			double pre;
			//			CalPre( tcv, tcv, pre );
			//			double actT = pre / tcv[0] * gamma * ppa->Ma * ppa->Ma;
			//			diff = m_abs( actT- tmcvT );

			// ������
			tmcv[0] = RIC(xx, yy, t, ppa);
			diff = m_abs(tcv[0] - tmcv[0]);
			// GetVar( nfvdof, pvar[e].divCV[1], QdPtBasTri[nqdord][np], tcv );
			// tmcv[0] = rhouicdy( xx, yy, t, ppa );
			// tmcv[0] = -0.2*PI2*sin(PI*xx)*sin(PI2*yy);
			// double pre;
			// CalPre2d( tcv, tcv, pre );
			// double T = pre/(ppa->R*tcv[0]);
			//??tmcv[0] = rhouic( xx, yy, t, ppa );
			//??diff = m_abs( tcv[1]/tcv[0] - tmcv[0] );
			/*double pre;
			const double* pcv = tcv;
			int ret = CalPre2d( tcv, pvar[e].wh[0], pre );
			if( ret==1 ) pcv = pvar[e].wh[0];
			diff = m_abs( pre * pow(pcv[0], -gamma )*ppa->Sinfn1 - 1.0 );*/
#endif
			ptLinf = m_max(ptLinf, diff);
			tmpL2 += diff * diff * coe;
			tmpL1 += diff * coe;
		} // element gauss quad

		avgL2 += tmpL2 / pele[e].volr;
		avgL1 += tmpL1 / pele[e].volr;
	} // end of elem

	avgL1 = avgL1 / tolvol;
	avgL2 = sqrt(avgL2 / tolvol);

	fprintf(stxt, " % 9.3e", ptLinf);
	fprintf(stxt, " % 9.3e", avgL2);
	fprintf(stxt, " % 9.3e\n", avgL1);

	fclose(stxt);
}

void CalMtxLR2d(
    double uh[NEQ], double *norm, double (&MtxR)[NEQ][NEQ], double (&MtxL)[NEQ][NEQ])
{
	double nx = norm[0];
	double ny = norm[1];
	double rr = uh[0];
	double ru = uh[1];
	double rv = uh[2];
	double ee = uh[3];

	double uu = ru / rr;
	double vv = rv / rr;
	double q2 = uu * uu + vv * vv;
	double f2 = 0.5 * gamm1 * q2; //
	double pp = gamm1 * (ee - 0.5 * rr * q2);
	double cc = sqrt(gammaa * pp / rr);
	double c2r = 1.0 / cc / cc;
	double c2rg = gamm1 * c2r;
	double c2rg1d2 = 0.5 * gamm1 * c2r;
	double hh = (ee + pp) / rr;

	double un = nx * uu + ny * vv;
	double nxcc = nx * cc;
	double nycc = ny * cc;
	double uncc = un * cc;

	MtxR[0][0] = 1.0;
	MtxR[0][1] = 0.0;
	MtxR[0][2] = 1.0;
	MtxR[0][3] = 1.0;
	MtxR[1][0] = uu;
	MtxR[1][1] = -ny;
	MtxR[1][2] = uu + nxcc;
	MtxR[1][3] = uu - nxcc;
	MtxR[2][0] = vv;
	MtxR[2][1] = nx;
	MtxR[2][2] = vv + nycc;
	MtxR[2][3] = vv - nycc;
	MtxR[3][0] = 0.5 * q2;
	MtxR[3][1] = nx * vv - ny * uu;
	MtxR[3][2] = hh + uncc;
	MtxR[3][3] = hh - uncc;

	MtxL[0][0] = 1.0 - f2 * c2r;
	MtxL[0][1] = uu * c2rg;
	MtxL[0][2] = vv * c2rg;
	MtxL[0][3] = -c2rg;
	MtxL[1][0] = ny * uu - nx * vv;
	MtxL[1][1] = -ny;
	MtxL[1][2] = nx;
	MtxL[1][3] = 0.0;
	MtxL[2][0] = (f2 - uncc) * gam1r;
	MtxL[2][1] = -(uu - nxcc * gam1r);
	MtxL[2][2] = -(vv - nycc * gam1r);
	MtxL[2][3] = 1.0;
	MtxL[3][0] = (f2 + uncc) * gam1r;
	MtxL[3][1] = -(uu + nxcc * gam1r);
	MtxL[3][2] = -(vv + nycc * gam1r);
	MtxL[3][3] = 1.0;

	for (unsigned int j = 0; j < NEQ; j++)
	{
		MtxL[2][j] *= c2rg1d2;
		MtxL[3][j] *= c2rg1d2;
	}

	// for( unsigned int j=0; j<NEQ; j++ )
	//	printf( "%f", MtxL[j] );

	//_getch();
}

inline void MtxMulVect(double Mtx[NEQ][NEQ], double *Vect, double *ret)
{
	for (unsigned int i = 0; i < NEQ; i++)
	{
		ret[i] = 0.0;
		for (unsigned int j = 0; j < NEQ; j++)
			ret[i] += Mtx[i][j] * Vect[j];
	}
}

// the "M_TVB"
void TVB_shu(
    unsigned int irk, int it, double t, double cputm,
    unsigned int nele, unsigned int nsid, unsigned int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa)
{
#pragma omp parallel for
	for (int e = 0; e < (int)nele; e++)
	{
		unsigned int ndgord = pvar[e].pdg;
		if (ndgord == 0)
			continue;
		unsigned int ndgdof = DOFofOrd2d[ndgord];
		// unsigned int nfvord = pvar[e].pfv;
		// if( nfvord==0 ) continue;
		// unsigned int nfvdof = DOFofOrd2d[nfvord];

		/*for( unsigned int j=0; j<NEQ; j++ )
			pvar[e].isDiscon[j] = 0.0;*/
		double MtxR[NEQ][NEQ], MtxL[NEQ][NEQ];
		double norm[NDIM], dl;
		int enb;

		// start of x-direction
		dl = psid[pele[e].Side[0]].dl;
		double Mdxdx = ppa->M_TVB * dl * dl;

		if (psid[pele[e].Side[1]].ofElem[1] == e)
			for (unsigned int i = 0; i < NDIM; i++)
				norm[i] = psid[pele[e].Side[1]].norm[i];
		else
			for (unsigned int i = 0; i < NDIM; i++)
				norm[i] = -psid[pele[e].Side[1]].norm[i];

		CalMtxLR2d(pvar[e].uh[irk][0], norm, MtxR, MtxL); // norm of Side[1] means quasi x direction

		// for( unsigned int i=0; i<NEQ; i++ )
		//{
		//	for( unsigned int j=0; j<NEQ; j++ )
		//	{
		//		double ret = 0.0;

		//		for( unsigned int k=0; k<NEQ; k++ )
		//		ret += MtxL[i*NEQ+k] * MtxR[k*NEQ+j];

		//	printf( "%f ", ret );
		//	}
		//}
		//_getch();

		// L. ux, ui+1 - ui, ui - ui-1
		// b.c???
		double ux[NEQ], duxp[NEQ], duxm[NEQ];
		for (unsigned int j = 0; j < NEQ; j++)
			ux[j] = pvar[e].uh[irk][1][j];

		enb = psid[pele[e].Side[1]].ofElem[0] + psid[pele[e].Side[1]].ofElem[1] - ((int)e);
		if (enb < 0) // need bc
			for (unsigned int j = 0; j < NEQ; j++)
				duxp[j] = 3.46411 * ux[j];
		else
			for (unsigned int j = 0; j < NEQ; j++)
				duxp[j] = pvar[enb].uh[irk][0][j] - pvar[e].uh[irk][0][j];

		enb = psid[pele[e].Side[3]].ofElem[0] + psid[pele[e].Side[3]].ofElem[1] - ((int)e);
		if (enb < 0) // need bc
			for (unsigned int j = 0; j < NEQ; j++)
				duxm[j] = 3.46411 * ux[j];
		else
			for (unsigned int j = 0; j < NEQ; j++)
				duxm[j] = pvar[e].uh[irk][0][j] - pvar[enb].uh[irk][0][j];

		double wx[NEQ], dwxp[NEQ], dwxm[NEQ];

		MtxMulVect(MtxL, ux, wx);
		MtxMulVect(MtxL, duxp, dwxp);
		MtxMulVect(MtxL, duxm, dwxm);

		unsigned int flag = 0;
		for (unsigned int j = 0; j < NEQ; j++)
		{
			if (m_abs(wx[j]) <= Mdxdx)
				continue;

			double tmp = minmod3(wx[j], 0.57735 * dwxp[j], 0.57735 * dwxm[j]);
			if (tmp != wx[j])
			{
				flag = 1;
				wx[j] = tmp;
			}
		}

		if (flag == 1) // means something changed!!!
		{
			MtxMulVect(MtxR, wx, ux);

			for (unsigned int j = 0; j < NEQ; j++)
			{
				pvar[e].uh[irk][1][j] = ux[j];
				pvar[e].isDiscon[j] = 1.0;

				for (unsigned int ndof = 3; ndof < ndgdof; ndof++)
					pvar[e].uh[irk][ndof][j] = 0.0;
			}
		}
		// end of x-direction

		// start of y-direction
		dl = psid[pele[e].Side[1]].dl;
		double Mdydy = ppa->M_TVB * dl * dl;

		if (psid[pele[e].Side[2]].ofElem[1] == e)
			for (unsigned int i = 0; i < NDIM; i++)
				norm[i] = psid[pele[e].Side[2]].norm[i];
		else
			for (unsigned int i = 0; i < NDIM; i++)
				norm[i] = -psid[pele[e].Side[2]].norm[i];

		CalMtxLR2d(pvar[e].uh[irk][0], norm, MtxR, MtxL); // norm of Side[2] means quasi y direction

		// L. uy, uj+1 - uj, uj-uj - 1
		// b.c???
		double uy[NEQ], duyp[NEQ], duym[NEQ];
		for (unsigned int j = 0; j < NEQ; j++)
			uy[j] = pvar[e].uh[irk][2][j];

		enb = psid[pele[e].Side[2]].ofElem[0] + psid[pele[e].Side[2]].ofElem[1] - ((int)e);
		if (enb < 0) // need bc
			for (unsigned int j = 0; j < NEQ; j++)
				duyp[j] = 3.46411 * uy[j];
		else
			for (unsigned int j = 0; j < NEQ; j++)
				duyp[j] = pvar[enb].uh[irk][0][j] - pvar[e].uh[irk][0][j];

		enb = psid[pele[e].Side[0]].ofElem[0] + psid[pele[e].Side[0]].ofElem[1] - ((int)e);
		if (enb < 0) // need bc
			for (unsigned int j = 0; j < NEQ; j++)
				duym[j] = 3.46411 * uy[j];
		else
			for (unsigned int j = 0; j < NEQ; j++)
				duym[j] = pvar[e].uh[irk][0][j] - pvar[enb].uh[irk][0][j];

		double wy[NEQ], dwyp[NEQ], dwym[NEQ];

		MtxMulVect(MtxL, uy, wy);
		MtxMulVect(MtxL, duyp, dwyp);
		MtxMulVect(MtxL, duym, dwym);

		flag = 0;
		for (unsigned int j = 0; j < NEQ; j++)
		{
			if (m_abs(wy[j]) <= Mdydy)
				continue;

			double tmp = minmod3(wy[j], 0.57735 * dwyp[j], 0.57735 * dwym[j]);
			if (tmp != wy[j])
			{
				flag = 1;
				wy[j] = tmp;
			}
		}

		if (flag == 1) // means something changed!!!
		{
			MtxMulVect(MtxR, wy, uy);

			for (unsigned int j = 0; j < NEQ; j++)
			{
				pvar[e].uh[irk][2][j] = uy[j];
				pvar[e].isDiscon[j] = 1.0;

				for (unsigned int ndof = 3; ndof < ndgdof; ndof++)
					pvar[e].uh[irk][ndof][j] = 0.0;
			}
		}
		// end of y-direction
	} // end of elem
}

void DisDetect(
    unsigned int irk, int it, double t, double cputm,
    unsigned int nele, unsigned int nsid, unsigned int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa)
{
	return;
	TVB_shu(irk, it, t, cputm, nele, nsid, nnod, pele, psid, pnod, pvar, ppa);

	// MomentP2( irk, it, t, cputm, nele, nsid, nnod, pele, psid, pnod, pvar, ppa );
	return;

#pragma omp parallel for
	for (int e = 0; e < (int)nele; e++)
	{
		for (unsigned int j = 0; j < NEQ; j++)
		{
			double tmp[MAXFVDOF];
			for (unsigned int i = 0; i < MAXFVDOF; i++)
				tmp[i] = 0.0;

			for (unsigned int k = 0; k < MAXFVDOF; k++)
				for (unsigned int i = 0; i < MAXFVDOF; i++)
					tmp[k] += pvar[e].wh[i][j] * osindiMatrix[i][k];

			pvar[e].isDiscon[j] = 0.0;
			for (unsigned int i = 0; i < MAXFVDOF; i++)
				pvar[e].isDiscon[j] += tmp[i] * pvar[e].wh[i][j];
		}
		continue;

		unsigned int ndgord = pvar[e].pdg;
		unsigned int ndgdof = DOFofOrd2d[ndgord];

		for (unsigned int j = 0; j < NEQ; j++)
			pvar[e].isDiscon[j] = 0.0;

		if (ndgord < 1)
			continue;

		unsigned int nfvord = pvar[e].pfv;
		unsigned int nfvdof = DOFofOrd2d[nfvord];

		// unsigned int nnof   = DOFofOrd2d[ndgord];
		// unsigned int nnofm1 = DOFofOrd2d[ndgord-1];
		unsigned int nnof = nfvdof;
		unsigned int nnofm1 = ndgdof;

		double tmp1[NEQ], tmp2[NEQ];
		for (unsigned int j = 0; j < NEQ; j++)
		{
			tmp1[j] = 0.0;
			for (unsigned int i = 0; i < nnofm1; i++)
				tmp1[j] += pvar[e].uh[irk][i][j] * pvar[e].uh[irk][i][j];

			tmp2[j] = 0.0;
			for (unsigned int i = nnofm1; i < nnof; i++)
				tmp2[j] += pvar[e].uh[irk][i][j] * pvar[e].uh[irk][i][j];

			double Se = (tmp2[j] + TOL) / (tmp1[j] + tmp2[j] + TOL + TOL);
			double se = log10(Se);
			double s0 = 1.0 / (ndgord * ndgord * ndgord * ndgord);
			double eps0 = pele[e].dh / ndgord;
			double kap = 4.0;

			pvar[e].isDiscon[j] = se;
		}

		continue;
		//////////////////////////////////////////////////////////////////////////////////////
		// super convergence!!!
		for (unsigned int j = 0; j < NEQ; j++)
			pvar[e].isDiscon[j] = TOL;

		if (ndgord < 1)
			continue;

		if (e == 142)
		{
			int i = 0;
			i++;
		}

		double dhpp1d2 = 1.0 / pow(pele[e].dh, (ndgord + 1.0) / 2.0);

		unsigned int nType = ELEMTYPE;
#ifdef QUAD
		int *inode = pele[e].Node;
		if (inode[3] == inode[0])
			nType = 3;
#endif
		for (unsigned int ts = 0; ts < nType; ts++)
		{
			unsigned int s = pele[e].Side[ts]; // ��ts�ߵ�ȫ�ֱ��

			int e0 = psid[s].ofElem[0];
			if (e0 < 0) // ����BC_INTERBLK�߽� ���ڱ�ı߽������Ŀǰ���账����ֱ����һ���ߡ�
			{
				unsigned int nbctype = GetBCType(e0);

				if (nbctype == BC_INTERBLK)
					e0 = GetElement(e0);
				else
					continue;
			}

			int ts0 = psid[s].lsdidx[0];
			int e1 = psid[s].ofElem[1];
			int ts1 = psid[s].lsdidx[1];

			double outsdnorm[NDIM];
			if (e1 == ((int)e))
				for (unsigned int nd = 0; nd < NDIM; nd++)
					outsdnorm[nd] = psid[s].norm[nd];
			else
				for (unsigned int nd = 0; nd < NDIM; nd++)
					outsdnorm[nd] = -psid[s].norm[nd];

#if EQUATIO != EULNSEQ
			double outnorm = outsdnorm[0] * ppa->a + outsdnorm[1] * ppa->b;
#endif
#if EQUATIO == ADVECTI
			if (outnorm > 0.0)
				continue;
#endif
			double u[NEQ], du[NEQ];
			for (unsigned int j = 0; j < NEQ; j++)
			{
				u[j] = 0.0;
				du[j] = 0.0;
			}

			unsigned int uqdpts = psid[s].nqdpt;
			// each gauss quad pts
			unsigned int qdpos1 = uqdpts - 1;
			unsigned int qdpos0 = 0;
			for (unsigned int np = 0; np < uqdpts; np++, qdpos1--, qdpos0++)
			{
				double quadcoe = GauQd1dCoe[uqdpts][np][2];
				double tu1[NEQ], tu0[NEQ];

				if (nType == 3)
				{
					GetVar(ndgdof, pvar[e1].uh[irk], SdQdPtBasTri[ts1][uqdpts][qdpos1], tu1);
					GetVar(ndgdof, pvar[e0].uh[irk], SdQdPtBasTri[ts0][uqdpts][qdpos0], tu0);
				}
				else
				{
					GetVar(ndgdof, pvar[e1].uh[irk], SdQdPtBasRect[ts1][uqdpts][qdpos1], tu1);
					GetVar(ndgdof, pvar[e0].uh[irk], SdQdPtBasRect[ts0][uqdpts][qdpos0], tu0);
				}

#if EQUATIO == EULNSEQ
				if ((
					outsdnorm[0] * 0.5 * (tu1[1] / tu1[0] + tu0[1] / tu0[0]) +
					outsdnorm[1] * 0.5 * (tu1[2] / tu1[0] + tu0[2] / tu0[0])) > TOL)
					break;
#elif EQUATIO == BURGERS
				if (outnorm * 0.5 * (tu1[0] + tu0[0]) > TOL)
					break;
#endif
				for (unsigned int j = 0; j < NEQ; j++)
				{
					du[j] += (tu0[j] - tu1[j]) * quadcoe; // ���ŵ�˳����Ҫ����Ϊ����Ҫ����ֵ��
					u[j] = m_max(u[j], m_abs(tu0[0]));
					u[j] = m_max(u[j], m_abs(tu1[0]));
				}
			} // end of no pts

			for (unsigned int j = 0; j < NEQ; j++)
			{
				double isdiscon = m_abs(du[j]) * dhpp1d2 / (u[j] + TOL);
				pvar[e].isDiscon[j] = m_max(pvar[e].isDiscon[j], isdiscon);

				if (isdiscon > 1.0)
				{
					pvar[e].isDiscon[j] = 1;
					break;
				}
			}
		} // end of inside s

		for (unsigned int j = 0; j < NEQ; j++)
		{
			if (pvar[e].isDiscon[j] > 1.0)
				pvar[e].isDiscon[j] = 1.0;
			else
				pvar[e].isDiscon[j] = TOL;
		}
	} // end of elem
}

// positivity correction procedure
int PosiCorreProc(
    unsigned int irk, int it, double t, double cputm,
    unsigned int nele, unsigned int nsid, unsigned int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa)
{
	return 0;
}
