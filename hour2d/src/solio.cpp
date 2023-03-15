// file io
// 2011 07.28 liming
// write sol file
// save intermediate file
// read intermediate file

#include <math.h>
#include <stdio.h>
//#include <conio.h>
#include <string.h>
#include <stdlib.h>

#include "ic.h"
#include "bc2d.h"
#include "misc.h"
#include "solio.h"
#include "poly1d.h"
#include "poly2d.h"
#include "quadrature.h"
#include "m_math_def.h"

void OutSol(
    unsigned int irk, int it, double t, double cputm,
    unsigned int nele, unsigned int nsid, unsigned int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa)
{
	double(*pvavg)[NEQ]; // average
	pvavg = (double(*)[NEQ])calloc(nnod, sizeof(double) * NEQ);

	if (pvavg == NULL)
	{
		PRINTERRORN("calloc error!\n");
	}

	for (unsigned int e = 0; e < nele; e++)
	{
		unsigned int pfv = pvar[e].pfv;
		unsigned int nfvdof = DOFofOrd2d[pfv];

		unsigned int nType = ELEMTYPE;
		int *inode = pele[e].Node;
#ifdef QUAD
		if (inode[3] == inode[0])
			nType = 3;
#endif

		double(*EmQdBas)[MAXFVDOF];
		if (nType == 3)
			EmQdBas = (double(*)[MAXFVDOF])QdPtBasTri[0]; // tri
		else
			EmQdBas = (double(*)[MAXFVDOF])QdPtBasRect[0]; // rect

		for (unsigned int tn = 0; tn < nType; tn++)
		{
			double tcv[NEQ];
			double *tpv = tcv;
			GetVar(nfvdof, pvar[e].wh, EmQdBas[tn], tcv);

#if EQUATIO == EULNSEQ
#if VARIABLE == CONSERV
			double tmp[NEQ];
			tpv = tmp;
			// CV2PV( tcv, pvar[e].wh[0], tpv );
			double rho = tcv[0];
			if (rho < RMIN)
				rho = RMIN;
			if (rho > RMAX)
				rho = RMAX;

			// calculate pressure
			double pre = gamm1 * (tcv[NEQ - 1] - 0.5 * (tcv[1] * tcv[1] + tcv[2] * tcv[2]) / rho);
			if (pre < PMIN)
				pre = PMIN;
			if (pre > PMAX)
				pre = PMAX;

			tpv[0] = rho;
			tpv[1] = tcv[1] / rho;
			tpv[2] = tcv[2] / rho;
			tpv[NEQ - 1] = pre;
#endif
#endif
			for (unsigned int j = 0; j < NEQ; j++)
				pvavg[inode[tn]][j] += tpv[j];
		} // end of node of an element
	}	  // end of element

	for (unsigned int n = 0; n < nnod; n++)
		for (unsigned int j = 0; j < NEQ; j++)
		{
			pvavg[n][j] /= pnod[n].nNElem;
		}

	FILE *ssolu = NULL;
	// errno_t err;
	char fname[256] = "sol";

	if (ppa->isteady == UNSTEADY)
	{
		char buf[16];
		snprintf(buf, 16, "%d", it); // sprintf_s convert to snprintf,linux mode
		strcat(fname, buf);
	}

	strcat(fname, ".dat");
	/*
	if( ( err = fopen_s( &ssolu, fname, "w" ) ) !=0 )
	{
		free( pvavg );
		PRINTERRORN( "The file 'solution.dat' was not opened\n" );
	}*/
	ssolu = fopen(fname, "w"); // linux file open mode
	if (ssolu == NULL)
	{
		free(pvavg);
		PRINTERRORN("The file 'solution.dat' was not opened\n");
	}
	fprintf(ssolu, "TITLE = '2D Unstructured Hybrid Mesh'\n");
	fprintf(ssolu, "FILETYPE = FULL\n");

#if EQUATIO != EULNSEQ
	fprintf(ssolu, "VARIABLES = x, y, u\n");
#else
	fprintf(ssolu, "VARIABLES = x, y, R, u, v, P\n");
#endif

	fprintf(ssolu, "ZONE T = \"Zone var\"\n");

#ifdef TRI
	fprintf(ssolu, "NODES = %d, ELEMENTS = %d, ZONETYPE = FETRIANGLE\n", nnod, nele);
#else
	fprintf(ssolu, "NODES = %d, ELEMENTS = %d, ZONETYPE = FEQUADRILATERAL\n", nnod, nele);
#endif

	fprintf(ssolu, "DATAPACKING = POINT\n");

	for (unsigned int n = 0; n < nnod; n++)
	{
#if EQUATIO != EULNSEQ
		fprintf(ssolu, "% 18.15e % 18.15e % 18.15e\n",
			pnod[n].xy[0], pnod[n].xy[1], pvavg[n][0]);
#else
		fprintf(ssolu, "% 18.15e % 18.15e % 18.15e % 18.15e % 18.15e % 18.15e\n",
			pnod[n].xy[0], pnod[n].xy[1],
			pvavg[n][0], pvavg[n][1], pvavg[n][2], pvavg[n][NEQ - 1]); /// ppa->Preinf );
#endif
	}
	free(pvavg);

#ifdef TRI
	for (unsigned int e = 0; e < nele; e++)
	{
		unsigned int *inode = pele[e].Node;
		fprintf(ssolu, "% d % d % d\n", inode[0] + 1, inode[1] + 1, inode[2] + 1);
	}
#else
	for (unsigned int e = 0; e < nele; e++)
	{
		int *inode = pele[e].Node;
		fprintf(ssolu, "% d % d % d % d\n", inode[0] + 1, inode[1] + 1, inode[2] + 1, inode[3] + 1);
	}
#endif

	fclose(ssolu);
}

void OutWallSol(
    int irk, int it, double t, double cputm,
    int nbosd, int *pbosid,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa)
{
#if EQUATIO == EULNSEQ
	FILE *swall;
	unsigned int err;

	swall = fopen("wall.dat", "w");
	fprintf(swall, "VARIABLES = theta, x, y, Cp, entropy, Cf\n");

	double Cpx = 0.0;
	double Cpy = 0.0;
	double Ctx = 0.0;
	double Cty = 0.0;

	for (unsigned int itbd = 0; itbd < nbosd; ++itbd)
	{

		unsigned int ibd = pbosid[itbd];
		// printf("%d  %d  %d  %d\n",itbd,ibd,psid[ibd].Node[0], psid[ibd].Node[1]);
		if (GetBCType(psid[ibd].ofElem[0]) != BC_SOLIDSURFACE)
		{
			continue;
		}
		unsigned int eL = psid[ibd].ofElem[1];
		unsigned int tsL = psid[ibd].lsdidx[1];

		unsigned int nqdpt = psid[ibd].nqdpt;

		unsigned int nType = ELEMTYPE;
		int *inode = pele[eL].Node;
#ifdef QUAD
		if (inode[3] == inode[0])
			nType = 3;
#endif

		double(*SdQdPtBasL)[MAXFVDOF];
		if (nType == 3)
			SdQdPtBasL = (double(*)[MAXFVDOF])SdQdPtBasTri[tsL + 3][1];
		else
			SdQdPtBasL = (double(*)[MAXFVDOF])SdQdPtBasRect[tsL + 4][1];

		// ʹ�ò���������� ���ע�͵��˶δ�����ʹ��ֱ�߶�
		unsigned int nfvord = pvar[eL].pfv;
		unsigned int nfvdof = DOFofOrd2d[nfvord];

		double tcvL[NEQ];
		double gradcvL[NDIM][NEQ];
		GetVar(nfvdof, pvar[eL].wh, SdQdPtBasL[0], tcvL);
		GetVar(nfvdof, pvar[eL].grad[0], SdQdPtBasL[0], gradcvL[0]);
		GetVar(nfvdof, pvar[eL].grad[1], SdQdPtBasL[0], gradcvL[1]);

		double tx = 0.5 * (pnod[psid[ibd].Node[0]].xy[0] + pnod[psid[ibd].Node[1]].xy[0]);
		double ty = 0.5 * (pnod[psid[ibd].Node[0]].xy[1] + pnod[psid[ibd].Node[1]].xy[1]);

		double nx = psid[ibd].norm[0];
		double ny = psid[ibd].norm[1];

		double *tcvR = tcvL;

#if VARIABLE == CONSERV
		double tpv[NEQ];
		CV2PV(tcvR, tcvR, tpv);
#if AUXVARIABLE == CONSERV
		double gradpvL[NDIM][NEQ];
		CV2GradPV(tcvR, gradcvL, ppa->Tcoe, gradpvL);
#elif AUXVARIABLE == PRIMITI
		double(&gradpvL)[NDIM][NEQ] = gradcvL;
#endif
#elif VARIABLE == PRIMITI
		double *tpv = tcvR;
		double(&gradpvL)[NDIM][NEQ] = gradcvL;
#endif
		double T = tpv[3] / (tpv[0] * ppa->R);
		double miu = ppa->Re * CalMu(T, ppa->c4suth);

		double pre = tpv[NEQ - 1];
		double Tau[NDIM][NDIM];
		CalTau(gradpvL, miu, Tau);

		double dl = psid[ibd].dl;
		Cpx += pre * nx * dl;
		Cpy += pre * ny * dl;
		Ctx -= (Tau[0][0] * nx + Tau[0][1] * ny) * dl;
		Cty -= (Tau[1][0] * nx + Tau[1][1] * ny) * dl;

		double theta = atan2(ty, tx) * 180.0 / PI;

		double Cp = (tpv[NEQ - 1] - ppa->Preinf) * ppa->Rhoq2d2n1;
		double en = tpv[NEQ - 1] * pow(tpv[0], -gammaa);
		en = en * ppa->Sinfn1 - 1.0;
		double Cf = miu * (gradpvL[1][1]) * ppa->Rhoq2d2n1;
		fprintf(swall, "% .18f % .18f % .18f % .18f % .18f % .18f\n", theta, tx, ty, Cp, en, Cf);
	}

	double CDp = (Cpx * ppa->cosAoA + Cpy * ppa->sinAoA) * ppa->Rhoq2d2n1;
	double CLp = (-Cpx * ppa->sinAoA + Cpy * ppa->cosAoA) * ppa->Rhoq2d2n1;
	double CDt = (Ctx * ppa->cosAoA + Cty * ppa->sinAoA) * ppa->Rhoq2d2n1;
	double CLt = (-Ctx * ppa->sinAoA + Cty * ppa->cosAoA) * ppa->Rhoq2d2n1;

	FILE *stxt;
	if ((stxt = fopen("CDCL.dat", "a")) != 0)
	{
		fprintf(stxt, "%d %f, % .18f % .18f  % .18f % .18f  % .18f % .18f\n",
			it, t, CDp, CLp, CDt, CLt, CDp + CDt, CLp + CLt);
		fclose(stxt);
	}

	fclose(swall);
#else //����Euler����
	PRINTERRORN("wall function for Euler/NS eqs only!\n");
#endif
}

// save the data for restart
void SaveData(
    unsigned int irk, int it, double t, double cputm,
    double *res10, double *res20, double *resx0,
    unsigned int nele, uvar *pvar, para *ppa)
{
	FILE *sout = NULL;
	// errno_t err;
	/*
	if( ( err = fopen_s( &sout, "init.bin", "wb" ) ) !=0 )
	{
		PRINTERRORN( "File 'init.bin' was not opened\n" );
	}*/
	sout = fopen("init.bin", "wb"); // linux file open mode
	if (sout == NULL)
	{
		PRINTERRORN("File 'init.bin' was not opened\n");
	}
	fwrite(&it, sizeof(int), 1, sout);
	fwrite(&t, sizeof(double), 1, sout);
	fwrite(&cputm, sizeof(double), 1, sout);
	fwrite(res10, sizeof(double), NEQ, sout);
	fwrite(res20, sizeof(double), NEQ, sout);
	fwrite(resx0, sizeof(double), NEQ, sout);

	for (unsigned int e = 0; e < nele; e++)
	{
		unsigned int pdg = pvar[e].pdg;
		unsigned int ndgdof = DOFofOrd2d[pdg];
		fwrite(&pdg, sizeof(unsigned int), 1, sout);
		fwrite(pvar[e].uh[irk][0], sizeof(double), ndgdof * NEQ, sout);
	}

	fclose(sout);
}

// read the data for restart
void ReadData(
    unsigned int irk, int &it0, double &t0, double &cputm0,
    double *res10, double *res20, double *resx0,
    unsigned int nele, uvar *pvar, para *ppa)
{
	FILE *sini = NULL;
	// errno_t err;
	/*
	if( ( err = fopen_s( &sini, "init.bin", "rb" ) ) !=0 )
	{
		PRINTERRORN( "File 'init.bin' was not opened\n" );
	}*/
	sini = fopen("init.bin", "rb"); // linux mode
	if (sini == NULL)
	{
		PRINTERRORN("File 'init.bin' was not opened\n");
	}

	fread(&it0, sizeof(int), 1, sini);
	fread(&t0, sizeof(double), 1, sini);
	fread(&cputm0, sizeof(double), 1, sini);
	fread(res10, sizeof(double), NEQ, sini);
	fread(res20, sizeof(double), NEQ, sini);
	fread(resx0, sizeof(double), NEQ, sini);

	for (unsigned int e = 0; e < nele; e++)
	{
		pvar[e].pdg = ppa->pDG;
		pvar[e].pfv = ppa->pFV;
		unsigned int pdgread;
		fread(&(pdgread), sizeof(unsigned int), 1, sini);
		if (pdgread > MAXDGORD)
		{
			PRINTERRORN("poly read exceed the maxdgord, adjust the def and rebuild it.\n");
		}
		unsigned int ndofread = DOFofOrd2d[pdgread];
		fread(pvar[e].uh[irk][0], sizeof(double), ndofread * NEQ, sini);
	}

	fclose(sini);
}
