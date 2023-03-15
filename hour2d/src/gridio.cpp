// gridio.cpp
// ��ȡ���ָ�ʽ�������ļ�����ʱ֧�֡���
// 2011 11 14 rewrite liming.

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bc2d.h"
#include "gridio.h"
#include "gridparam.h"
#include "m_math_def.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "rcm.h"
#include "omp.h"
using namespace std;

/* read grid
   construct grid connections
   ��element_typeΪquad ����ʵ����������ʱ������ĸ���͵�һ��������ͬ
   return 1 fail! 0 ok����  */
int *ReadGridTetrex(
    const char *sGridFile,
    int &nele, int &nsid, int &nnod,
    elem *&pele, side *&psid, node *&pnod, uvar *&pvar,
    int &nbsd, int &nisd, int &ncvbd,
    int *&pbsd, int *&pisd, cvbd2d *&pcvbd, int rcm)
{
	FILE *sgrid = NULL;
	// errno_t err;

	// Open for read (will fail if file does not exist)
	/*
	if( (err = fopen_s( &sgrid, sGridFile, "r" )) !=0 )
	{
		printf( "%s", sGridFile );
		PRINTERRORE( " was not opened\n" );
	}*/
	sgrid = fopen(sGridFile, "r");
	if (sgrid == NULL)
	{
		printf("%s", sGridFile);
		PRINTERRORE(" was not opened\n");
	}

	char stemp[256];
	// ��ȡ���������Ϣ
	// 1
	// Tetrex Grid File exported by Gridgen 15.15R1 Fri Apr 27 10:06:55 2012
	// 1
	// ndimn    nzone   npoin      nvp
	for (unsigned int i = 0; i < 4; ++i)
		fgets(stemp, 256, sgrid);
	unsigned int ndimn, nzone, npoin, nvp;
	// ndimn 2d 2 3d 3,	// nzone number of zones
	// npoin number of pts,	// nvp ÿ����Ԫ�ı���Ŀ������3*TriEm+4*QuadEm,���ظ�
	fscanf(sgrid, "%u %u %u %u \n", &ndimn, &nzone, &npoin, &nvp);

	// support 2d only now
	if (ndimn != 2)
	{
		printf("%s", sGridFile);
		PRINTERRORE(" was not 2d\n");
	}

	nnod = npoin;
	nele = 0;
	nbsd = 0;
	// zone name     ncell   nbsurf
	fgets(stemp, 256, sgrid);
	for (unsigned int i = 0; i < nzone; ++i)
	{
		char zonename[64];
		unsigned int ncell, nbsurf;
		fscanf(sgrid, "%s", zonename, 64);
		fscanf(sgrid, "%u %u \n", &ncell, &nbsurf);

		nele += ncell;
		nbsd += nbsurf;
	}

	nsid = (nbsd + nvp) / 2; // ���������ϵı���Ŀ
	nisd = nsid - nbsd;

	printf("TETREX Grid Info:\n");
	printf("Number of Elem:       % 8d\n", nele);
	printf("Number of Side:       % 8d\n", nsid);
	printf("Number of Inner Side: % 8d\n", nisd);
	printf("Number of Bound Side: % 8d\n", nbsd);
	printf("Number of Node:       % 8d\n", nnod);
	// �����ڴ�
	pele = (elem *)calloc(nele, sizeof(elem));
	psid = (side *)calloc(nsid, sizeof(side));
	pnod = (node *)calloc(nnod, sizeof(node));
	pvar = (uvar *)calloc(nele, sizeof(uvar));

	if ((pele == NULL) || (psid == NULL) || (pnod == NULL) || (pvar == NULL))
	{
		PRINTERRORE("Insufficient memory available readgridtetrex\n");
	}

	// ��ȡ��Ԫ�ĵ��б�
	printf("\n2. read node list of elements\n");
	// volume number   volume type    point list
	fgets(stemp, 256, sgrid);
	for (unsigned int e = 0; e < nele; e++)
	{
		unsigned int volnum, voltyp;
		fscanf(sgrid, "%u %u", &volnum, &voltyp);
		int *inode = pele[volnum - 1].Node;

		if (voltyp == TETREXQUAD) // �ı���
		{
#ifdef TRI
			{
				PRINTERRORE("the elem read is quad but the program just for tri\n");
			}
#else
			// ˳ʱ�뻹����ʱ��ȡ����������������
			fscanf(sgrid, "%u %u %u %u\n", inode, inode + 1, inode + 2, inode + 3);

			for (unsigned int i = 0; i < 4; i++)
				inode[i]--;
#endif
#ifdef RECT
			if (inode[3] == inode[0])
			{
				PRINTERRORE("the elem read is nd3==nd0 but the program rect\n");
			}
#endif
		}
		else if (voltyp == TETREXTRI) // ������
		{
#ifdef RECT
			{
				PRINTERRORE("the elem read is tri but the program just for rect\n");
			}
#else
			// ˳ʱ�뻹����ʱ��ȡ����������������
			fscanf(sgrid, "%u %u %u\n", inode, inode + 1, inode + 2);

			for (unsigned int i = 0; i < 3; i++)
				inode[i]--;
#endif
#ifdef QUAD
			inode[3] = inode[0];
#endif
		}
		else
		{
			printf("%d ", voltyp);
			PRINTERRORE("unsupported elem type!\n");
		}
	}

	//��ȡ������λ����Ϣ
	printf("   read nodes' coordinates\n");
	// coordinates
	fgets(stemp, 256, sgrid);
	for (unsigned int n = 0; n < nnod; n++)
	{
		fscanf(sgrid, "%le %le\n", pnod[n].xy, pnod[n].xy + 1);
	}

	int *elemind = (int *)calloc(nele, sizeof(int));
	int *ofelemIndex = (int *)calloc(nsid * 2, sizeof(int));
	int *pbsd_inv = (int *)calloc(nsid, sizeof(int));

	if ((rcm == 1) || (rcm == 2) || (rcm == 3) || (rcm == 5))
	{
		printf("   Using RCM to renumber the node \n");
		GridRenumber_node(nele, nnod, pele, pnod);
		// printf("   Using RCM to renumber the elem \n");
		// elemind = GridRenumber_elem( nele, nnod, nsid, pele, pnod, psid );
		// GridRenumber_side( nele, nnod, nsid, pele, pnod, psid );
	}

	// �������
	printf("\n3. cal elem vol and xyc\n");
#pragma omp parallel for
	for (int e = 0; e < (int)nele; e++)
	{
		unsigned int nType = ELEMTYPE;
		int *inode = pele[e].Node;
#ifdef QUAD
		if (inode[3] == inode[0])
			nType = 3;
#endif
		double xx[4], yy[4];
		for (unsigned int i = 0; i < nType; ++i)
		{
			xx[i] = pnod[inode[i]].xy[0];
			yy[i] = pnod[inode[i]].xy[1];
		}

		// double area = area_of_tri( xx[0], yy[0], xx[1], yy[1], xx[2], yy[2] );
		double area = area_of_tri(xx, yy);

		if (nType == 3) // tri
		{
			pele[e].xyc[0] = (xx[0] + xx[1] + xx[2]) * 0.33333333333333333333333333333333;
			pele[e].xyc[1] = (yy[0] + yy[1] + yy[2]) * 0.33333333333333333333333333333333;
		}
		else // quad
		{
			double area1 = area;
			double area2 = area_of_tri(xx[0], yy[0], xx[2], yy[2], xx[3], yy[3]);
			double xc1 = (xx[0] + xx[1] + xx[2]) * 0.33333333333333333333333333333333;
			double yc1 = (yy[0] + yy[1] + yy[2]) * 0.33333333333333333333333333333333;
			double xc2 = (xx[0] + xx[2] + xx[3]) * 0.33333333333333333333333333333333;
			double yc2 = (yy[0] + yy[2] + yy[3]) * 0.33333333333333333333333333333333;
			area += area2;
			pele[e].xyc[0] = (xc1 * area1 + xc2 * area2) / area;
			pele[e].xyc[1] = (yc1 * area1 + yc2 * area2) / area;
		}

		if (area < TOL)
		{
			printf("elem %d %f\n", e, area);
			// PRINTERRORE( "area<tol\n" );
			getchar();

			// printf( "elem %d, area<tol\n", e );	getchar();return 1;
		}

		pele[e].volr = 1.0 / area;
	}

	printf("\n4. adjust node sequence of elem\n");
	unsigned int *pFirstNode = (unsigned int *)calloc(nele, sizeof(unsigned int));
	if (pFirstNode == NULL)
	{
		PRINTERRORE("Insufficient memory available readgridtetrex\n");
	}

// ��������Ԫ��Ŀ����
#pragma omp parallel for
	for (int n = 0; n < (int)nnod; ++n)
		// for( unsigned int n=0; n<nnod; n++ )
		pnod[n].nNElem = 0;

	// ��������Ԫ��Ŀ
	for (unsigned int e = 0; e < nele; ++e)
	{
		int *inode = pele[e].Node;
		unsigned int nType = ELEMTYPE;
#ifdef QUAD
		if (inode[3] == inode[0])
			nType = 3;
#endif
		for (unsigned int n = 0; n < nType; ++n)
			pnod[inode[n]].nNElem++;
	}

	// ����ofSide��ʱ�ڴ�
	for (unsigned int n = 0; n < nnod; ++n)
	{
		pnod[n].nNSide = 0;
		pnod[n].ofSide = (unsigned int *)calloc(pnod[n].nNElem + 1, sizeof(unsigned int));
		if (pnod[n].ofSide == NULL)
		{
			PRINTERRORE("Insufficient memory in readgridtetrex2d ofSide\n");
		}
	}

	// �������Ϣ
	unsigned int tsid = 0;
	int count = 0;
	if ((rcm == 2) || (rcm == 5))
	{

		for (unsigned int e = 0; e < nele; ++e)
		{
			unsigned int nType1 = ELEMTYPE;
			int *inode = pele[e].Node;
#ifdef QUAD
			if (inode[3] == inode[0])
				nType1 = 3;
#endif
			unsigned int tndlist1[8];
			for (unsigned int i = 0; i < nType1; ++i)
				tndlist1[i] = inode[i];
			tndlist1[nType1] = tndlist1[0];

			// ����e�ĸ�����
			for (unsigned int ts = 0; ts < nType1; ++ts)
			{
				unsigned int nd01 = tndlist1[ts];
				unsigned int nd11 = tndlist1[ts + 1];

				// ������nd0������
				int softs = 1; // û���ҵ��ı�־
				unsigned int idx = 0;
				for (; idx < pnod[nd01].nNSide; ++idx)
				{
					unsigned int sdofnd = pnod[nd01].ofSide[idx];
					unsigned int nd = psid[sdofnd].Node[1];
					if (nd == nd01)
						nd = psid[sdofnd].Node[0];
					if (nd == nd11)
					{ // ˵����(nd0,nd1)�Ѿ�����
						softs = 0;
						break;
					}
				}

				if (softs == 1)
				{				   // û�ҵ��� ts, then create side
					psid[tsid].Node[0] = nd11; // ע�⣬����ѡ��ʹ�ñ߷���͵�Ԫ���������
					psid[tsid].Node[1] = nd01;
					psid[tsid].ofElem[0] = -1;
					psid[tsid].ofElem[1] = e;
					pnod[nd11].ofSide[pnod[nd11].nNSide++] = tsid; // ���ӵ�nd1�Ա�tsid��������ϵ
					++tsid;
				}
				else
				{ // �ҵ��� ts ˵���˱��Ѿ�����һ����Ԫ��һ���ߣ������Ѿ�����ʼ��
					unsigned int s = pnod[nd01].ofSide[idx];

					psid[s].ofElem[0] = psid[s].ofElem[1]; // te1;
					psid[s].ofElem[1] = e;
					unsigned int tmpnd = psid[s].Node[0];
					psid[s].Node[0] = psid[s].Node[1];
					psid[s].Node[1] = tmpnd;
				}
			}
		}

		for (unsigned int i = 0; i < nsid; i++)
		{
			if (psid[i].ofElem[0] == -1)
			{
				ofelemIndex[count * 2] = psid[i].Node[0];
				ofelemIndex[count * 2 + 1] = psid[i].Node[1];
				count++;
			}
		}
		for (unsigned int n = 0; n < nnod; ++n)
		{
			pnod[n].nNSide = 0;
			for (unsigned int i = 0; i < pnod[n].nNElem + 1; i++)
			{
				pnod[n].ofSide[i] = 0;
			}
		}
		for (unsigned int i = 0; i < nsid; i++)
		{
			psid[i].ofElem[0] = 0;
			psid[i].ofElem[1] = 0;
			psid[i].Node[0] = 0;
			psid[i].Node[1] = 0;
		}
		Elemordernode(nele, nnod, nsid, pele, psid, elemind);
	}

	// �������Ϣ
	printf("\n5. creat sides\n");
	tsid = 0;
	for (unsigned int e = 0; e < nele; ++e)
	{
		unsigned int nType = ELEMTYPE;
		int *inode = pele[e].Node;
#ifdef QUAD
		if (inode[3] == inode[0])
			nType = 3;
#endif
		unsigned int tndlist[8];
		for (unsigned int i = 0; i < nType; ++i)
			tndlist[i] = inode[i];
		tndlist[nType] = tndlist[0];

		// ����e�ĸ�����
		for (unsigned int ts = 0; ts < nType; ++ts)
		{
			unsigned int nd0 = tndlist[ts];
			unsigned int nd1 = tndlist[ts + 1];

			// ������nd0������
			int softs = 1; // û���ҵ��ı�־
			unsigned int idx = 0;
			for (; idx < pnod[nd0].nNSide; ++idx)
			{
				unsigned int sdofnd = pnod[nd0].ofSide[idx];
				unsigned int nd = psid[sdofnd].Node[1];
				if (nd == nd0)
					nd = psid[sdofnd].Node[0];
				// unsigned int nd = psid[sdofnd].Node[0] + psid[sdofnd].Node[1] - nd0;

				// if( (nd0+nd1) == ( psid[sdofnd].Node[0]+ psid[sdofnd].Node[1]) )
				if (nd == nd1) // ˵����(nd0,nd1)�Ѿ�����
				{
					softs = 0;
					break;
				}
			}

			if (softs == 1) // û�ҵ��� ts, then create side
			{
				psid[tsid].Node[0] = nd1; // ע�⣬����ѡ��ʹ�ñ߷���͵�Ԫ���������
				psid[tsid].Node[1] = nd0;

				psid[tsid].ofElem[0] = -1;
				psid[tsid].ofElem[1] = e;
				psid[tsid].lsdidx[1] = ts;

				// pnod[nd0].ofSide[pnod[nd0].nNSide++] = tsid; // ���ӵ�nd0�Ա�tsid��������ϵ
				pnod[nd1].ofSide[pnod[nd1].nNSide++] = tsid; // ���ӵ�nd1�Ա�tsid��������ϵ

				pele[e].Side[ts] = tsid; //  ��Ԫe�ĵ�ts����Ϊ tsid
				++tsid;
			}
			else // �ҵ��� ts ˵���˱��Ѿ�����һ����Ԫ��һ���ߣ������Ѿ�����ʼ��
			{
				unsigned int s = pnod[nd0].ofSide[idx];

				psid[s].ofElem[0] = psid[s].ofElem[1]; // te1;
				psid[s].lsdidx[0] = psid[s].lsdidx[1];
				psid[s].ofElem[1] = e;
				psid[s].lsdidx[1] = ts;
				unsigned int tmpnd = psid[s].Node[0];
				psid[s].Node[0] = psid[s].Node[1];
				psid[s].Node[1] = tmpnd;

				pele[e].Side[ts] = s; // ��Ԫe�ĵ�ts����Ϊ s
			}
		} // end of ����e�ĸ����� ts

#ifdef QUAD
		if (inode[3] == inode[0])
			pele[e].Side[3] = pele[e].Side[0];
#endif
	} // end of e

	if ((rcm == 2) || (rcm == 5))
	{
		for (unsigned int i = 0; i < count; i++)
		{
			int a = ofelemIndex[i * 2];
			int b = ofelemIndex[i * 2 + 1];
			for (unsigned int j = 0; j < nsid; j++)
			{
				if (a == psid[j].Node[0])
				{
					if (b == psid[j].Node[1])
					{
						pbsd_inv[i] = j;
						break;
					}
				}
			}
			// printf("%d  %d\n",i,pbsd_inv[i]);
		}
	}
	free(ofelemIndex);

	if ((rcm == 3) || (rcm == 5))
	{
		printf("   Reorder the side at rank of node \n");
		SideReorder(nsid, nele, nnod, psid, pele, pnod);
	}

	/*printf("\n\n\nside      ofElem        Node\n");
	for(int i=0;i<nsid;i++){
		printf("%d     %d %d    %d %d\n",i+1,psid[i].ofElem[0],psid[i].ofElem[1],psid[i].Node[0],psid[i].Node[1]);
	}
	printf("\n\n\nelem      Side   \n");
	for(int i=0;i<nele;i++){
		printf( "%d      %d   %d   %d   %d\n",i+1,pele[i].Side[0]+1,pele[i].Side[1]+1,pele[i].Side[2]+1,pele[i].Side[3]+1 );
	}
	printf("\n\n\nelem      Node   \n");
	for(int i=0;i<nele;i++){
		printf( "%d      %d   %d   %d   %d\n",i+1,pele[i].Node[0]+1,pele[i].Node[1]+1,pele[i].Node[2]+1,pele[i].Node[3]+1 );
	}*/

	if (tsid != nsid)
	{
		PRINTERRORE("error!\n");
	}

	// ����߽����� �����߽�����
	printf("\n6. read bcs and modi bcs\n");
	ReadBCTetrex(sgrid, nele, nsid, nnod, pele, psid, pnod, pvar, nbsd, nisd, ncvbd, pbsd, pisd, pcvbd, elemind, rcm);
	free(pFirstNode);
	free(elemind);
	printf("   close tetrex grid file\n");
	fclose(sgrid);

	// ��������������Ϣ ������
	printf("\n7. modi node's lists and sort it\n");
	SortListofNd(nele, nsid, nnod, pele, psid, pnod);

	// CalElemWhtToNd( nele, nsid, nnod, pele, psid, pnod, pvar );

	printf("leave readgridtetrex\n");

	return pbsd_inv;
}

// ����߽�����
void ReadBCTetrex(
    FILE *sgrid,
    int nele, int nsid, int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar,
    int nbsd, int nisd, int ncvbd,
    int *&pbsd, int *&pisd, cvbd2d *pcvbd, int *elemind, int rcm)
{
	char stemp[256];
	// bc name   volume number   local surface number
	fgets(stemp, 256, sgrid);

	int numofBCs[N_BCTYPE];
	for (unsigned int i = 0; i < N_BCTYPE; ++i)
		numofBCs[i] = 0;

	for (unsigned int ts = 0; ts < nbsd; ++ts)
	{
		char sbctype[16];
		unsigned int volnum, locsur;
		fscanf(sgrid, "%s %u %u\n", sbctype, &volnum, &locsur);
		if ((rcm == 2) || (rcm == 5))
		{
			volnum = elemind[volnum - 1] + 1;
		}

		int nbctype = 0;
		unsigned int idx = 0;
		for (; idx < N_BCTYPE; ++idx)
		{
			if (strstr(strTetrexBC[idx], sbctype) != NULL)
			// if( strstr( sbctype, strTetrexBC[idx] ) != NULL )
			{
				nbctype = nTetrexBC[idx];
				break;
			}
		}

		if (idx == 0)
		{
			printf("%s %u %u\n", sbctype, volnum - 1, locsur - 1);
		}

		if (idx == N_BCTYPE)
		{
			/*printf( "unknown bc types %s\n", sbctype );
			getchar();
			return;*/
			printf("%s ", sbctype);
			PRINTERRORE("unknown bc types\n");
		}

#if EQUATIO == ADVECTI
		// ��������ʹ��inflow��outflow���������ڱ߽�����BC_INTERBLK BndCond_6
		// 1) inflow ���� g(x,y,t)��outflow �Ϸ����� �Զ��ж���in����out
		// 2) ����
		// nbctype = BC_INTERBLK;
		nbctype = BC_FARFIELD;
#elif EQUATIO == BURGERS
		// BURGERS����ʹ��inflow��outflow���������ڱ߽�����BC_INTERBLK BndCond_6
		// 1) inflow ���� g(x,y,t)��outflow �Ϸ����� �Զ��ж���in����out
		//    �˺���һ���ǵ��������ģ�����ʱ�侫�ȿ������Ե����ף���Ҫ��߽�ʱ����׵�����
		// 2) ����
		nbctype = BC_INTERBLK;
#elif EQUATIO == EULNSEQ
		// ����
		// if( nbctype != BC_SOLIDSURFACE )
		//		if( nbctype==BC_SYMMETRY )
		//			nbctype = BC_INFLOW;//BC_OUTFLOW;//BC_FARFIELD; //;
		// nbctype = BC_FARFIELD;
		if (nbctype == BC_INFLOW)
			nbctype = BC_FARFIELD;
		else if (nbctype == BC_OUTFLOW)
			nbctype = BC_FARFIELD;
			//		if( nbctype == BC_INFLOW )
			//			nbctype = BC_GENERIC1;
			//		if( nbctype == BC_OUTFLOW )
			//			nbctype = BC_SOLIDSURFACE;

			// ����������
			//		nbctype = BC_INTERBLK;
			//		nbctype = BC_OUTFLOW;
			// if( nbctype == BC_INFLOW )
			//	nbctype = BC_FARFIELD;
			// if( nbctype == BC_OUTFLOW )
			//	nbctype = BC_FARFIELD;

			// if( nbctype == BC_GENERIC1 )
			//	nbctype = BC_FARFIELD;

			// Couette
			/*if( nbctype == BC_INFLOW )
				nbctype = BC_INTERBLK;
			if( nbctype == BC_OUTFLOW )
				nbctype = BC_INTERBLK;
			if( nbctype == BC_SYMMETRY )
				nbctype = BC_GENERIC1;*/
			// if( nbctype == BC_SOLIDSURFACE )
			//	nbctype = BC_GENERIC1;

			// double mach
			// if( nbctype==BC_FARFIELD )
			//	nbctype = BC_INFLOW;
#endif

		// if( nbctype==BC_NONE ) { PRINTERRORE( "BC_NONE!\n" ); }

		++numofBCs[nbctype];

		int elemIdx = volnum - 1;
		unsigned int bdlocsd = locsur - 1;
		unsigned int nType = ELEMTYPE;

#ifdef QUAD
		int *inode = pele[elemIdx].Node;
		if (inode[3] == inode[0])
			nType = 3;
#endif

		for (unsigned int ts = 0; ts < nType; ++ts)
		{
			if (psid[pele[elemIdx].Side[ts]].ofElem[0] == -1)
			{
				// psid[pele[elemIdx].Side[ts]].ofElem[0] = SetBCType( elemIdx, nbctype );
				psid[pele[elemIdx].Side[ts]].ofElem[0] = SetBCType(0, nbctype);
				// break;
			}
		}
	}

	for (unsigned int i = 0; i < 9; ++i)
	{
		printf("%8d %s\n", numofBCs[i], strTetrexBC[i]);
	}

	pbsd = (int *)calloc(nbsd, sizeof(int));
	pisd = (int *)calloc(nisd, sizeof(int));

	if ((pbsd == NULL) || (pisd == NULL))
	{
		PRINTERRORE("Insufficient memory available pbsd pisd\n");
	}

	unsigned int tnbsd = 0;
	unsigned int tnisd = 0;
	for (unsigned int s = 0; s < nsid; ++s)
	{
		int e0 = psid[s].ofElem[0];

		if (e0 < 0)
			pbsd[tnbsd++] = s;
		else
			pisd[tnisd++] = s;
	}

	printf("bdsd %d =? %d, insd %d =? %d\n", nbsd, tnbsd, nisd, tnisd);
}

inline double sign_of_area_elem1(elem *pele, side *psid, node *pnod, unsigned int s, unsigned int nd0)
{
	unsigned int nd1 = psid[s].Node[0];
	if (nd1 == nd0)
		nd1 = psid[s].Node[1];

	int elem1 = psid[s].ofElem[1];

	unsigned int nd2;
	for (unsigned int i = 0; i < 3; ++i)
	{
		nd2 = pele[elem1].Node[i];
		if ((nd2 != nd0) && (nd2 != nd1))
			break;
	}

	return (pnod[nd1].xy[0] - pnod[nd0].xy[0]) * (pnod[nd2].xy[1] - pnod[nd0].xy[1]) -
	       (pnod[nd1].xy[1] - pnod[nd0].xy[1]) * (pnod[nd2].xy[0] - pnod[nd0].xy[0]);
}

// �Ե������߼���Ԫ����������ʱ�뷽��
// ���ڱ߽�㣬����ֻ����nNSide=nNElem+1
void SortListofNd(
    int nele, int nsid, int nnod,
    elem *pele, side *psid, node *pnod)
{
	// ���ԭ����Ϣ
	for (unsigned int n = 0; n < nnod; ++n)
	{
		pnod[n].nNSide = 0;
		pnod[n].nNElem = 0;
		if (pnod[n].ofSide != NULL)
		{
			free(pnod[n].ofSide);
			pnod[n].ofSide = NULL;
		}
		if (pnod[n].ofElemWht != NULL)
		{
			free(pnod[n].ofElemWht);
			pnod[n].ofElemWht = NULL;
		}
	}

	// �ҵ���������Ԫ��Ŀ
	for (unsigned int e = 0; e < nele; ++e)
	{
		unsigned int nType = ELEMTYPE;
		int *inode = pele[e].Node;
#ifdef QUAD
		if (inode[3] == inode[0])
			nType = 3;
#endif
		for (unsigned int n = 0; n < nType; ++n)
			pnod[inode[n]].nNElem++;
	}
	// �ҵ�����������Ŀ
	for (unsigned int s = 0; s < nsid; ++s)
	{
		pnod[psid[s].Node[0]].nNSide++;
		pnod[psid[s].Node[1]].nNSide++;
	}

	// �����ڴ棬�������ߺ͵�������Ԫ�ڴ�һ������
	// ������ָ��ͷ��������Ԫָ���м�ĳ��λ��
	// ���ڱ߽�㣬����ֻ����nNSide=nNElem+1������������
	for (unsigned int n = 0; n < nnod; ++n)
	{
		pnod[n].ofSide = (unsigned int *)
		    calloc(pnod[n].nNSide + pnod[n].nNElem, sizeof(unsigned int));
		pnod[n].ofElemWht = (double *)calloc(pnod[n].nNElem, sizeof(double));

		if ((pnod[n].ofSide == NULL) || (pnod[n].ofElemWht == NULL))
		{
			PRINTERRORE("Insufficient memory available sortndlists\n");
		}

		pnod[n].ofElem = pnod[n].ofSide + pnod[n].nNSide;
	}

	// �ٴ����
	for (unsigned int n = 0; n < nnod; ++n)
		//	{
		pnod[n].nNSide = 0;

	// �ҵ���������
	for (unsigned int s = 0; s < nsid; ++s)
	{
		pnod[psid[s].Node[0]].ofSide[pnod[psid[s].Node[0]].nNSide++] = s;
		pnod[psid[s].Node[1]].ofSide[pnod[psid[s].Node[1]].nNSide++] = s;
	}

	// each node's elem. list
	for (unsigned int n = 0; n < nnod; ++n)
	{
		unsigned int tsdlist[32];
		int temlist[32];

		if ((pnod[n].nNSide) == (pnod[n].nNElem)) // inner node i.c.
		{
			tsdlist[0] = pnod[n].ofSide[0]; // ���һ������Ϊ��ʼ
			double tarea2 = sign_of_area_elem1(pele, psid, pnod, tsdlist[0], n);

			if (tarea2 > TOL)
				temlist[0] = psid[tsdlist[0]].ofElem[1];
			else
				temlist[0] = psid[tsdlist[0]].ofElem[0]; // �˴��޷���֤ > TOL
		}
		else // boundary node i.c.
		{
			unsigned int k = 0;
			for (; k < pnod[n].nNSide; ++k)
			{
				unsigned int s = pnod[n].ofSide[k];
				if (psid[s].ofElem[0] >= 0)
					continue; // ��s�߲��Ǳ߽��

				double tarea2 = sign_of_area_elem1(pele, psid, pnod, s, n);

				if (tarea2 > TOL)
				{
					tsdlist[0] = s;
					temlist[0] = psid[s].ofElem[1];
					break;
				}
			} // ������б�

			if (k == pnod[n].nNSide)
			{
				printf("sort list Error at boundary point %d\n", n);
				getchar();
				return;
			}
		} // end of boundary node i.c.

		for (unsigned int k = 0; k < pnod[n].nNElem; ++k) // start fill the sorted list !!
		{
			if (k == (pnod[n].nNSide - 1))
				break; // �ڵ�ʱ���һ��k�Ƕ��

			unsigned int nType = ELEMTYPE;
#ifdef QUAD
			int *inode = pele[temlist[k]].Node;
			if (inode[3] == inode[0])
				nType = 3;
#endif
			unsigned int idx = 0;
			for (; idx < nType; ++idx) // each side of element
			{
				unsigned int s = pele[temlist[k]].Side[idx];

				if ((s != tsdlist[k]) && ((psid[s].Node[0] == n) || (psid[s].Node[1] == n))) // find another side of element to node n;
				{
					tsdlist[k + 1] = s;
					temlist[k + 1] = psid[s].ofElem[1];
					if (temlist[k + 1] == temlist[k])
						temlist[k + 1] = psid[s].ofElem[0];
					break;
				}
			}

			if (idx == nType)
			{
				printf("Error at boundary point %d\n", n);
				getchar();
				return;
			}
		}

		for (unsigned int k = 0; k < pnod[n].nNSide; ++k)
			pnod[n].ofSide[k] = tsdlist[k];
		for (unsigned int k = 0; k < pnod[n].nNElem; ++k)
			pnod[n].ofElem[k] = temlist[k];
	} // end of each node in mesh.
}

// return 1 false
// return 0 ok!
int CheckGrid(
    int nele, int nsid, int nnod,
    const elem *pele, const side *psid, const node *pnod)
{
	for (unsigned int s = 0; s < nsid; ++s) // ����ÿ����
	{
		for (unsigned int tn = 0; tn < 2; ++tn) // �����
		{
			unsigned int n = psid[s].Node[tn];
			unsigned int tns = pnod[n].nNSide;

			if (tns == 0)
				PRINTERROR1("��û�кͱ�����!\n");

			unsigned int k = 0;
			for (; k < tns; ++k)
			{
				if (pnod[n].ofSide[k] == s)
					break;
			}

			if (k == tns)
				PRINTERROR1("��ͱߵ�������ϵ����!\n");
		} // end of n of side

		for (unsigned int te = 0; te < 2; ++te) // ������
		{
			int e = psid[s].ofElem[te];

			if (e < 0)
				continue;

			unsigned int nType = ELEMTYPE;
#ifdef QUAD
			const int *inode = pele[e].Node;
			if (inode[3] == inode[0])
				nType = 3;
#endif
			unsigned int k = 0;
			for (; k < nType; ++k)
			{
				if (pele[e].Side[k] == s)
					break;
			}

			if (k == nType)
				PRINTERROR1("��ͱ����ڹ�ϵ����!\n");
		} // end of e of side
	}	  // end of each side

	return 0;
}
