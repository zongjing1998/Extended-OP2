// left hand side
// 2011 07 28 liming
// matrix and inverse of matrix etc

#include <math.h>
#include <time.h>
//#include <conio.h>
#include <stdio.h>
#include <stdlib.h>

#include "bc2d.h"
#include "def2d.h"
#include "gridparam.h"
#include "quadrature.h"

#include "m_math_def.h"

// cal side para dl norm check something
void CalSidePara(
    unsigned int nele, unsigned int nsid, unsigned int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa)
{
	printf("calculate side parameters\n");

	for (unsigned int s = 0; s < nsid; ++s)
	{
		psid[s].nqdpt = 1 + ppa->pFV; // ppa->PDG;//

		unsigned int nd0 = psid[s].Node[0];
		unsigned int nd1 = psid[s].Node[1];

		double dxy[NDIM];
		for (unsigned int tn = 0; tn < NDIM; ++tn)
		{
			dxy[tn] = pnod[nd1].xy[tn] - pnod[nd0].xy[tn];
		}

		psid[s].dl = sqrt(dxy[0] * dxy[0] + dxy[1] * dxy[1]);

		psid[s].norm[0] = -dxy[1] / (psid[s].dl + TOL);
		psid[s].norm[1] = dxy[0] / (psid[s].dl + TOL);

		int e1 = psid[s].ofElem[1];
		int e0 = psid[s].ofElem[0];

		unsigned int nType1 = ELEMTYPE;
#ifdef QUAD
		int *inode1 = pele[e1].Node;
		if (inode1[3] == inode1[0])
			nType1 = 3;
#endif
		if (e0 < 0)
		{
			if (nType1 == 3)
				psid[s].etaf = 1.5; // tri
			else
				psid[s].etaf = 2.0; // rect
		}
		else
		{
			unsigned int nType0 = ELEMTYPE;
#ifdef QUAD
			int *inode0 = pele[e0].Node;
			if (inode0[3] == inode0[0])
				nType0 = 3;
#endif
			double ar1 = 1.0 / pele[e1].volr;
			double ar0 = 1.0 / pele[e0].volr;

			double tmp;
			if (nType1 == 3)
				tmp = ar1 * 2.0 / 3.0;
			else
				tmp = ar1 * 0.5;
			if (nType0 == 3)
				tmp += ar0 * 2.0 / 3.0;
			else
				tmp += ar0 * 0.5;

			psid[s].etaf = 4.0 * ar1 * ar0 / (ar1 + ar0) / tmp;
		}

		// ���¶�����ע�͵�
		double nx = psid[s].norm[0];
		double ny = psid[s].norm[1];

		unsigned int tmnd;
		// �ҵ������λ��ı��ε�Ԫps[s].ofElem[1]�Ǵ�����֮��ĵ�
		for (unsigned int tn = 0; tn < 3; ++tn)
		{
			tmnd = pele[psid[s].ofElem[1]].Node[tn];

			if ((tmnd != nd0) && (tmnd != nd1))
				break;
		}

		// ʹ�÷�����ofElem[1]���ⷨ��
		if ((
			nx * (pnod[tmnd].xy[0] - pnod[nd0].xy[0]) +
			ny * (pnod[tmnd].xy[1] - pnod[nd0].xy[1])) > TOL)
		{
			printf("calsidepara norm error!\n");
			getchar();
			psid[s].norm[0] = -nx;
			psid[s].norm[1] = -ny;
		}

		// int e1 = psid[s].ofElem[1];
		unsigned int nType = ELEMTYPE;
#ifdef QUAD
		int *inode = pele[e1].Node;
		if (inode[3] == inode[0])
			nType = 3;
#endif
		unsigned int ts1 = 0;
		// �ҵ��˱���1��Ԫ�ľֲ����
		for (; ts1 < nType; ts1++)
		{
			if (pele[e1].Side[ts1] == s)
				break;
		}

		if (ts1 != psid[s].lsdidx[1]) // ||(ts1==nType) )
		{
			printf("calsidepara lsdidx1 error!\n");
			getchar();
		}

		// int e0 = psid[s].ofElem[0];
		if (e0 < 0)
			continue;
		nType = ELEMTYPE;
#ifdef QUAD
		inode = pele[e0].Node;
		if (inode[3] == inode[0])
			nType = 3;
#endif
		unsigned int ts0 = 0;
		// �ҵ��˱���0��Ԫ�ľֲ����
		for (; ts0 < nType; ts0++)
		{
			if (pele[e0].Side[ts0] == s)
				break;
		}

		if (ts0 != psid[s].lsdidx[0])
		{
			printf("calsidepara lsdidx0 error!\n");
			getchar();
		}
	} // end side s
}

// cal element para
void CalElemPara(
    unsigned int nele, unsigned int nsid, unsigned int nnod,
    elem *pele, side *psid, node *pnod, uvar *pvar, para *ppa)
{
	printf("calculate element parameters\n");

	for (unsigned int e = 0; e < nele; ++e)
	{
		int *inode = pele[e].Node;
		unsigned int nType = ELEMTYPE;
#ifdef QUAD
		if (inode[3] == inode[0])
			nType = 3;
#endif
		unsigned int nd0 = inode[0];
		unsigned int nd1 = inode[1];
		unsigned int nd23;
		double twodn1;

		if (nType == 3)
		{
			nd23 = inode[2];
			twodn1 = 0.5 * pele[e].volr;
		}
		else
		{
			nd23 = inode[3];
			twodn1 = 2.0 * pele[e].volr;
		}

		pele[e].LocdGlb[0][0] = (pnod[nd23].xy[1] - pnod[nd0].xy[1]) * twodn1;
		pele[e].LocdGlb[0][1] = -(pnod[nd23].xy[0] - pnod[nd0].xy[0]) * twodn1;
		pele[e].LocdGlb[1][0] = -(pnod[nd1].xy[1] - pnod[nd0].xy[1]) * twodn1;
		pele[e].LocdGlb[1][1] = (pnod[nd1].xy[0] - pnod[nd0].xy[0]) * twodn1;

		// max dl of element triangle
		double maxdl = psid[pele[e].Side[0]].dl;
		pele[e].dxy[0] = m_abs(pnod[nd23].xy[0] - pnod[nd0].xy[0]);
		pele[e].dxy[1] = m_abs(pnod[nd23].xy[1] - pnod[nd0].xy[1]);

		for (unsigned int i = 1; i < nType; i++)
		{
			if (psid[pele[e].Side[i]].dl > maxdl)
				maxdl = psid[pele[e].Side[i]].dl;

			double tmp = m_abs(pnod[inode[i]].xy[0] - pnod[inode[i - 1]].xy[0]);

			if (tmp > pele[e].dxy[0])
				pele[e].dxy[0] = tmp;

			tmp = m_abs(pnod[inode[i]].xy[1] - pnod[inode[i - 1]].xy[1]);

			if (tmp > pele[e].dxy[1])
				pele[e].dxy[1] = tmp;
		}

		if (nType == 3)
			pele[e].dh = 2.0 / (pele[e].volr * maxdl);
		else
			pele[e].dh = 1.0 / (pele[e].volr * maxdl);
	}
}
