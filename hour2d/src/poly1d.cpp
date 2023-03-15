// poly 1d
// 2011 07 28 liming
// polynomials and its derivative

#include "poly1d.h"

#include <math.h>

////lp0
inline double lp0(double xi)
{
	return 1.0;
}

inline double lp0d1(double xi)
{
	return 0.0;
}
////lp1
inline double lp1(double xi)
{
	return xi;
}

inline double lp1d1(double xi)
{
	return 1.0;
}
////lp2
inline double lp2(double xi, double xi2)
{
	// return 0.5*(3.0*xi*xi-1.0);
	// return 1.5*xi*xi-0.5;
	return 1.5 * xi2 - 0.5;
}

inline double lp2d1(double xi, double xi2)
{
	return 3.0 * xi;
}
////lp3
inline double lp3(double xi, double xi2)
{
	// return 0.5*xi*(5.0*xi2-3.0);
	return xi * (2.5 * xi2 - 1.5);
}

inline double lp3d1(double xi, double xi2)
{
	// return (7.5*xi*xi-1.5);
	return (7.5 * xi2 - 1.5);
}
////lp4
inline double lp4(double xi, double xi2)
{
	// return 0.125*( (35.0*xi2-30.0)*xi2 + 3.0 );
	return ((4.375 * xi2 - 3.75) * xi2 + 0.375);
}

inline double lp4d1(double xi, double xi2)
{
	return xi * (17.5 * xi2 - 7.5);
}

double BasisPoly1d(unsigned int nord, double xi)
{
	double xi2;
	switch (nord)
	{
	case 0:
		return 1.0;
		break;
	case 1:
		return xi * 1.7320508075688772935274463415059;
		break;
	case 2:
		xi2 = xi * xi;
		return lp2(xi, xi2) * 2.2360679774997896964091736687313;
		break;
	case 3:
		xi2 = xi * xi;
		return lp3(xi, xi2) * 2.6457513110645905905016157536393;
		break;
	case 4:
		xi2 = xi * xi;
		return lp4(xi, xi2) * 3.0;
		break;
	default:
		return 0.0;
	}
}

double BasisPolyDe1d(unsigned int nord, double xi)
{
	double xi2;
	switch (nord)
	{
	case 0:
		return 0.0;
		break;
	case 1:
		return 1.7320508075688772935274463415059;
		break;
	case 2:
		xi2 = xi * xi;
		return lp2d1(xi, xi2) * 2.2360679774997896964091736687313;
		break;
	case 3:
		xi2 = xi * xi;
		return lp3d1(xi, xi2) * 2.6457513110645905905016157536393;
		break;
	case 4:
		xi2 = xi * xi;
		return lp4d1(xi, xi2) * 3.0;
		break;
	default:
		return 0.0;
	}
}

//             use 0-MAXSEGQDPTS:MAXSEGQDPTS :8
double QdPtCoe1d[MAXSEGQDPTS + 1][MAXSEGQDPTS][2];
double QdPtCoe1d3[MAXSEGQDPTS + 1][MAXSEGQDPTS][3];
double QdPtBas1d[MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVORD + 1];
double QdPtBasDe1d[MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVORD + 1];
int QdNPt2NOrd1d[MAXSEGQDPTS + 1];
int QdNOrd2NPt1d[MAXSEGQDPTS + MAXSEGQDPTS];

//                use 0-MAXSEGQDPTS :MAXSEGQDPTS*MAXSEGQDPTS
double QdPtCoeRect[MAXSEGQDPTS + 1][MAXSEGQDPTS * MAXSEGQDPTS][4];
double QdPtBasRect[MAXSEGQDPTS + 1][MAXSEGQDPTS * MAXSEGQDPTS][MAXFVDOF];
double QdPtBasDxiRect[MAXSEGQDPTS + 1][MAXSEGQDPTS * MAXSEGQDPTS][MAXFVDOF];
double QdPtBasDetRect[MAXSEGQDPTS + 1][MAXSEGQDPTS * MAXSEGQDPTS][MAXFVDOF];
//                     use 1-MAXSDQDPTS  :MAXSEGQDPTS
double SdQdPtBasRect[8][MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVDOF];
double SdQdPtBasDxiRect[8][MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVDOF];
double SdQdPtBasDetRect[8][MAXSEGQDPTS + 1][MAXSEGQDPTS][MAXFVDOF];
double SdQdPtRectXi[8][MAXSEGQDPTS + 1][MAXSEGQDPTS];
double SdQdPtRectEt[8][MAXSEGQDPTS + 1][MAXSEGQDPTS];

void CalBasArr1d(para *ppa)
{
	QdNPt2NOrd1d[0] = -1; // not used
	for (unsigned int nqdpt = 1; nqdpt <= MAXSEGQDPTS; nqdpt++)
	{
		int qdord = nqdpt + nqdpt - 1;
		QdNPt2NOrd1d[nqdpt] = qdord;
	}

	for (int i = 0; i < MAXSEGQDPTS; i++)
	{
		int i2 = i + i;
		QdNOrd2NPt1d[i2] = QdNOrd2NPt1d[i2 + 1] = i + 1;
	}

	for (unsigned int nqdpt = 1; nqdpt <= MAXSEGQDPTS; nqdpt++)
		for (unsigned int qd = 0; qd < nqdpt; qd++) // at each gauss points
		{
			QdPtCoe1d[nqdpt][qd][0] = GauQd1dCoe2[nqdpt][qd][0];
			QdPtCoe1d[nqdpt][qd][1] = GauQd1dCoe2[nqdpt][qd][1];

			double xi = QdPtCoe1d[nqdpt][qd][0];

			QdPtCoe1d3[nqdpt][qd][0] = 0.5 - 0.5 * xi;
			QdPtCoe1d3[nqdpt][qd][1] = 0.5 + 0.5 * xi;
			QdPtCoe1d3[nqdpt][qd][2] = QdPtCoe1d[nqdpt][qd][1];

			for (unsigned int i = 0; i <= MAXFVORD; i++)
			{
				QdPtBas1d[nqdpt][qd][i] = BasisPoly1d(i, xi);
				QdPtBasDe1d[nqdpt][qd][i] = BasisPolyDe1d(i, xi);
			}
		} // end of qd

	// nqdpt==0
	QdPtCoe1d[0][0][0] = -1.0;
	QdPtCoe1d[0][0][1] = 0.5;
	QdPtCoe1d[0][1][0] = 1.0;
	QdPtCoe1d[0][1][1] = 0.5;
	// QdPtCoe1d[0][2][0] = -1.0; QdPtCoe1d[0][2][1] =  0.5;
	QdPtCoe1d[0][3][0] = 1.0; // QdPtCoe1d[0][3][1] =  0.5; for 1d program only

	QdPtCoe1d3[0][0][0] = 1.0;
	QdPtCoe1d3[0][0][1] = 0.0;
	QdPtCoe1d3[0][0][2] = 0.5;
	QdPtCoe1d3[0][1][0] = 0.0;
	QdPtCoe1d3[0][1][1] = 1.0;
	QdPtCoe1d3[0][1][2] = 0.5;

	// GetBas1d  ( MAXFVORD+1, -1.0, QdPtBas1d  [0][0] );
	// GetBas1d  ( MAXFVORD+1,  1.0, QdPtBas1d  [0][1] );
	// GetBasDe1d( MAXFVORD+1, -1.0, QdPtBasDe1d[0][0] );
	// GetBasDe1d( MAXFVORD+1,  1.0, QdPtBasDe1d[0][1] );
	for (unsigned int i = 0; i <= MAXFVORD; i++)
	{
		QdPtBas1d[0][0][i] = BasisPoly1d(i, -1.0);
		QdPtBas1d[0][1][i] = BasisPoly1d(i, 1.0);
		QdPtBasDe1d[0][0][i] = BasisPolyDe1d(i, -1.0);
		QdPtBasDe1d[0][1][i] = BasisPolyDe1d(i, 1.0);
		// QdPtBas1d  [0][2][i] = BasisPoly1d  ( i, -1.0 );
		// QdPtBas1d  [0][3][i] = BasisPoly1d  ( i,  1.0 ) * 0.5;
		// QdPtBasDe1d[0][2][i] = BasisPolyDe1d( i, -1.0 );
		// QdPtBasDe1d[0][3][i] = BasisPolyDe1d( i,  1.0 ) * 0.5;
		QdPtBas1d[0][3][i] = QdPtBas1d[0][1][i] * 0.5; // for 1d program only
		QdPtBasDe1d[0][3][i] = QdPtBasDe1d[0][1][i] * 0.5;
	}
}

void CalBasArrRect()
{
	// rect side qd 0(first node for out sol) 1 2 ...
	for (unsigned int nqdpt = 0; nqdpt <= MAXSEGQDPTS; nqdpt++)
	{
		unsigned int npnum = nqdpt;
		if (npnum == 0)
			npnum = 1;
		// at each gauss points
		for (unsigned int qd = 0; qd < npnum; qd++)
		{
			// double xiet = GauQd1d Coe2[nqdpt][qd][0];
			double xiet = QdPtCoe1d[nqdpt][qd][0];

			SdQdPtRectXi[0][nqdpt][qd] = xiet;
			SdQdPtRectEt[0][nqdpt][qd] = -1.0;
			SdQdPtRectXi[4][nqdpt][qd] = -xiet;
			SdQdPtRectEt[4][nqdpt][qd] = -1.0;
			SdQdPtRectXi[1][nqdpt][qd] = 1.0;
			SdQdPtRectEt[1][nqdpt][qd] = xiet;
			SdQdPtRectXi[5][nqdpt][qd] = 1.0;
			SdQdPtRectEt[5][nqdpt][qd] = -xiet;
			SdQdPtRectXi[2][nqdpt][qd] = -xiet;
			SdQdPtRectEt[2][nqdpt][qd] = 1.0;
			SdQdPtRectXi[6][nqdpt][qd] = xiet;
			SdQdPtRectEt[6][nqdpt][qd] = 1.0;
			SdQdPtRectXi[3][nqdpt][qd] = -1.0;
			SdQdPtRectEt[3][nqdpt][qd] = -xiet;
			SdQdPtRectXi[7][nqdpt][qd] = -1.0;
			SdQdPtRectEt[7][nqdpt][qd] = xiet;

			for (int j = 0; j < 8; j++)
			{
				Get1dTensorArbi(
				    MAXFVORD, SdQdPtRectXi[j][nqdpt][qd], SdQdPtRectEt[j][nqdpt][qd],
				    SdQdPtBasRect[j][nqdpt][qd], SdQdPtBasDxiRect[j][nqdpt][qd], SdQdPtBasDetRect[j][nqdpt][qd]);
			}
		} // end of qd
	}

	// rect vol qd
	for (unsigned int nqdpt = 1; nqdpt <= MAXSEGQDPTS; nqdpt++)
		for (unsigned int npet = 0; npet < nqdpt; npet++) // at each gauss points
			for (unsigned int npxi = 0; npxi < nqdpt; npxi++)
			{
				unsigned int pos = npet * nqdpt + npxi;
				QdPtCoeRect[nqdpt][pos][0] = QdPtCoe1d[nqdpt][npxi][0];
				QdPtCoeRect[nqdpt][pos][1] = QdPtCoe1d[nqdpt][npet][0];
				QdPtCoeRect[nqdpt][pos][3] = QdPtCoe1d[nqdpt][npxi][1] * QdPtCoe1d[nqdpt][npet][1];

				Get1dTensor(MAXFVORD,
					    QdPtBas1d[nqdpt][npxi], QdPtBasDe1d[nqdpt][npxi], QdPtBas1d[nqdpt][npet], QdPtBasDe1d[nqdpt][npet],
					    QdPtBasRect[nqdpt][pos], QdPtBasDxiRect[nqdpt][pos], QdPtBasDetRect[nqdpt][pos]);
			} // end of qd// end of nqdpt

	// rect vol qd 0: first node ...
	for (unsigned int npet = 0; npet < 2; npet++) // at each gauss points
		for (unsigned int npxi = 0; npxi < 2; npxi++)
		{
			unsigned int pos = npet * 2 + npxi;
			QdPtCoeRect[0][pos][0] = QdPtCoe1d[0][npxi][0];
			QdPtCoeRect[0][pos][1] = QdPtCoe1d[0][npet][0];
			QdPtCoeRect[0][pos][3] = QdPtCoe1d[0][npxi][1] * QdPtCoe1d[0][npet][1];

			Get1dTensor(MAXFVORD,
				    QdPtBas1d[0][npxi], QdPtBasDe1d[0][npxi], QdPtBas1d[0][npet], QdPtBasDe1d[0][npet],
				    QdPtBasRect[0][pos], QdPtBasDxiRect[0][pos], QdPtBasDetRect[0][pos]);
		} // end of qd// end of nqdpt
}

void Get1dTensor(
    unsigned int nord,
    double *Bas1d_Xi, double *Basde1d_Xi, double *Bas1d_Et, double *Basde1d_Et,
    double *BasRect, double *BasDxiRect, double *BasDetRect)
{
	unsigned int idx = 0;
	for (int i = 0; i <= (int)nord; i++)
		for (int k = i; k >= 0; k--)
		{
			int i_k = i - k;

			BasRect[idx] = Bas1d_Xi[k] * Bas1d_Et[i_k];
			BasDxiRect[idx] = Basde1d_Xi[k] * Bas1d_Et[i_k];
			BasDetRect[idx] = Bas1d_Xi[k] * Basde1d_Et[i_k];

			idx++;
		}
}

void Get1dTensorArbi(
    unsigned int nord, double xi, double et,
    double *BasRect, double *BasDxiRect, double *BasDetRect)
{
	double Bas1d_Xi[MAXFVORD + 1], Basde1d_Xi[MAXFVORD + 1];
	double Bas1d_Et[MAXFVORD + 1], Basde1d_Et[MAXFVORD + 1];

	for (unsigned int i = 0; i <= nord; i++)
	{
		Bas1d_Xi[i] = BasisPoly1d(i, xi);
		Bas1d_Et[i] = BasisPoly1d(i, et);
		Basde1d_Xi[i] = BasisPolyDe1d(i, xi);
		Basde1d_Et[i] = BasisPolyDe1d(i, et);
	}

	Get1dTensor(
	    nord, Bas1d_Xi, Basde1d_Xi, Bas1d_Et, Basde1d_Et, BasRect, BasDxiRect, BasDetRect);
}
