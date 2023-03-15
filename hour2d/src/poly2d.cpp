// poly 2d
// 2011 10 22 liming
// polynomials and its derivative
// 2d orthogonal basis
// for 2d poly order p means Np=(p+1)(p+2)/2 basis ( max p is 4,then Np=15, 0--14 )
// see the paper or the ppt for detail
#include "poly1d.h"
#include "poly2d.h"

////d2bs0
inline double d2bs0( double xi, double et )
{
	return 1.0;
}
//////////////////////////////////
////d2bs1
inline double d2bs1( double xi, double et )
{
	return 3.0*xi-1.0;
}
////d2bs2
inline double d2bs2( double xi, double et )
{
	return xi+2.0*et-1.0;
}
//////////////////////////////////////////
////d2bs3
inline double d2bs3( double xi, double et )
{
	return (10.0*xi-8.0)*xi+1.0;
}
////d2bs4
inline double d2bs4( double xi, double et )
{
	return (5.0*xi-1.0) * (xi+2.0*et-1.0);
}
////d2bs5
inline double d2bs5( double xi, double et )
{
	double xim1=xi-1.0;
//	return (xim1 + (3.0+1.7320508075688772935274463415059)*et) * 
//		   (xim1 + (3.0-1.7320508075688772935274463415059)*et);
	return (xim1 + 4.7320508075688772935274463415059*et) *
		   (xim1 + 1.2679491924311227064725536584941*et);
}
////////////////////////////////////////
////d2bs6
inline double d2bs6( double xi, double et )
{
	return ( (35.0*xi-45.0)*xi+15.0 )*xi-1.0;
}
////d2bs7
inline double d2bs7( double xi, double et )
{
	return ( (21.0*xi-12.0)*xi+1.0 ) * (xi+2.0*et-1.0);
}
////d2bs8
inline double d2bs8( double xi, double et )
{
	double xim1=xi-1.0;
	return (7.0*xi-1.0) * 
		(xim1 + 4.7320508075688772935274463415059*et) *
		(xim1 + 1.2679491924311227064725536584941*et);
}
////d2bs9
inline double d2bs9( double xi, double et )
{
	double xim1=xi-1.0;
	double et2=et*et;
	return ( (xim1+12.0*et)*xim1+30.0*et2 )*xim1 + 20.0*et*et2;
}
///////////////////////////////


////d2bs0
inline double d2bs0d10( double xi, double et )
{
	return 0.0;
}
//////////////////////////////////
////d2bs1
inline double d2bs1d10( double xi, double et )
{
	return 3.0;
}
////d2bs2
inline double d2bs2d10( double xi, double et )
{
	return 1.0;
}
//////////////////////////////////
////d2bs3
inline double d2bs3d10( double xi, double et )
{
	return 20.0*xi-8.0;
}
////d2bs4
inline double d2bs4d10( double xi, double et )
{
	return 10.0*(xi+et)-6.0;
}
////d2bs5
inline double d2bs5d10( double xi, double et )
{
	return 2.0*xi+6.0*et-2.0;
}
//////////////////////////////////
////d2bs6
inline double d2bs6d10( double xi, double et )
{
	return (105.0*xi-90.0)*xi+15.0;
}
////d2bs7
inline double d2bs7d10( double xi, double et )
{
	return (63.0*xi+84.0*et-66.0)*xi-24.0*et+13.0;
}
////d2bs8
inline double d2bs8d10( double xi, double et )
{
	double tmp = xi+et-0.57142857142857142857142857142857;
	return 42.0*tmp*tmp - (21.0*xi-18.0)*xi-
		4.7142857142857142857142857142857;
}
////d2bs9
inline double d2bs9d10( double xi, double et )
{
	double xim1=xi-1.0;
	return 3.0 * 
		   (xim1 + (4.0+2.4494897427831780981972840747059)*et) * 
		   (xim1 + (4.0-2.4494897427831780981972840747059)*et);
}
//////////////////////////////////////


////d2bs0
inline double d2bs0d11( double xi, double et )
{
	return 0.0;
}
//////////////////////////////////
////d2bs1
inline double d2bs1d11( double xi, double et )
{
	return 0.0;
}
////d2bs2
inline double d2bs2d11( double xi, double et )
{
	return 2.0;
}
//////////////////////////////////
////d2bs3
inline double d2bs3d11( double xi, double et )
{
	return 0.0;
}
////d2bs4
inline double d2bs4d11( double xi, double et )
{
	return 10.0*xi-2.0;
}
////d2bs5
inline double d2bs5d11( double xi, double et )
{
	return 6.0*(xi+2.0*et-1.0);
}
//////////////////////////////////
////d2bs6
inline double d2bs6d11( double xi, double et )
{
	return 0.0;
}
////d2bs7
inline double d2bs7d11( double xi, double et )
{
	return (42.0*xi-24.0)*xi+2.0;
}
////d2bs8
inline double d2bs8d11( double xi, double et )
{
	return (42.0*xi-6.0)*(xi+2.0*et-1.0);
}
////d2bs9
inline double d2bs9d11( double xi, double et )
{
	double xim1=xi-1.0;
	return 12.0 * 
		   (xim1 + (2.5+1.1180339887498948482045868343656)*et) * 
		   (xim1 + (2.5-1.1180339887498948482045868343656)*et);
}

////////////////////////////////
//dof which basis 
double BasisPolyTri( unsigned int dof, double xi, double et )
{
	switch( dof )
	{
		case 0:
			return d2bs0(xi,et);
            break;
		case 1:
            return d2bs1(xi,et)*1.4142135623730950488016887242097;
            break;
		case 2:
			return d2bs2(xi,et)*2.4494897427831780981972840747059;
			break;
		case 3:
			return d2bs3(xi,et)*1.7320508075688772935274463415059;
			break;
		case 4:
			return d2bs4(xi,et)*3.0;
			break;
		case 5:
			return d2bs5(xi,et)*3.8729833462074168851792653997824;
			break;
		case 6:
			return d2bs6(xi,et)*2.0;
			break;
		case 7:
			return d2bs7(xi,et)*3.4641016151377545870548926830117;
			break;
		case 8:
			return d2bs8(xi,et)*4.4721359549995793928183473374626;
            break;
		case 9:
            return d2bs9(xi,et)*5.2915026221291811810032315072785;
            break;
         default:
            return 0.0;
	}
}
////////////////////////////////
//dof which basis 
//twh (to which?) derivative to xi et?
double BasisPolyTriDe( unsigned int dof, double xi, double et, unsigned int twh )
{
	switch( twh )
	{
	case 0:
		switch( dof )
		{
		case 0:
			return d2bs0d10(xi,et);
			break;
		case 1:
			return d2bs1d10(xi,et)*1.4142135623730950488016887242097;
			break;
		case 2:
			return d2bs2d10(xi,et)*2.4494897427831780981972840747059;
			break;
		case 3:
			return d2bs3d10(xi,et)*1.7320508075688772935274463415059;
			break;
		case 4:
			return d2bs4d10(xi,et)*3.0;
			break;
		case 5:
			return d2bs5d10(xi,et)*3.8729833462074168851792653997824;
			break;
		case 6:
			return d2bs6d10(xi,et)*2.0;
			break;
		case 7:
			return d2bs7d10(xi,et)*3.4641016151377545870548926830117;
			break;
		case 8:
			return d2bs8d10(xi,et)*4.4721359549995793928183473374626;
            break;
		case 9:
            return d2bs9d10(xi,et)*5.2915026221291811810032315072785;
            break;
		default:
			return 0.0;
		}
		break;
	case 1:
		switch( dof )
		{
		case 0:
			return d2bs0d11(xi,et);
			break;
		case 1:
			return d2bs1d11(xi,et)*1.4142135623730950488016887242097;
			break;
		case 2:
			return d2bs2d11(xi,et)*2.4494897427831780981972840747059;
			break;
		case 3:
			return d2bs3d11(xi,et)*1.7320508075688772935274463415059;
			break;
		case 4:
			return d2bs4d11(xi,et)*3.0;
			break;
		case 5:
			return d2bs5d11(xi,et)*3.8729833462074168851792653997824;
			break;
		case 6:
			return d2bs6d11(xi,et)*2.0;
			break;
		case 7:
			return d2bs7d11(xi,et)*3.4641016151377545870548926830117;
			break;
		case 8:
			return d2bs8d11(xi,et)*4.4721359549995793928183473374626;
            break;
		case 9:
            return d2bs9d11(xi,et)*5.2915026221291811810032315072785;
            break;
		default:
			return 0.0;
		}
		break;
	default:
		return 0.0;
	}
}


//               use 0-MAXQDORDTRI :MAXQDPTSTRI
double QdPtBasTri   [MAXQDORDTRI+1][MAXQDPTSTRI][MAXFVDOF];
double QdPtBasDxiTri[MAXQDORDTRI+1][MAXQDPTSTRI][MAXFVDOF];
double QdPtBasDetTri[MAXQDORDTRI+1][MAXQDPTSTRI][MAXFVDOF];
//                    use 1-MAXSEGQDPTS :MAXSEGQDPTS
double SdQdPtTriXi    [6][MAXSEGQDPTS+1][MAXSEGQDPTS];
double SdQdPtTriEt    [6][MAXSEGQDPTS+1][MAXSEGQDPTS];
double SdQdPtBasTri   [6][MAXSEGQDPTS+1][MAXSEGQDPTS][MAXFVDOF];
double SdQdPtBasDxiTri[6][MAXSEGQDPTS+1][MAXSEGQDPTS][MAXFVDOF];
double SdQdPtBasDetTri[6][MAXSEGQDPTS+1][MAXSEGQDPTS][MAXFVDOF];

void CalBasArrTri()
{
	// tri side qd 0(first node for out sol) 1 2 ...
	for( unsigned int ts=0; ts<3; ts++ )
	{
		double xist = XiofTri[ts  ];
		double xied = XiofTri[ts+1];
		double etst = EtofTri[ts  ];
		double eted = EtofTri[ts+1];
		
		for( unsigned int nqdpt=0; nqdpt<=MAXSEGQDPTS; nqdpt++ )
		{
			unsigned int npnum = nqdpt;
			if( npnum==0 ) npnum = 1;
			// at each gauss points
			for( unsigned int np=0; np<npnum; np++ )
			{
				//double tst = GauQd1dCoe[nqdpt][np][0];
				//double ted = GauQd1dCoe[nqdpt][np][1];
				double tst = QdPtCoe1d3[nqdpt][np][0];
				double ted = QdPtCoe1d3[nqdpt][np][1];
				double xi0 = xied + (xist-xied) * tst;
				double et0 = eted + (etst-eted) * tst;
				SdQdPtTriXi[ts  ][nqdpt][np] = xi0;
				SdQdPtTriEt[ts  ][nqdpt][np] = et0;
				double xi1 = xied + (xist-xied) * ted;
				double et1 = eted + (etst-eted) * ted;
				SdQdPtTriXi[ts+3][nqdpt][np] = xi1;
				SdQdPtTriEt[ts+3][nqdpt][np] = et1;

				for( unsigned int i=0; i<MAXFVDOF; i++ )
				{
					SdQdPtBasTri   [ts  ][nqdpt][np][i] = BasisPolyTri  ( i, xi0, et0 );
					SdQdPtBasDxiTri[ts  ][nqdpt][np][i] = BasisPolyTriDe( i, xi0, et0, 0 );
					SdQdPtBasDetTri[ts  ][nqdpt][np][i] = BasisPolyTriDe( i, xi0, et0, 1 );
					SdQdPtBasTri   [ts+3][nqdpt][np][i] = BasisPolyTri  ( i, xi1, et1 );
					SdQdPtBasDxiTri[ts+3][nqdpt][np][i] = BasisPolyTriDe( i, xi1, et1, 0 );
					SdQdPtBasDetTri[ts+3][nqdpt][np][i] = BasisPolyTriDe( i, xi1, et1, 1 );
				}
			}// end of np
		}
	}// end of ts

	// tri vol qd 0: first node .. middle pt..
	for( unsigned int nqdord=0; nqdord<=MAXQDORDTRI; nqdord++ )
	{
		unsigned int nqdpt = SymQdTriDim[nqdord];
		if( nqdord==0 ) nqdpt = 6;
		
		// at each gauss points
		for( unsigned int np=0; np<nqdpt; np++ )
		{
			double xi = SymQdTriCoe[nqdord][np][0];
			double et = SymQdTriCoe[nqdord][np][1];

			for( unsigned int i=0; i<MAXFVDOF; i++ )
			{
				QdPtBasTri   [nqdord][np][i] = BasisPolyTri  ( i, xi, et );
				QdPtBasDxiTri[nqdord][np][i] = BasisPolyTriDe( i, xi, et, 0 );
				QdPtBasDetTri[nqdord][np][i] = BasisPolyTriDe( i, xi, et, 1 );
			}
		}// end of np
	}// end of nqdord
}
