#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <string>
#include <vector>
#include <algorithm>
#include "def2d.h"
#include "rcm.h"
using namespace std;

int adj_bandwidth(int node_num, int adj_num, int adj_row[], int adj[])

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
//
//  Discussion:
//
//    Thanks to Man Yuan, of Southeast University, China, for pointing out
//    an inconsistency in the indexing of ADJ_ROW, 06 June 2011.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 June 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Output, int ADJ_BANDWIDTH, the bandwidth of the adjacency
//    matrix.
//
{
	int band_hi;
	int band_lo;
	int col;
	int i;
	int j;
	int value;

	band_lo = 0;
	band_hi = 0;

	for (i = 0; i < node_num; i++)
	{
		for (j = adj_row[i]; j <= adj_row[i + 1] - 1; j++)
		{
			col = adj[j - 1] - 1;
			band_lo = i4_max(band_lo, i - col);
			band_hi = i4_max(band_hi, col - i);
		}
	}

	value = band_lo + 1 + band_hi;

	return value;
}

int adj_perm_bandwidth(int node_num, int adj_num, int adj_row[], int adj[],
		       int perm[], int perm_inv[])

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
//
//  Discussion:
//
//    The matrix is defined by the adjacency information and a permutation.
//
//    The routine also computes the bandwidth and the size of the envelope.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, int PERM[NODE_NUM], PERM_INV(NODE_NUM), the permutation
//    and inverse permutation.
//
//    Output, int ADJ_PERM_BANDWIDTH, the bandwidth of the permuted
//    adjacency matrix.
//
{
	int band_hi;
	int band_lo;
	int bandwidth;
	int col;
	int i;
	int j;

	band_lo = 0;
	band_hi = 0;

	for (i = 0; i < node_num; i++)
	{
		for (j = adj_row[perm[i] - 1]; j <= adj_row[perm[i]] - 1; j++)
		{
			col = perm_inv[adj[j - 1] - 1];
			band_lo = i4_max(band_lo, i - col);
			band_hi = i4_max(band_hi, col - i);
		}
	}

	bandwidth = band_lo + 1 + band_hi;

	return bandwidth;
}
//****************************************************************************80

void adj_perm_show(int node_num, int adj_num, int adj_row[], int adj[],
		   int perm[], int perm_inv[])

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_PERM_SHOW displays a symbolic picture of a permuted adjacency matrix.
//
//  Discussion:
//
//    The matrix is defined by the adjacency information and a permutation.
//
//    The routine also computes the bandwidth and the size of the envelope.
//
//    If no permutation has been done, you must set PERM(I) = PERM_INV(I) = I
//    before calling this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, int PERM[NODE_NUM], PERM_INV[NODE_NUM], the permutation
//    and inverse permutation.
//
{
	char *band;
	int band_lo;
	int col;
	int i;
	int j;
	int k;
	int nonzero_num;

	band = new char[node_num];

	band_lo = 0;
	nonzero_num = 0;

	cout << "\n";
	cout << "  Nonzero structure of matrix:\n";
	cout << "\n";

	for (i = 0; i < node_num; i++)
	{
		for (k = 0; k < node_num; k++)
		{
			band[k] = '.';
		}

		band[i] = 'D';

		for (j = adj_row[perm[i] - 1]; j <= adj_row[perm[i]] - 1; j++)
		{
			col = perm_inv[adj[j - 1] - 1] - 1;

			if (col < i)
			{
				nonzero_num = nonzero_num + 1;
			}

			band_lo = i4_max(band_lo, i - col);

			if (col != i)
			{
				band[col] = 'X';
			}
		}
		cout << "  " << setw(8) << i + 1 << " ";
		for (j = 0; j < node_num; j++)
		{
			cout << band[j];
		}
		cout << "\n";
	}
	cout << "\n";
	cout << "  Lower bandwidth = " << band_lo << "\n";
	cout << "  Lower envelope contains " << nonzero_num << " nonzeros.\n";

	return;
}

//****************************************************************************80

void adj_show(int node_num, int adj_num, int adj_row[], int adj[])

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_SHOW displays a symbolic picture of an adjacency matrix.
//
//  Discussion:
//
//    The matrix is defined by the adjacency information and a permutation.
//
//    The routine also computes the bandwidth and the size of the envelope.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
{
	char *band;
	int band_lo;
	int col;
	int i;
	int j;
	int k;
	int nonzero_num;

	band = new char[node_num];

	band_lo = 0;
	nonzero_num = 0;

	cout << "\n";
	cout << "  Nonzero structure of matrix:\n";
	cout << "\n";

	for (i = 0; i < node_num; i++)
	{
		for (k = 0; k < node_num; k++)
		{
			band[k] = '.';
		}

		band[i] = 'D';

		for (j = adj_row[i]; j <= adj_row[i + 1] - 1; j++)
		{
			col = adj[j - 1] - 1;
			if (col < i)
			{
				nonzero_num = nonzero_num + 1;
			}
			band_lo = max(band_lo, i - col);
			band[col] = 'X';
		}
		cout << "  " << setw(8) << i + 1 << " ";
		for (j = 0; j < node_num; j++)
		{
			cout << band[j];
		}
		cout << "\n";
	}

	cout << "\n";
	cout << "  Lower bandwidth = " << band_lo << "\n";
	cout << "  Lower envelope contains " << nonzero_num << " nonzeros.\n";

	delete[] band;

	return;
}
//****************************************************************************80

void degree(int root, int adj_num, int adj_row[], int adj[], int mask[],
	    int deg[], int *iccsze, int ls[], int node_num)

//****************************************************************************80
//
//  Purpose:
//
//    DEGREE computes the degrees of the nodes in the connected component.
//
//  Discussion:
//
//    The connected component is specified by MASK and ROOT.
//    Nodes for which MASK is zero are ignored.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int ROOT, the node that defines the connected component.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, int MASK[NODE_NUM], is nonzero for those nodes which are
//    to be considered.
//
//    Output, int DEG[NODE_NUM], contains, for each  node in the connected
//    component, its degree.
//
//    Output, int *ICCSIZE, the number of nodes in the connected component.
//
//    Output, int LS[NODE_NUM], stores in entries 1 through ICCSIZE the nodes
//    in the connected component, starting with ROOT, and proceeding
//    by levels.
//
//    Input, int NODE_NUM, the number of nodes.
//
{
	int i;
	int ideg;
	int j;
	int jstop;
	int jstrt;
	int lbegin;
	int lvlend;
	int lvsize;
	int nbr;
	int node;
	//
	//  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
	//
	ls[0] = root;
	adj_row[root - 1] = -adj_row[root - 1];
	lvlend = 0;
	*iccsze = 1;
	//
	//  LBEGIN is the pointer to the beginning of the current level, and
	//  LVLEND points to the end of this level.
	//
	for (;;)
	{
		lbegin = lvlend + 1;
		lvlend = *iccsze;
		//
		//  Find the degrees of nodes in the current level,
		//  and at the same time, generate the next level.
		//
		for (i = lbegin; i <= lvlend; i++)
		{
			node = ls[i - 1];
			jstrt = -adj_row[node - 1];
			jstop = abs(adj_row[node]) - 1;
			ideg = 0;

			for (j = jstrt; j <= jstop; j++)
			{
				nbr = adj[j - 1];

				if (mask[nbr - 1] != 0)
				{
					ideg = ideg + 1;

					if (0 <= adj_row[nbr - 1])
					{
						adj_row[nbr - 1] = -adj_row[nbr - 1];
						*iccsze = *iccsze + 1;
						ls[*iccsze - 1] = nbr;
					}
				}
			}
			deg[node - 1] = ideg;
		}
		//
		//  Compute the current level width.
		//
		lvsize = *iccsze - lvlend;
		//
		//  If the current level width is nonzero, generate another level.
		//
		if (lvsize == 0)
		{
			break;
		}
	}
	//
	//  Reset ADJ_ROW to its correct sign and return.
	//
	for (i = 0; i < *iccsze; i++)
	{
		node = ls[i] - 1;
		adj_row[node] = -adj_row[node];
	}

	return;
}
//****************************************************************************80

void genrcm(int node_num, int adj_num, int *adj_row, int *adj, int *perm)

//****************************************************************************80
//
//  Purpose:
//
//    GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
//
//  Discussion:
//
//    For each connected component in the graph, the routine obtains
//    an ordering by calling RCM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int  ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Output, int  PERM[NODE_NUM], the RCM ordering.
//
//  Local Parameters:
//
//    Local, int  LEVEL_ROW[NODE_NUM+1], the index vector for a level
//    structure.  The level structure is stored in the currently unused
//    spaces in the permutation vector PERM.
//
//    Local, int MASK[NODE_NUM], marks variables that have been numbered.
//
{
	int i;
	int iccsze;
	int level_num;
	int *level_row;
	int *mask;
	int num;
	int root;

	level_row = new int[node_num + 1];
	mask = new int[node_num];

	for (i = 0; i < node_num; i++)
	{
		mask[i] = 1;
	}

	num = 1;

	for (i = 0; i < node_num; i++)
	{
		//
		//  For each masked connected component...
		//
		if (mask[i] != 0)
		{
			root = i + 1;
			//
			//  Find a pseudo-peripheral node ROOT.  The level structure found by
			//  ROOT_FIND is stored starting at PERM(NUM).
			//
			root_find(&root, adj_num, adj_row, adj, mask, &level_num,
				  level_row, perm + num - 1, node_num);
			//
			//  RCM orders the component using ROOT as the starting node.
			//
			rcm(root, adj_num, adj_row, adj, mask, perm + num - 1, &iccsze,
			    node_num);

			num = num + iccsze;
			//
			//  We can stop once every node is in one of the connected components.
			//
			if (node_num < num)
			{
				free(level_row);
				free(mask);
				return;
			}
		}
	}

	free(level_row);
	free(mask);

	return;
}
//****************************************************************************80

int i4_max(int i1, int i2)

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
	if (i2 < i1)
	{
		return i1;
	}
	else
	{
		return i2;
	}
}
//****************************************************************************80

void i4vec_reverse(int n, int a[])

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_REVERSE reverses the elements of an I4VEC.
//
//  Example:
//
//    Input:
//
//      N = 5,
//      A = ( 11, 12, 13, 14, 15 ).
//
//    Output:
//
//      A = ( 15, 14, 13, 12, 11 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A(N), the array to be reversed.
//
{
	int i;
	int j;

	for (i = 0; i < n / 2; i++)
	{
		j = a[i];
		a[i] = a[n - 1 - i];
		a[n - 1 - i] = j;
	}

	return;
}
//****************************************************************************80

void level_set(int root, int adj_num, int adj_row[], int adj[], int mask[],
	       int *level_num, int level_row[], int level[], int node_num)

//****************************************************************************80
//
//  Purpose:
//
//    LEVEL_SET generates the connected level structure rooted at a given node.
//
//  Discussion:
//
//    Only nodes for which MASK is nonzero will be considered.
//
//    The root node chosen by the user is assigned level 1, and masked.
//    All (unmasked) nodes reachable from a node in level 1 are
//    assigned level 2 and masked.  The process continues until there
//    are no unmasked nodes adjacent to any node in the current level.
//    The number of levels may vary between 2 and NODE_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int ROOT, the node at which the level structure
//    is to be rooted.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input/output, int MASK[NODE_NUM].  On input, only nodes with nonzero
//    MASK are to be processed.  On output, those nodes which were included
//    in the level set have MASK set to 1.
//
//    Output, int *LEVEL_NUM, the number of levels in the level
//    structure.  ROOT is in level 1.  The neighbors of ROOT
//    are in level 2, and so on.
//
//    Output, int LEVEL_ROW[NODE_NUM+1], LEVEL[NODE_NUM], the rooted
//    level structure.
//
//    Input, int NODE_NUM, the number of nodes.
//
{
	int i;
	int iccsze;
	int j;
	int jstop;
	int jstrt;
	int lbegin;
	int lvlend;
	int lvsize;
	int nbr;
	int node;

	mask[root - 1] = 0;
	level[0] = root;
	*level_num = 0;
	lvlend = 0;
	iccsze = 1;
	//
	//  LBEGIN is the pointer to the beginning of the current level, and
	//  LVLEND points to the end of this level.
	//
	for (;;)
	{
		lbegin = lvlend + 1;
		lvlend = iccsze;
		*level_num = *level_num + 1;
		level_row[*level_num - 1] = lbegin;
		//
		//  Generate the next level by finding all the masked neighbors of nodes
		//  in the current level.
		//
		for (i = lbegin; i <= lvlend; i++)
		{
			node = level[i - 1];
			jstrt = adj_row[node - 1];
			jstop = adj_row[node] - 1;

			for (j = jstrt; j <= jstop; j++)
			{
				nbr = adj[j - 1];

				if (mask[nbr - 1] != 0)
				{
					iccsze = iccsze + 1;
					level[iccsze - 1] = nbr;
					mask[nbr - 1] = 0;
				}
			}
		}
		//
		//  Compute the current level width (the number of nodes encountered.)
		//  If it is positive, generate the next level.
		//
		lvsize = iccsze - lvlend;

		if (lvsize <= 0)
		{
			break;
		}
	}

	level_row[*level_num] = lvlend + 1;
	//
	//  Reset MASK to 1 for the nodes in the level structure.
	//
	for (i = 0; i < iccsze; i++)
	{
		mask[level[i] - 1] = 1;
	}

	return;
}
//****************************************************************************80

void perm_inverse3(int n, int perm[], int perm_inv[])

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE3 produces the inverse of a given permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items permuted.
//
//    Input, int PERM[N], a permutation.
//
//    Output, int PERM_INV[N], the inverse permutation.
//
{
	int i;

	for (i = 0; i < n; i++)
	{
		perm_inv[perm[i] - 1] = i;
	}

	return;
}
//****************************************************************************80

void rcm(int root, int adj_num, int adj_row[], int adj[], int mask[],
	 int perm[], int *iccsze, int node_num)

//****************************************************************************80
//
//  Purpose:
//
//    RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
//
//  Discussion:
//
//    The connected component is specified by a node ROOT and a mask.
//    The numbering starts at the root node.
//
//    An outline of the algorithm is as follows:
//
//    X(1) = ROOT.
//
//    for ( I = 1 to N-1)
//      Find all unlabeled neighbors of X(I),
//      assign them the next available labels, in order of increasing degree.
//
//    When done, reverse the ordering.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 August 2013
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//  Parameters:
//
//    Input, int ROOT, the node that defines the connected component.
//    It is used as the starting point for the RCM ordering.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input/output, int MASK[NODE_NUM], a mask for the nodes.  Only
//    those nodes with nonzero input mask values are considered by the
//    routine.  The nodes numbered by RCM will have their mask values
//    set to zero.
//
//    Output, int PERM[NODE_NUM], the RCM ordering.
//
//    Output, int *ICCSZE, the size of the connected component
//    that has been numbered.
//
//    Input, int NODE_NUM, the number of nodes.
//
//  Local Parameters:
//
//    Workspace, int DEG[NODE_NUM], a temporary vector used to hold
//    the degree of the nodes in the section graph specified by mask and root.
//
{
	int *deg;
	int fnbr;
	int i;
	int j;
	int jstop;
	int jstrt;
	int k;
	int l;
	int lbegin;
	int lnbr;
	int lperm;
	int lvlend;
	int nbr;
	int node;
	//
	//  If node_num out of bounds, something is wrong.
	//
	if (node_num < 1)
	{
		cerr << "\n";
		cerr << "RCM - Fatal error!\n";
		cerr << "  Unacceptable input value of NODE_NUM = " << node_num << "\n";
		exit(1);
	}
	//
	//  If the root is out of bounds, something is wrong.
	//
	if (root < 1 || node_num < root)
	{
		cerr << "\n";
		cerr << "RCM - Fatal error!\n";
		cerr << "  Unacceptable input value of ROOT = " << root << "\n";
		cerr << "  Acceptable values are between 1 and " << node_num << ", inclusive.\n";
		exit(1);
	}
	//
	//  Allocate memory for the degree array.
	//
	deg = new int[node_num];
	//
	//  Find the degrees of the nodes in the component specified by MASK and ROOT.
	//
	degree(root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num);
	//
	//  If the connected component size is less than 1, something is wrong.
	//
	if (*iccsze < 1)
	{
		cerr << "\n";
		cerr << "RCM - Fatal error!\n";
		cerr << "  Connected component size ICCSZE returned from DEGREE as "
		     << *iccsze << "\n";
		exit(1);
	}
	//
	//  Set the mask value for the root.
	//
	mask[root - 1] = 0;
	//
	//  If the connected component is a singleton, there is no ordering necessary.
	//
	if (*iccsze == 1)
	{
		delete[] deg;
		return;
	}
	//
	//  Carry out the reordering.
	//
	//  LBEGIN and LVLEND point to the beginning and
	//  the end of the current level respectively.
	//
	lvlend = 0;
	lnbr = 1;

	while (lvlend < lnbr)
	{
		lbegin = lvlend + 1;
		lvlend = lnbr;

		for (i = lbegin; i <= lvlend; i++)
		{
			//
			//  For each node in the current level...
			//
			node = perm[i - 1];
			jstrt = adj_row[node - 1];
			jstop = adj_row[node] - 1;
			//
			//  Find the unnumbered neighbors of NODE.
			//
			//  FNBR and LNBR point to the first and last neighbors
			//  of the current node in PERM.
			//
			fnbr = lnbr + 1;

			for (j = jstrt; j <= jstop; j++)
			{
				nbr = adj[j - 1];

				if (mask[nbr - 1] != 0)
				{
					lnbr = lnbr + 1;
					mask[nbr - 1] = 0;
					perm[lnbr - 1] = nbr;
				}
			}
			//
			//  If no neighbors, skip to next node in this level.
			//
			if (lnbr <= fnbr)
			{
				continue;
			}
			//
			//  Sort the neighbors of NODE in increasing order by degree.
			//  Linear insertion is used.
			//
			k = fnbr;

			while (k < lnbr)
			{
				l = k;
				k = k + 1;
				nbr = perm[k - 1];

				while (fnbr < l)
				{
					lperm = perm[l - 1];

					if (deg[lperm - 1] <= deg[nbr - 1])
					{
						break;
					}

					perm[l] = lperm;
					l = l - 1;
				}
				perm[l] = nbr;
			}
		}
	}
	//
	//  We now have the Cuthill-McKee ordering.
	//  Reverse it to get the Reverse Cuthill-McKee ordering.
	//
	i4vec_reverse(*iccsze, perm);
	//
	//  Free memory.
	//
	delete[] deg;

	return;
}
//****************************************************************************80

void root_find(int *root, int adj_num, int adj_row[], int adj[], int mask[],
	       int *level_num, int level_row[], int level[], int node_num)

//****************************************************************************80
//
//  Purpose:
//
//    ROOT_FIND finds a pseudo-peripheral node.
//
//  Discussion:
//
//    The diameter of a graph is the maximum distance (number of edges)
//    between any two nodes of the graph.
//
//    The eccentricity of a node is the maximum distance between that
//    node and any other node of the graph.
//
//    A peripheral node is a node whose eccentricity equals the
//    diameter of the graph.
//
//    A pseudo-peripheral node is an approximation to a peripheral node;
//    it may be a peripheral node, but all we know is that we tried our
//    best.
//
//    The routine is given a graph, and seeks pseudo-peripheral nodes,
//    using a modified version of the scheme of Gibbs, Poole and
//    Stockmeyer.  It determines such a node for the section subgraph
//    specified by MASK and ROOT.
//
//    The routine also determines the level structure associated with
//    the given pseudo-peripheral node; that is, how far each node
//    is from the pseudo-peripheral node.  The level structure is
//    returned as a list of nodes LS, and pointers to the beginning
//    of the list of nodes that are at a distance of 0, 1, 2, ...,
//    NODE_NUM-1 from the pseudo-peripheral node.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2007
//
//  Author:
//
//    Original FORTRAN77 version by Alan George, Joseph Liu.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Alan George, Joseph Liu,
//    Computer Solution of Large Sparse Positive Definite Systems,
//    Prentice Hall, 1981.
//
//    Norman Gibbs, William Poole, Paul Stockmeyer,
//    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
//    SIAM Journal on Numerical Analysis,
//    Volume 13, pages 236-250, 1976.
//
//    Norman Gibbs,
//    Algorithm 509: A Hybrid Profile Reduction Algorithm,
//    ACM Transactions on Mathematical Software,
//    Volume 2, pages 378-387, 1976.
//
//  Parameters:
//
//    Input/output, int *ROOT.  On input, ROOT is a node in the
//    the component of the graph for which a pseudo-peripheral node is
//    sought.  On output, ROOT is the pseudo-peripheral node obtained.
//
//    Input, int ADJ_NUM, the number of adjacency entries.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
//    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
//
//    Input, int ADJ[ADJ_NUM], the adjacency structure.
//    For each row, it contains the column indices of the nonzero entries.
//
//    Input, int MASK[NODE_NUM], specifies a section subgraph.  Nodes
//    for which MASK is zero are ignored by FNROOT.
//
//    Output, int *LEVEL_NUM, is the number of levels in the level structure
//    rooted at the node ROOT.
//
//    Output, int LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the
//    level structure array pair containing the level structure found.
//
//    Input, int NODE_NUM, the number of nodes.
//
{
	int iccsze;
	int j;
	int jstrt;
	int k;
	int kstop;
	int kstrt;
	int level_num2;
	int mindeg;
	int nabor;
	int ndeg;
	int node;
	//
	//  Determine the level structure rooted at ROOT.
	//
	level_set(*root, adj_num, adj_row, adj, mask, level_num,
		  level_row, level, node_num);
	//
	//  Count the number of nodes in this level structure.
	//
	iccsze = level_row[*level_num] - 1;
	//
	//  Extreme case:
	//    A complete graph has a level set of only a single level.
	//    Every node is equally good (or bad).
	//
	if (*level_num == 1)
	{
		return;
	}
	//
	//  Extreme case:
	//    A "line graph" 0--0--0--0--0 has every node in its only level.
	//    By chance, we've stumbled on the ideal root.
	//
	if (*level_num == iccsze)
	{
		return;
	}
	//
	//  Pick any node from the last level that has minimum degree
	//  as the starting point to generate a new level set.
	//
	for (;;)
	{
		mindeg = iccsze;

		jstrt = level_row[*level_num - 1];
		*root = level[jstrt - 1];

		if (jstrt < iccsze)
		{
			for (j = jstrt; j <= iccsze; j++)
			{
				node = level[j - 1];
				ndeg = 0;
				kstrt = adj_row[node - 1];
				kstop = adj_row[node] - 1;

				for (k = kstrt; k <= kstop; k++)
				{
					nabor = adj[k - 1];
					if (0 < mask[nabor - 1])
					{
						ndeg = ndeg + 1;
					}
				}

				if (ndeg < mindeg)
				{
					*root = node;
					mindeg = ndeg;
				}
			}
		}
		//
		//  Generate the rooted level structure associated with this node.
		//
		level_set(*root, adj_num, adj_row, adj, mask, &level_num2,
			  level_row, level, node_num);
		//
		//  If the number of levels did not increase, accept the new ROOT.
		//
		if (level_num2 <= *level_num)
		{
			break;
		}

		*level_num = level_num2;
		//
		//  In the unlikely case that ROOT is one endpoint of a line graph,
		//  we can exit now.
		//
		if (iccsze <= *level_num)
		{
			break;
		}
	}

	return;
}

//****************************************************************************
// extern running function to hour2d code
typedef struct elemsort
{
	unsigned int Node[4];
	unsigned int Side[4];
	double xyc[2];
	double volr;
	int number;
	int idx;
} esort;

typedef struct idxpele
{
	unsigned int Node[2];
	int ofElem[2];		// 1: left elem; 0: right element  包含边界条件 0 < 1
	unsigned int lsdidx[2]; // local side index of 1: left elem; 0: right elem
	int number;
	int idx;
} iele;

void GridRenumber_node(unsigned int nele, unsigned int nnod, elem *pele, node *pnod)
{

	int iftera = 0; // 判断三角形或四边形
	int adj_num, *adj, *adj_row, bandwidth;
	int *perm, *perm_inv, node_num;

	if (pele[0].Node[3] == pele[0].Node[0])
	{ //三角形
		iftera = 3;
	}
	else
	{ //四边形
		iftera = 4;
	}

	node_num = nnod;
	adj_row = (int *)malloc((node_num + 1) * sizeof(int));
	perm = (int *)malloc(node_num * sizeof(int));
	perm_inv = (int *)malloc(node_num * sizeof(int));

	//网格点格式转化
	std::vector<int> adjvec;
	adjvec = GridToRCM_node(node_num, adj_row, iftera, nele, pele);
	adj_num = adjvec.size();
	adj = (int *)calloc(adj_num, sizeof(int));
	for (unsigned int i = 0; i < adj_num; i++)
	{
		adj[i] = adjvec[i];
	}

	// RCM算法
	genrcm(node_num, adj_num, adj_row, adj, perm);
	// 重排序
	perm_inverse3(node_num, perm, perm_inv);
	// adj_show( node_num, adj_num, adj_row, adj);
	for (unsigned int i = 0; i < nele; i++)
	{
		// printf("%d  ",i+1);
		for (unsigned int j = 0; j < iftera; j++)
		{
			unsigned int indexN = pele[i].Node[j];
			pele[i].Node[j] = perm_inv[indexN];
			// printf( "%d  ",pele[i].Node[j] + 1 );
		}
		// printf("\n");
	}
	if (iftera == 3)
	{
		for (unsigned int i = 0; i < nele; i++)
		{
			pele[i].Node[3] = pele[i].Node[0];
		}
	}

	// adj_perm_show( node_num, adj_num, adj_row, adj, perm, perm_inv );

	int origin_bandwidth = adj_bandwidth(node_num, adj_num, adj_row, adj);
	int reversed_bandwidth = adj_perm_bandwidth(node_num, adj_num, adj_row, adj, perm, perm_inv);
	printf("   Original bandwidth:%d , Renumbered bandwidth:%d\n", origin_bandwidth, reversed_bandwidth);

	double *xx = (double *)malloc(sizeof(double) * nnod);
	double *yy = (double *)malloc(sizeof(double) * nnod);

	for (unsigned int i = 0; i < nnod; i++)
	{
		xx[i] = pnod[i].xy[0];
		yy[i] = pnod[i].xy[1];
	}
	for (unsigned int i = 0; i < nnod; i++)
	{
		pnod[perm_inv[i]].xy[0] = xx[i];
		pnod[perm_inv[i]].xy[1] = yy[i];
		// printf("%d %.9e %.9e\n", i+1, pnod[i].xy[0], pnod[i].xy[1]);
	}

	std::vector<int> temp = adjvec;
	adjvec.swap(temp);
	free(adj_row);
	free(adj);
	free(perm);
	free(perm_inv);
	free(xx);
	free(yy);
	return;
}
vector<int> GridToRCM_node(int node_num, int *adj_row, int iftera, unsigned int nele, elem *pele)
{

	adj_row[0] = 1;
	int adj_num = 0;
	std::vector<int> adjvec;

	for (int i = 0; i < node_num; i++)
	{
		std::vector<int> temp;
		int count = 0;
		if (iftera == 3)
		{
			for (int j = 0; j < nele; j++)
			{
				for (int k = 0; k < iftera; k++)
				{
					if (i == pele[j].Node[k])
					{
						temp.insert(temp.end(), pele[j].Node[0] + 1);
						temp.insert(temp.end(), pele[j].Node[1] + 1);
						temp.insert(temp.end(), pele[j].Node[2] + 1);
						count += iftera;
						break;
					}
				}
			}
		}
		else if (iftera == 4)
		{
			for (int j = 0; j < nele; j++)
			{
				for (int k = 0; k < iftera; k++)
				{
					if (i == pele[j].Node[k])
					{
						switch (k)
						{
						case 0:
							temp.insert(temp.end(), pele[j].Node[0] + 1);
							temp.insert(temp.end(), pele[j].Node[1] + 1);
							temp.insert(temp.end(), pele[j].Node[3] + 1);
							break;
						case 1:
							temp.insert(temp.end(), pele[j].Node[0] + 1);
							temp.insert(temp.end(), pele[j].Node[1] + 1);
							temp.insert(temp.end(), pele[j].Node[2] + 1);
							break;
						case 2:
							temp.insert(temp.end(), pele[j].Node[1] + 1);
							temp.insert(temp.end(), pele[j].Node[2] + 1);
							temp.insert(temp.end(), pele[j].Node[3] + 1);
							break;
						case 3:
							temp.insert(temp.end(), pele[j].Node[2] + 1);
							temp.insert(temp.end(), pele[j].Node[3] + 1);
							temp.insert(temp.end(), pele[j].Node[0] + 1);
						}

						count += iftera - 1;
						break;
					}
				}
			}
		}

		//元素去重
		for (int ii = 0; ii < count - 1; ii++)
		{
			for (int jj = ii + 1; jj < count; jj++)
			{
				if ((temp[ii] == temp[jj]))
				{
					vector<int>::iterator ite = temp.begin();
					ite = temp.erase(ite + jj);
					jj--;
					count--;
				}
			}
		}
		for (int ii = 0; ii < count; ii++)
		{
			if (temp[ii] == i + 1)
			{
				vector<int>::iterator ite = temp.begin();
				ite = temp.erase(ite + ii);
				ii--;
				count--;
			}
		}
		//元素升序排序
		sort(temp.begin(), temp.end());

		for (int ii = 0; ii < count; ii++)
		{
			adjvec.insert(adjvec.end(), temp[ii]);
		}

		adj_row[i + 1] = adj_row[i] + count;
		std::vector<int> tmp1 = temp; // free the vector
		temp.swap(tmp1);
	}
	// adj_num = adjvec.size();
	return adjvec;
}

int *GridRenumber_elem(unsigned int nele, unsigned int nnod, unsigned int nsid, elem *pele, node *pnod, side *psid)
{

	int iftera = 0; // 判断三角形或四边形
	int adj_num, *adj, *adj_row, bandwidth;
	int *perm, *perm_inv, node_num;

	if (pele[0].Node[3] == pele[0].Node[0])
	{ //三角形
		iftera = 3;
	}
	else
	{ //四边形
		iftera = 4;
	}

	node_num = nele;
	adj_row = (int *)malloc((node_num + 1) * sizeof(int));
	perm = (int *)malloc(node_num * sizeof(int));
	perm_inv = (int *)malloc(node_num * sizeof(int));

	//网格点格式转化
	std::vector<int> adjvec;
	adjvec = GridToRCM_elem(node_num, adj_row, iftera, nele, pele);
	adj_num = adjvec.size();
	adj = (int *)calloc(adj_num, sizeof(int));
	for (unsigned int i = 0; i < adj_num; i++)
	{
		adj[i] = adjvec[i];
	}
	std::vector<int> temp = adjvec;
	adjvec.swap(temp);
	// adj_show( node_num, adj_num, adj_row, adj);

	// RCM算法
	genrcm(node_num, adj_num, adj_row, adj, perm);
	perm_inverse3(node_num, perm, perm_inv);

	// adj_perm_show( node_num, adj_num, adj_row, adj, perm, perm_inv );
	int origin_bandwidth = adj_bandwidth(node_num, adj_num, adj_row, adj);
	int reversed_bandwidth = adj_perm_bandwidth(node_num, adj_num, adj_row, adj, perm, perm_inv);
	printf("   Original bandwidth:%d , Renumbered bandwidth:%d\n", origin_bandwidth, reversed_bandwidth);

	free(adj_row);
	free(adj);
	free(perm);
	// free(perm_inv);

	return perm_inv;
}
vector<int> GridToRCM_elem(int node_num, int *adj_row, int iftera, unsigned int nele, elem *pele)
{
	adj_row[0] = 1;
	int adj_num = 0;
	std::vector<int> adjvec;

	for (int i = 0; i < node_num; i++)
	{
		int count = 1;
		std::vector<int> temp;
		temp.insert(temp.end(), i + 1);
		for (int j = 0; j < iftera; j++)
		{
			int elemflag = 0;
			int tp1 = pele[i].Node[j];
			int tp2 = pele[i].Node[(j + 1) % iftera];
			for (int k = 0; k < nele; k++)
			{
				if (i == k)
				{
					continue;
				}
				for (unsigned int l = 0; l < iftera; l++)
				{
					if ((tp1 == pele[k].Node[l]) || (tp2 == pele[k].Node[l]))
					{
						elemflag++;
					}
				}
				if (elemflag == 2)
				{
					count++;
					temp.insert(temp.end(), k + 1);
					break;
				}
				elemflag = 0;
			}
		}

		//元素升序排序
		sort(temp.begin(), temp.end());

		// printf("%d  ",i+1);
		for (int ii = 0; ii < count; ii++)
		{
			adjvec.insert(adjvec.end(), temp[ii]);
			// printf("%d  ",temp[ii]);
		}
		// printf("\n");
		adj_row[i + 1] = adj_row[i] + count;
		std::vector<int> tmp1 = temp; // free the vector
		temp.swap(tmp1);
	}
	// adj_num = adjvec.size();
	return adjvec;
}

int *GridRenumber_side(unsigned int nele, unsigned int nnod, unsigned int nsid, elem *pele, node *pnod, side *psid)
{

	int iftera = 0; // 判断三角形或四边形
	int adj_num, *adj, *adj_row, bandwidth;
	int *perm, *perm_inv, node_num;

	if (pele[0].Node[3] == pele[0].Node[0])
	{ //三角形
		iftera = 3;
	}
	else
	{ //四边形
		iftera = 4;
	}

	node_num = nsid;
	adj_row = (int *)malloc((node_num + 1) * sizeof(int));
	perm = (int *)malloc(node_num * sizeof(int));
	perm_inv = (int *)malloc(node_num * sizeof(int));

	//网格点格式转化
	std::vector<int> adjvec;
	adjvec = GridToRCM_elem(node_num, adj_row, iftera, nele, pele);
	adj_num = adjvec.size();
	adj = (int *)calloc(adj_num, sizeof(int));
	for (unsigned int i = 0; i < adj_num; i++)
	{
		adj[i] = adjvec[i];
	}
	std::vector<int> temp = adjvec;
	adjvec.swap(temp);
	// adj_show( node_num, adj_num, adj_row, adj);

	// RCM算法
	genrcm(node_num, adj_num, adj_row, adj, perm);
	perm_inverse3(node_num, perm, perm_inv);

	// adj_perm_show( node_num, adj_num, adj_row, adj, perm, perm_inv );
	int origin_bandwidth = adj_bandwidth(node_num, adj_num, adj_row, adj);
	int reversed_bandwidth = adj_perm_bandwidth(node_num, adj_num, adj_row, adj, perm, perm_inv);
	printf("   Original bandwidth:%d , Renumbered bandwidth:%d\n", origin_bandwidth, reversed_bandwidth);

	std::vector<iele> ievec;
	ievec.reserve(node_num);
	for (unsigned int i = 0; i < node_num; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			ievec[i].Node[j] = psid[i].Node[j];
			ievec[i].ofElem[j] = psid[i].ofElem[j];
			ievec[i].lsdidx[j] = psid[i].lsdidx[j];
		}
	}
	for (unsigned int i = 0; i < node_num; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			psid[perm_inv[i]].Node[j] = ievec[i].Node[j];
			psid[perm_inv[i]].ofElem[j] = ievec[i].ofElem[j];
			psid[perm_inv[i]].lsdidx[j] = ievec[i].lsdidx[j];
		}
	}
	for (unsigned int i = 0; i < nele; i++)
	{
		for (unsigned int j = 0; j < iftera; j++)
		{
			pele[i].Side[j] = perm_inv[pele[i].Side[j]];
		}
	}
	for (unsigned int i = 0; i < nnod; i++)
	{
		for (unsigned int j = 0; j < pnod[i].nNElem; j++)
		{
			pnod[i].ofSide[j] = perm_inv[pnod[i].ofSide[j]];
		}
	}
	free(adj_row);
	free(adj);
	free(perm);
	// free(perm_inv);

	return perm_inv;
}
vector<int> GridToRCM_side(int node_num, int *adj_row, int iftera, unsigned int nele, elem *pele, side *psid)
{
	adj_row[0] = 1;
	int adj_num = 0;
	std::vector<int> adjvec;

	for (int i = 0; i < node_num; i++)
	{
		int count = 1;
		std::vector<int> temp;
		temp.insert(temp.end(), i + 1);
		for (int j = 0; j < 2; j++)
		{
			int elemflag = 0;
			int tp1 = psid[i].Node[j];
			for (int k = 0; k < node_num; k++)
			{
				if (i == k)
				{
					continue;
				}
				for (unsigned int l = 0; l < 2; l++)
				{
					if ((tp1 == psid[k].Node[l]))
					{
						elemflag++;
					}
				}
				if (elemflag == 1)
				{
					count++;
					temp.insert(temp.end(), k + 1);
					break;
				}
				elemflag = 0;
			}
		}

		//元素升序排序
		sort(temp.begin(), temp.end());

		// printf("%d  ",i+1);
		for (int ii = 0; ii < count; ii++)
		{
			adjvec.insert(adjvec.end(), temp[ii]);
			// printf("%d  ",temp[ii]);
		}
		// printf("\n");
		adj_row[i + 1] = adj_row[i] + count;
		std::vector<int> tmp1 = temp; // free the vector
		temp.swap(tmp1);
	}
	// adj_num = adjvec.size();
	return adjvec;
}

int *SideReorder(unsigned int nsid, unsigned int nele, unsigned int nnod,
		 side *psid, elem *pele, node *pnod)
{ //调整排序，按照ofElem[0] >< 0界限以及大小顺序排序

	std::vector<iele> ielevec;
	ielevec.reserve(nsid);
	int *sideindex = (int *)calloc(nsid, sizeof(int));
	int pbsdnum = 0;
	unsigned int iftera;
	if (pele[0].Node[3] == pele[0].Node[0])
	{ //三角形
		iftera = 3;
	}
	else
	{ //四边形
		iftera = 4;
	}

	//-------------- The sorting algorithm for SIDE data accroding to ELEM -----------//

	for (unsigned int i = 0; i < nsid; i++)
	{
		int maxvalue = 0;
		for (unsigned int j = 0; j < 2; j++)
		{
			ielevec[i].Node[j] = psid[i].Node[j];
			ielevec[i].ofElem[j] = psid[i].ofElem[j];
			ielevec[i].lsdidx[j] = psid[i].lsdidx[j];
			maxvalue = max(maxvalue, psid[i].ofElem[j]);
		}

		ielevec[i].number = maxvalue;
		ielevec[i].idx = i;
	}

	// reorder according to the psid.ofelem >< 0
	/*for(unsigned int i=0; i<nsid-1; i++){
		if( ielevec[i].ofElem[0] == -1 ){
			pbsdnum++;
			continue;
		}
		for(unsigned int j=i+1; j<nsid; j++){
			if( ielevec[j].ofElem[0] == -1 ){
				pbsdnum++;
				swap( ielevec[i] , ielevec[j] );
				break;
			}
		}
	}
	for(unsigned int i=1;i<pbsdnum;i++){
		if( ielevec[i].Node[0] == 0 ){
			if( ielevec[i].Node[1] < ielevec[0].Node[1] ){
				swap( ielevec[0] , ielevec[i] );
			}
		}
	}
	for(unsigned int i=0; i<pbsdnum-1; i++){
		for(unsigned int j=i+1; j<pbsdnum; j++){
		    if( ielevec[i].Node[1] == ielevec[j].Node[0] ){
				if( i+1 != j )  { swap( ielevec[i+1] , ielevec[j] );}
				break;
			}
		}
	}*/

	//  min-Node renumber
	for (unsigned int i = pbsdnum; i < nsid - 1; i++)
	{
		int minnumber = ielevec[i].Node[0]; // min( ielevec[i].Node[0] , ielevec[i].Node[1] );
		for (unsigned int j = i + 1; j < nsid; j++)
		{
			int tempmin = ielevec[j].Node[0]; // min( ielevec[j].Node[0] , ielevec[j].Node[1] );
			if (minnumber > tempmin)
			{
				swap(ielevec[i], ielevec[j]);
				minnumber = tempmin;
			}
		}
	}
	for (unsigned int i = pbsdnum; i < nsid - 1; i++)
	{
		int minnumber = ielevec[i].Node[1]; // min( ielevec[i].Node[0] , ielevec[i].Node[1] );
		for (unsigned int j = i + 1; j < nsid; j++)
		{
			if (ielevec[i].Node[0] != ielevec[j].Node[0])
			{
				break;
			};
			int tempmin = ielevec[j].Node[1]; // min( ielevec[j].Node[0] , ielevec[j].Node[1] );
			if (minnumber > tempmin)
			{
				swap(ielevec[i], ielevec[j]);
				minnumber = tempmin;
			}
		}
	}

	// write back to psid
	for (unsigned int i = 0; i < nsid; i++)
	{
		for (unsigned int j = 0; j < 2; j++)
		{
			psid[i].Node[j] = ielevec[i].Node[j];
			psid[i].ofElem[j] = ielevec[i].ofElem[j];
			psid[i].lsdidx[j] = ielevec[i].lsdidx[j];
		}
		sideindex[ielevec[i].idx] = i;
	}
	//-------------- The sorting algorithm for SIDE data accroding to ELEM -----------//

	for (unsigned int i = 0; i < nele; i++)
	{ //重写pele中面相邻边side的数据
		for (unsigned int j = 0; j < 4; j++)
		{
			pele[i].Side[j] = sideindex[pele[i].Side[j]];
		}
	}
	for (unsigned int i = 0; i < nnod; i++)
	{ //重写pnod中点连接边ofside的数据
		for (unsigned int j = 0; j < pnod[i].nNElem; j++)
		{
			pnod[i].ofSide[j] = sideindex[pnod[i].ofSide[j]];
		}
	}
	std::vector<iele> tmp = ielevec; // free the vector
	ielevec.swap(tmp);

	return sideindex;
}

void Elemordernode(unsigned int nele, unsigned int nnod, unsigned int nsid, elem *pele, side *psid, int *elemind)
{

	std::vector<esort> evec;
	evec.reserve(nele);
	unsigned int iftera;
	if (pele[0].Node[3] == pele[0].Node[0])
	{ //三角形
		iftera = 3;
	}
	else
	{ //四边形
		iftera = 4;
	}
	for (unsigned int i = 0; i < nele; i++)
	{
		for (unsigned int j = 0; j < 4; j++)
		{
			evec[i].Node[j] = pele[i].Node[j];
		}
		evec[i].volr = pele[i].volr;
		evec[i].xyc[0] = pele[i].xyc[0];
		evec[i].xyc[1] = pele[i].xyc[1];
	}

	for (unsigned int i = 0; i < nele; i++)
	{
		int minvalue = 999999;
		int minindex = -1;
		for (unsigned int j = 0; j < iftera; j++)
		{
			minvalue = min(minvalue, pele[i].Node[j]);
			if (minvalue == pele[i].Node[j])
			{
				minindex = j;
			}
		}
		if (iftera == 3)
		{
			switch (minindex)
			{
			case 1:
			{
				evec[i].Node[0] = pele[i].Node[1];
				evec[i].Node[1] = pele[i].Node[2];
				evec[i].Node[2] = pele[i].Node[0];
				evec[i].Node[3] = pele[i].Node[1];
				break;
			}
			case 2:
			{
				evec[i].Node[0] = pele[i].Node[2];
				evec[i].Node[1] = pele[i].Node[0];
				evec[i].Node[2] = pele[i].Node[1];
				evec[i].Node[3] = pele[i].Node[2];
			}
			}
		}
		else
		{
			switch (minindex)
			{
			case 1:
			{
				evec[i].Node[0] = pele[i].Node[1];
				evec[i].Node[1] = pele[i].Node[2];
				evec[i].Node[2] = pele[i].Node[3];
				evec[i].Node[3] = pele[i].Node[0];
				break;
			}
			case 2:
			{
				evec[i].Node[0] = pele[i].Node[2];
				evec[i].Node[1] = pele[i].Node[3];
				evec[i].Node[2] = pele[i].Node[0];
				evec[i].Node[3] = pele[i].Node[1];
				break;
			}
			case 3:
			{
				evec[i].Node[0] = pele[i].Node[3];
				evec[i].Node[1] = pele[i].Node[0];
				evec[i].Node[2] = pele[i].Node[1];
				evec[i].Node[3] = pele[i].Node[2];
			}
			}
		}
	}

	//  max-value renumber
	for (unsigned int i = 0; i < nele; i++)
	{
		int maxvalue = 0;
		for (unsigned int j = 0; j < 4; j++)
		{
			maxvalue = max(maxvalue, pele[i].Node[j]);
		}
		evec[i].number = maxvalue;
		evec[i].idx = i;
	}

	for (unsigned int i = 0; i < nele - 1; i++)
	{
		for (unsigned int j = i + 1; j < nele; j++)
		{
			if (evec[i].number > evec[j].number)
			{
				swap(evec[i], evec[j]);
			}
		}
	}
	for (unsigned int i = 0; i < nele; i++)
	{
		for (unsigned int j = 0; j < 4; j++)
		{
			pele[i].Node[j] = evec[i].Node[j];
		}
	}

	//  min-value renumber
	for (unsigned int i = 0; i < nele; i++)
	{
		int minvalue = 999999;
		for (unsigned int j = 0; j < 4; j++)
		{
			minvalue = min(minvalue, pele[i].Node[j]);
		}
		evec[i].number = minvalue;
	}
	for (unsigned int i = 0; i < nele - 1; i++)
	{
		for (unsigned int j = i + 1; j < nele; j++)
		{
			if (evec[i].number > evec[j].number)
			{
				swap(evec[i], evec[j]);
			}
		}
	}

	for (unsigned int i = 0; i < nele; i++)
	{
		for (unsigned int j = 0; j < 4; j++)
		{
			pele[i].Node[j] = evec[i].Node[j];
		}
		pele[i].volr = evec[i].volr;
		pele[i].xyc[0] = evec[i].xyc[0];
		pele[i].xyc[1] = evec[i].xyc[1];
		elemind[evec[i].idx] = i;
		// printf("%d  %d\n",i,evec[i].idx);
	}
	/*for(unsigned int i=0;i<nsid;i++){
		psid[i].ofElem[1] = elemind[ psid[i].ofElem[1] ];
		if( psid[i].ofElem[0] != -1 ){
		    psid[i].ofElem[0] = elemind[ psid[i].ofElem[0] ];
		}
	}*/
	std::vector<esort> tmp = evec; // free the vector
	evec.swap(tmp);

	return;
}

// Coloring algorithm according to the edge
int **colorEdge(unsigned int nele, unsigned int nsid, unsigned int colorNum, elem *pele, side *psid)
{

	colorNum = 1;
	int *colorIdx = (int *)calloc(nsid, sizeof(int));
	for (unsigned int i = 0; i < nsid; i++)
	{
		// unsigned int colorItem = 1;
		bool colorflag = false;
		int elem0 = psid[i].ofElem[0];
		int elem1 = psid[i].ofElem[1];
		for (unsigned int colorItem = 0; colorItem < colorNum; colorItem++)
		{
			bool faceflag = true;
			for (unsigned int k = 0; k < i; k++)
			{
				if (colorItem != colorIdx[k])
				{
					continue;
				}
				if ((elem0 < 0) && (psid[k].ofElem[0] >= 0))
				{
					if ((elem1 == psid[k].ofElem[0]) || (elem1 == psid[k].ofElem[1]))
					{
						faceflag = false;
						break;
					}
				}
				else if ((elem0 >= 0) && (psid[k].ofElem[0] < 0))
				{
					if ((elem0 == psid[k].ofElem[1]) || (elem1 == psid[k].ofElem[1]))
					{
						faceflag = false;
						break;
					}
				}
				else if ((elem0 < 0) && (psid[k].ofElem[0] < 0))
				{
					if (elem1 == psid[k].ofElem[1])
					{
						faceflag = false;
						break;
					}
				}
				else
				{
					if ((elem0 == psid[k].ofElem[0]) || (elem0 == psid[k].ofElem[1]) || (elem1 == psid[k].ofElem[0]) || (elem1 == psid[k].ofElem[1]))
					{
						faceflag = false;
						break;
					}
				}
			}
			if (faceflag)
			{
				colorIdx[i] = colorItem;
				colorflag = true;
				break;
			}
		}
		if (!colorflag)
		{
			colorIdx[i] = colorNum;
			colorNum++;
		}
		// printf( "%d  %d  %d  %d\n",i+1,elem0,elem1,colorIdx[i] );
	}

	int *colortemp = (int *)calloc(colorNum, sizeof(int));
	for (unsigned int i = 0; i < nsid; i++)
	{
		colortemp[colorIdx[i]]++;
	}
	int **coloredSide = (int **)malloc(colorNum * sizeof(int *));
	for (unsigned int i = 0; i < colorNum; i++)
	{
		coloredSide[i] = (int *)calloc(colortemp[i] + 1, sizeof(int));
		coloredSide[i][0] = colortemp[i];
		colortemp[i] = 0;
	}
	for (unsigned int i = 0; i < nsid; i++)
	{
		coloredSide[colorIdx[i]][colortemp[colorIdx[i]] + 1] = i;
		colortemp[colorIdx[i]]++;
	}
	/*for(unsigned int i=0; i<colorNum; i++){
		for(unsigned int j=0; j<coloredSide[i][0]; j++){
			printf("%d  %d  %d\n",i+1,j+1,coloredSide[i][j+1]);
		}
	}*/
	free(colorIdx);
	free(colortemp);
	return coloredSide;
}