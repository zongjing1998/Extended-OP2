#include "def2d.h"
#include <vector>
#include <map>
using namespace std;

int adj_bandwidth(int node_num, int adj_num, int adj_row[], int adj[]);
bool adj_contains_ij(int node_num, int adj_num, int adj_row[], int adj[],
	int i, int j);
void adj_insert_ij(int node_num, int adj_max, int *adj_num, int adj_row[],
	int adj[], int i, int j);
int adj_perm_bandwidth(int node_num, int adj_num, int adj_row[], int adj[],
	int perm[], int perm_inv[]);
void adj_perm_show(int node_num, int adj_num, int adj_row[], int adj[],
	int perm[], int perm_inv[]);
void adj_show(int node_num, int adj_num, int adj_row[], int adj[]);
void degree(int root, int adj_num, int adj_row[], int adj[], int mask[],
	int deg[], int *iccsze, int ls[], int node_num);
void genrcm(int node_num, int adj_num, int adj_row[], int adj[], int perm[]);          
void level_set(int root, int adj_num, int adj_row[], int adj[], int mask[],
	int *level_num, int level_row[], int level[], int node_num);
void perm_inverse3(int n, int perm[], int perm_inv[]);
void rcm(int root, int adj_num, int adj_row[], int adj[], int mask[],
	int perm[], int *iccsze, int node_num);
void root_find(int *root, int adj_num, int adj_row[], int adj[], int mask[],
	int *level_num, int level_row[], int level[], int node_num);
int i4_max(int i1, int i2);

void GridRenumber_node(unsigned int nele, unsigned int nnod, elem *pele, node *pnod);
vector<int> GridToRCM_node( int node_num, int *adj_row, int iftera, unsigned int nele, elem *pele );

int* ElemReorder( unsigned int nele, unsigned int nsid, elem *pele, side *psid );
int* SideReorder( unsigned int nsid, unsigned int nele, unsigned int nnod, side *psid, elem *pele, node *pnod );
void Elemordernode( unsigned int nele , unsigned int nnod, unsigned int nsid, elem *pele, side *psid, int *elemind );

int* GridRenumber_elem(unsigned int nele, unsigned int nnod, unsigned int nsid, elem *pele, node *pnod, side *psid);
vector<int> GridToRCM_elem( int node_num, int *adj_row, int iftera, unsigned int nele, elem *pele );

int* GridRenumber_side(unsigned int nele, unsigned int nnod, unsigned int nsid, elem *pele, node *pnod, side *psid);
vector<int> GridToRCM_side( int node_num, int *adj_row, int iftera, unsigned int nele, elem *pele, side *psid );

int** colorEdge( unsigned int nele, unsigned int nsid, unsigned int colorNum, elem *pele, side *psid );